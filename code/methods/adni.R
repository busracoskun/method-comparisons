# ADNI integration-method wrappers, evaluation metrics, and the bootstrap driver
# Sourced into the dataset environment by load_dataset() (code/load_dataset.R)

# RESULT SCHEMA AND MULTIPLE-IMPUTATION SUMMARIES

metric_columns <- c(
  "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag", "rand_miss", "ari_miss", "chi_miss",
  "db_miss", "pve_miss", "R2"
)

mean_or_na <- function(x) {
  if (!length(x) || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

# Average metrics over imputations; a record is "partial" if only some imputations fit
summarise_mi_results <- function(raw) {
  out <- raw %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(metric_columns), mean_or_na),
      model_file = {
        files <- unique(model_file[!is.na(model_file) & nzchar(model_file)])
        if (length(files)) paste(files, collapse = ";") else NA_character_
      },
      n_attempted = sum(n_attempted, na.rm = TRUE),
      n_successful = sum(n_successful, na.rm = TRUE),
      n_failed = sum(n_failed, na.rm = TRUE),
      error_message = {
        messages <- unique(error_message[!is.na(error_message) & nzchar(error_message)])
        if (length(messages)) paste(messages, collapse = " | ") else NA_character_
      },
      elapsed_seconds = sum(elapsed_seconds, na.rm = TRUE),
      .groups = "drop"
    )

  # Classify each record by how many imputations were fitted
  out %>%
    dplyr::mutate(
      status = dplyr::case_when(
        n_failed == 0L ~ "success", n_successful > 0L ~ "partial", TRUE ~ "failed"
      )
    )
}

# BOOTSTRAP DRIVER

# Replicate b resamples patients with set.seed(seed + b), so replicates are reproducible
# individually and independent of mc.cores. Resampled rows get unique "<id>_b<b>_<i>" names
bootstrap_pipe_parallel <- function(
  X_stand, layer_vec, diag_labels, miss_labels, mmse_vec = NULL,
  style = c("CC", "CC_SF", "SI", "SI_SF", "MI", "MI_SF", "NA"),
  B = 30, M = 5, models, n_factors = 3, seed = 100, boot_start = 1L, mc.cores = 20,
  imputed_list = NULL, save_path = NULL, save_format = "both", save_models = FALSE,
  model_dir = CLEAN_MODEL_DIR, verbose = TRUE
) {
  # Validate options and prepare a single imputation when an SI strategy needs one
  style <- match.arg(style)
  save_format <- match.arg(save_format)
  ids_all <- rownames(X_stand)
  if (style %in% c("SI", "SI_SF") && is.null(imputed_list)) {
    if (verbose) cat("Computing single imputation before bootstrap...\n")
    set.seed(seed)
    X_si_full <- missForest::missForest(X_stand)$ximp
  } else if (style %in% c("SI", "SI_SF")) {
    X_si_full <- imputed_list[[1]]
  }

  # Fit all requested models on one bootstrap replicate
  run_one_boot <- function(b) {
    # Resample patients and carry labels along with the resampled rows
    set.seed(seed + b)
    ids_b <- sample(ids_all, length(ids_all), replace = TRUE)
    X_b <- X_stand[ids_b, , drop = FALSE]
    rownames(X_b) <- paste0(ids_b, "_b", b, "_", seq_along(ids_b))
    diag_b <- factor(diag_labels[ids_b])
    names(diag_b) <- rownames(X_b)
    miss_b <- factor(miss_labels[ids_b])
    names(miss_b) <- rownames(X_b)

    mmse_b <- if (!is.null(mmse_vec)) setNames(mmse_vec[ids_b], rownames(X_b)) else NULL
    is_cc <- ids_b %in% common_ids
    fs_diag_b <- setNames(diag_labels[ids_b], rownames(X_b))

    # Feature screening (*_SF) is repeated inside every replicate
    res <- switch(style,
      CC = {
        # Complete cases with fully observed features
        X_cc_all <- X_b[is_cc, , drop = FALSE]
        diag_cc <- diag_b[is_cc]
        miss_cc <- miss_b[is_cc]
        mmse_cc <- if (!is.null(mmse_b)) mmse_b[is_cc] else NULL
        X_cc_all <- X_cc_all[, colSums(is.na(X_cc_all)) == 0, drop = FALSE]
        run_complete(
          X_cc_all, layer_vec, diag_cc, miss_cc, mmse_cc, models, n_factors,
          save_models = save_models, model_dir = model_dir, style = style, boot = b
        )
      },
      CC_SF = {
        # Complete cases followed by SIS screening
        X_cc <- X_b[is_cc, , drop = FALSE]
        X_cc <- X_cc[, colSums(is.na(X_cc)) == 0, drop = FALSE]
        sel <- sel_SIS(X_cc, fs_diag_b[rownames(X_cc)])
        X_sf <- X_cc[, sel, drop = FALSE]

        run_complete(
          X_sf, layer_vec[sel], diag_b[rownames(X_sf)], miss_b[rownames(X_sf)],
          if (!is.null(mmse_b)) mmse_b[rownames(X_sf)] else NULL, models, n_factors,
          save_models = save_models, model_dir = model_dir, style = style, boot = b
        )
      },
      SI = {
        # Bootstrap from the precomputed single imputation
        X_si <- X_si_full[ids_b, , drop = FALSE]
        rownames(X_si) <- rownames(X_b)
        run_complete(
          X_si, layer_vec, diag_b, miss_b, mmse_b, models, n_factors,
          save_models = save_models, model_dir = model_dir, style = style, boot = b
        )
      },
      SI_SF = {
        # Bootstrap from the precomputed single imputation
        X_si <- X_si_full[ids_b, , drop = FALSE]
        rownames(X_si) <- rownames(X_b)
        sel <- sel_SIS(X_si, fs_diag_b[rownames(X_si)])
        X_sf <- X_si[, sel, drop = FALSE]
        run_complete(
          X_sf, layer_vec[sel], diag_b[rownames(X_sf)], miss_b[rownames(X_sf)],
          if (!is.null(mmse_b)) mmse_b[rownames(X_sf)] else NULL, models, n_factors,
          save_models = save_models, model_dir = model_dir, style = style, boot = b
        )
      },
      MI = {
        # Resample every imputed matrix with the same patient draw
        imp_list_b <- lapply(imputed_list, function(Xi) {
          Xi_b <- Xi[ids_b, , drop = FALSE]
          rownames(Xi_b) <- rownames(X_b)
          Xi_b
        })

        raw <- run_multiple_complete(
          imp_list_b, layer_vec,
          fs_diag_labels = fs_diag_b, cluster_diag_labels = diag_b,
          cluster_miss_labels = miss_b, clin_df = data.frame(MMSE = mmse_b), models = models,
          fs = FALSE, n_factors = n_factors, save_models = save_models, model_dir = model_dir,
          style = style, boot = b
        )
        summarise_mi_results(raw)
      },
      MI_SF = {
        # Resample every imputed matrix with the same patient draw
        imp_list_b <- lapply(imputed_list, function(Xi) {
          Xi_b <- Xi[ids_b, , drop = FALSE]
          rownames(Xi_b) <- rownames(X_b)
          Xi_b
        })

        raw <- run_multiple_complete(
          imp_list_b, layer_vec,
          fs_diag_labels = fs_diag_b, cluster_diag_labels = diag_b,
          cluster_miss_labels = miss_b, clin_df = data.frame(MMSE = mmse_b), models = models,
          fs = TRUE, n_factors = n_factors, save_models = save_models, model_dir = model_dir,
          style = style, boot = b
        )
        summarise_mi_results(raw)
      },
      {
        # Incomplete matrix for methods that accept missing values
        run_NA(
          X_b, layer_vec, diag_b, miss_b, mmse_b, models, n_factors,
          save_models = save_models, model_dir = model_dir, style = style, boot = b
        )
      }
    )

    if (!is.data.frame(res)) {
      stop("run_one_boot iteration ", b, " did not return a data frame")
    }
    dplyr::mutate(res, boot = b)
  }

  # Run the replicates in parallel, turning errors into try-error records
  if (verbose) cat("Starting parallel bootstrap with", mc.cores, "cores...\n")
  boot_ids <- seq.int(boot_start, boot_start + B - 1L)
  results_list <- parallel::mclapply(
    boot_ids,
    function(b) {
      tryCatch(
        run_one_boot(b),
        error = function(e) {
          msg <- paste0(
            "BOOT ", b, " FAILED\n",
            "Message: ", conditionMessage(e), "\n",
            paste(utils::capture.output(traceback()), collapse = "\n")
          )
          structure(msg, class = "try-error")
        }
      )
    },
    mc.cores = mc.cores, mc.preschedule = FALSE, mc.silent = FALSE
  )

  # Failed replicates are reported and dropped; they are not re-run or imputed
  errors <- sapply(results_list, function(x) inherits(x, "try-error"))
  if (any(errors)) {
    n_errors <- sum(errors)
    warning(n_errors, " bootstrap iterations failed")
    cat("First error:\n")
    print(results_list[[which(errors)[1]]])
    results_list <- results_list[!errors]
  }

  if (length(results_list) == 0) {
    stop("All bootstrap iterations failed")
  }

  # Combine replicates and optionally save them
  final_results <- dplyr::bind_rows(results_list)
  if (!is.null(save_path)) {
    if (verbose) cat("Saving results to:", save_path, "\n")
    if (save_format %in% c("rds", "both")) {
      rds_path <- sub("\\.(rds|csv)$", "", save_path)
      rds_path <- paste0(rds_path, ".rds")
      saveRDS(final_results, rds_path)
      if (verbose) cat("  Saved RDS:", rds_path, "\n")
    }

    if (save_format %in% c("csv", "both")) {
      csv_path <- sub("\\.(rds|csv)$", "", save_path)
      csv_path <- paste0(csv_path, ".csv")
      write.csv(final_results, csv_path, row.names = FALSE)
      if (verbose) cat("  Saved CSV:", csv_path, "\n")
    }
  }

  if (verbose) cat("Bootstrap completed successfully!\n")
  return(final_results)
}

# BIDIFAC NUMERICAL SAFEGUARDS

# BIDIFAC repeatedly calls softSVD(). LAPACK's dgesdd occasionally fails on bootstrap
# samples, so this softSVD() replaces the vendored one: it keeps the same update but falls
# back to an eigendecomposition when svd() errors. Calls and fallbacks are counted
bidifac_svd_counts <- new.env(parent = emptyenv())

bidifac_svd_counts$total <- 0L
bidifac_svd_counts$fallback <- 0L
bidifac_svd_counts$total_by_context <- integer()
bidifac_svd_counts$fallback_by_context <- integer()

# Record one SVD call, and whether it needed the fallback
add_bidifac_svd_count <- function(context = "unknown", fallback = FALSE) {
  bidifac_svd_counts$total <- bidifac_svd_counts$total + 1L
  old_total <- bidifac_svd_counts$total_by_context[context]
  if (is.na(old_total)) old_total <- 0L
  bidifac_svd_counts$total_by_context[context] <- old_total + 1L

  if (isTRUE(fallback)) {
    bidifac_svd_counts$fallback <- bidifac_svd_counts$fallback + 1L
    old_fallback <- bidifac_svd_counts$fallback_by_context[context]
    if (is.na(old_fallback)) old_fallback <- 0L
    bidifac_svd_counts$fallback_by_context[context] <- old_fallback + 1L
  }
  invisible(NULL)
}

# SVD with an eigendecomposition fallback when LAPACK fails
safe_svd_for_bidifac <- function(X, context = "unknown") {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop("safe_svd_for_bidifac: matrix contains NA/NaN/Inf before SVD.")
  }

  add_bidifac_svd_count(context, fallback = FALSE)
  tryCatch(
    svd(X),
    error = function(e) {
      add_bidifac_svd_count(context, fallback = TRUE)
      message("BIDIFAC svd() failed in ", context, "; using eigen fallback.")

      # Wide matrix: decompose X X^T
      if (nrow(X) <= ncol(X)) {
        eg <- eigen(tcrossprod(X), symmetric = TRUE)
        d <- sqrt(pmax(eg$values, 0))
        keep <- d > 1e-10
        if (!any(keep)) {
          return(list(u = matrix(0, nrow(X), 0), d = numeric(0), v = matrix(0, ncol(X), 0)))
        }

        u <- eg$vectors[, keep, drop = FALSE]
        d <- d[keep]
        v <- t(X) %*% u %*% diag(1 / d, nrow = length(d))
      } else {
        # Tall matrix: decompose X^T X
        eg <- eigen(crossprod(X), symmetric = TRUE)
        d <- sqrt(pmax(eg$values, 0))
        keep <- d > 1e-10
        if (!any(keep)) {
          return(list(u = matrix(0, nrow(X), 0), d = numeric(0), v = matrix(0, ncol(X), 0)))
        }

        v <- eg$vectors[, keep, drop = FALSE]
        d <- d[keep]
        u <- X %*% v %*% diag(1 / d, nrow = length(d))
      }
      list(u = u, d = d, v = v)
    }
  )
}

# Soft-thresholded SVD used by the BIDIFAC updates
softSVD <- function(X, lambda) {
  svdX <- safe_svd_for_bidifac(X, context = "internal_softSVD")
  nuc <- pmax(svdX$d - lambda, 0)
  if (length(nuc) == 0 || all(nuc <= 0)) {
    out <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  } else {
    keep <- nuc > 0
    out <- svdX$u[, keep, drop = FALSE] %*% diag(nuc[keep], nrow = sum(keep)) %*%
      t(svdX$v[, keep, drop = FALSE])
  }
  list(out = out, nuc = sum(nuc))
}

# MODEL WRAPPERS

# Fit one method and score its latent space. Every wrapper returns Z (samples x factors),
# W (features x factors, or NULL), and the fitted object; scaling inside the packages is
# disabled because X is already standardized
run_model <- function(
  model_name, X, layer_vec, y_diag, y_miss = NULL, mmse_vec = NULL, n_factors = 3,
  save_models = FALSE, model_dir = CLEAN_MODEL_DIR, style = "NA", boot = NA_integer_, imp = NULL
) {
  out <- switch(tolower(model_name),
    # MOFA2: one features-by-samples matrix per view; view and group scaling disabled
    mofa = {
      message("Starting MOFA")
      views <- lapply(unique(layer_vec[colnames(X)]), function(lay) {
        cols <- names(layer_vec)[layer_vec == lay]
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(views) <- unique(layer_vec[colnames(X)])

      # Configure and train MOFA without rescaling the views
      m <- create_mofa_from_matrix(views)
      data_opts <- get_default_data_options(m)
      data_opts$scale_views <- FALSE
      data_opts$scale_groups <- FALSE
      data_opts$center_groups <- FALSE
      opts <- get_default_model_options(m)
      opts$num_factors <- n_factors
      train_opts <- get_default_training_options(m)
      m <- prepare_mofa(
        m, data_options = data_opts, model_options = opts, training_options = train_opts
      )
      m <- run_mofa(m, use_basilisk = TRUE)

      # Extract factors and stack the view weights in feature order
      Z <- get_factors(m, factors = "all")$group1
      W_list <- get_weights(m, views = "all", factors = "all")
      W_mat <- do.call(rbind, W_list)
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(Z = Z, W = W_mat, fitted_model = m)
    },
    # iClusterPlus: Gaussian gene and metabolite blocks with fixed lasso penalties (0.03)
    icluster = {
      message("Starting iCluster")
      gene_mat <- as.matrix(X[, layer_vec == "gene", drop = FALSE])
      meta_mat <- as.matrix(X[, layer_vec == "meta", drop = FALSE])
      fit <- iClusterPlus(
        dt1 = gene_mat, dt2 = meta_mat,
        type = c("gaussian", "gaussian"), lambda = c(0.03, 0.03), K = n_factors
      )
      Z <- fit$meanZ
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      rownames(fit$beta[[1]]) <- colnames(gene_mat)
      rownames(fit$beta[[2]]) <- colnames(meta_mat)
      W_mat <- rbind(fit$beta[[1]], fit$beta[[2]])
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(Z = Z, W = W_mat, fitted_model = fit)
    },
    # FABIA: features-by-samples input; centering, normalization, and scaling disabled
    fabia = {
      message("Starting FABIA")

      f <- fabia(as.matrix(t(X)), p = n_factors, nL = n_factors, center = 0, norm = 0, scale = 0.0)

      Z <- t(f@Z)
      rownames(Z) <- colnames(f@Z)
      W_mat <- f@L
      rownames(W_mat) <- colnames(X)
      list(Z = Z, W = W_mat, fitted_model = f)
    },
    # SLIDE requires complete data. Following the package example, missing cells are
    # mean-imputed, columns centered, each block scaled to unit Frobenius norm, and the
    # structure selected by bi-cross-validation before the final fit
    slide = {
      message("Starting SLIDE")

      # Tuning settings, overridable from the calling environment
      n_lambda <- get0("slide_n_lambda", ifnotfound = 50L, inherits = TRUE)
      lambda_min <- get0("slide_lambda_min", ifnotfound = 0.01, inherits = TRUE)
      n_fold <- get0("slide_n_fold", ifnotfound = 3L, inherits = TRUE)
      p_fold <- get0("slide_p_fold", ifnotfound = 3L, inherits = TRUE)
      center <- get0("slide_center", ifnotfound = TRUE, inherits = TRUE)
      k_max <- get0("slide_k_max", ifnotfound = 5000L, inherits = TRUE)
      eps <- get0("slide_eps", ifnotfound = 1e-6, inherits = TRUE)
      ratio_max <- get0("slide_ratio_max", ifnotfound = NULL, inherits = TRUE)

      # Arrange columns by block and mean-impute missing cells
      blocks <- unique(layer_vec[colnames(X)])
      cols_list <- lapply(blocks, function(b) {
        colnames(X)[layer_vec[colnames(X)] == b]
      })
      pvec <- as.integer(vapply(cols_list, length, integer(1)))
      X_slide <- do.call(cbind, lapply(cols_list, function(cols) {
        as.matrix(X[, cols, drop = FALSE])
      }))
      for (j in seq_len(ncol(X_slide))) {
        if (anyNA(X_slide[, j])) {
          mu <- mean(X_slide[, j], na.rm = TRUE)
          if (!is.finite(mu)) mu <- 0
          X_slide[is.na(X_slide[, j]), j] <- mu
        }
      }
      # Center columns and normalize each modality block by its Frobenius norm
      X_std <- X_slide
      d <- length(pvec)
      norms <- rep(0, d)
      svec <- rep(0, d)
      pcum <- c(0L, cumsum(pvec))
      if (center) {
        Xmean <- colMeans(X_std)
        X_std <- X_std - matrix(Xmean, nrow(X_std), ncol(X_std), byrow = TRUE)
      } else {
        Xmean <- rep(0, ncol(X_std))
      }
      for (i in seq_len(d)) {
        idx <- (pcum[i] + 1L):pcum[i + 1L]
        norms[i] <- sum(X_std[, idx, drop = FALSE]^2)
        if (!is.finite(norms[i]) || norms[i] <= 0) {
          stop("SLIDE: at least one block has zero or non-finite norm.")
        }
        X_std[, idx] <- X_std[, idx, drop = FALSE] / sqrt(norms[i])
        svec[i] <- max(svd(X_std[, idx, drop = FALSE], nu = 1)$d)
      }

      # Select the structure by bi-cross-validation and fit it
      if (is.null(ratio_max)) {
        ratio_max <- min((n_fold - 1) / n_fold, (p_fold - 1) / p_fold)
      }
      out_struct <- SLIDE::create_structure_list(
        X = X_std, pvec = pvec, lambda_max = max(svec), standardized = TRUE,
        lambda_min = lambda_min, n_lambda = n_lambda, ratio_max = ratio_max
      )
      out_bcv <- SLIDE::slide_BCV(
        X = X_std, pvec = pvec, structure_list = out_struct$Slist, n_fold = n_fold, p_fold = p_fold,
        k_max = k_max, eps = eps, center = center
      )
      fit <- SLIDE::slide_givenS(
        X = X_slide, pvec = pvec, S = out_bcv$structure_min, k_max = k_max, eps = eps,
        standardized = TRUE
      )

      # Keep the first n_factors scores and loadings
      U <- as.matrix(fit$U)
      V <- as.matrix(fit$V)
      r <- min(n_factors, ncol(U))
      if (r < 1L) stop("SLIDE: no latent factors were estimated.")
      Z <- U[, seq_len(r), drop = FALSE]
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      W_mat <- V[, seq_len(r), drop = FALSE]
      rownames(W_mat) <- colnames(X_slide)
      colnames(W_mat) <- colnames(Z)
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(
        Z = Z, W = W_mat,
        fitted_model = list(
          standardized_X = X_std, structure_list = out_struct$Slist,
          selected_structure = out_bcv$structure_min, model = fit
        )
      )
    },
    # JIVE: given-rank mode with joint and individual ranks equal to n_factors; the joint
    # scores and loadings are returned
    jive = {
      message("Starting JIVE")
      layers <- unique(layer_vec[colnames(X)])
      views <- lapply(layers, function(lay) {
        cols <- which(layer_vec[colnames(X)] == lay)
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(views) <- layers
      jv <- r.jive::jive(
        data = views, method = "given", rankJ = n_factors, rankA = rep(n_factors, length(views)),
        center = FALSE, scale = FALSE
      )

      # Joint scores and loadings
      pred <- r.jive::jive.predict(data.new = views, jive.output = jv)
      S_mat <- as.matrix(pred$joint.scores)
      L_mat <- as.matrix(pred$joint.load)
      Z <- t(S_mat)
      W_mat <- L_mat
      rownames(Z) <- rownames(X)
      rownames(W_mat) <- colnames(X)
      list(Z = Z, W = W_mat, fitted_model = jv)
    },
    # MCIA: features-by-samples blocks; synthetic variables are the shared scores
    mcia = {
      message("Starting MCIA")
      views <- lapply(unique(layer_vec[colnames(X)]), function(lay) {
        cols <- which(layer_vec[colnames(X)] == lay)
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(views) <- unique(layer_vec[colnames(X)])
      mcia_res <- omicade4::mcia(
        df.list = views, cia.nf = n_factors, cia.scan = FALSE, nsc = TRUE, svd = TRUE
      )

      Z <- as.matrix(mcia_res$mcoa$SynVar)
      rownames(Z) <- rownames(X)
      W_mat <- as.matrix(mcia_res$mcoa$axis)
      rownames(W_mat) <- colnames(X)
      list(Z = Z, W = W_mat, fitted_model = mcia_res)
    },
    # RGCCA: fully connected design, n_factors components per block, scaling disabled. Block
    # scores are concatenated and block loadings zero-padded to all features
    rgcca = {
      message("Starting RGCCA")
      blocks <- lapply(unique(layer_vec), function(lay) {
        as.matrix(X[, layer_vec == lay, drop = FALSE])
      })
      names(blocks) <- unique(layer_vec)
      J <- length(blocks)
      Cmat <- 1 - diag(J)
      fit <- RGCCA::rgcca(
        blocks = blocks, connection = Cmat, ncomp = rep(n_factors, J), scale = FALSE,
        scale_block = FALSE, quiet = TRUE, verbose = FALSE
      )

      # Concatenate block scores and zero-pad block loadings to all features
      Z_list <- fit$Y
      Z <- do.call(cbind, Z_list)
      rownames(Z) <- rownames(X)
      all_feats <- colnames(X)
      W_list <- vector("list", length(fit$a))
      for (i in seq_along(fit$a)) {
        a_i <- fit$a[[i]]
        rgcca_features <- rownames(a_i)
        W_i <- matrix(
          0, nrow = length(all_feats), ncol = ncol(a_i), dimnames = list(all_feats, colnames(a_i))
        )
        W_i[rgcca_features, ] <- a_i
        W_list[[i]] <- W_i
      }
      W_mat <- do.call(cbind, W_list)
      list(Z = Z, W = W_mat, fitted_model = fit)
    },
    # IntNMF needs nonnegative data, so the standardized X is back-transformed to the raw
    # scale of X_all. Block-specific H matrices are stacked as loadings
    intnmf = {
      message("Starting IntNMF")

      # Back-transform to the raw scale
      ref <- X_all[, colnames(X), drop = FALSE]
      mu <- colMeans(ref, na.rm = TRUE)
      sigma <- apply(ref, 2, sd, na.rm = TRUE)
      X_raw <- as.matrix(X)
      for (j in seq_len(ncol(X_raw))) {
        X_raw[, j] <- X_raw[, j] * sigma[j] + mu[j]
      }

      # Split into modality blocks
      input_features <- colnames(X)
      if (
        is.null(input_features) || anyNA(input_features) || any(!nzchar(input_features)) ||
          anyDuplicated(input_features)
      ) {
        stop("IntNMF: input feature names must be present and unique")
      }
      feature_layers <- layer_vec[input_features]
      if (anyNA(feature_layers)) {
        stop("IntNMF: every input feature must have a modality label")
      }
      blocks <- unique(feature_layers)
      dat_list <- lapply(blocks, function(b) {
        cols <- input_features[feature_layers == b]
        as.matrix(X_raw[, cols, drop = FALSE])
      })
      names(dat_list) <- blocks

      fit <- nmf.mnnals(dat = dat_list, k = n_factors)
      Z <- fit$W
      rownames(Z) <- rownames(X_raw)

      if (length(fit$H) != length(dat_list)) {
        stop("IntNMF: fitted loading blocks do not match the input modality blocks")
      }
      W_list <- lapply(seq_along(dat_list), function(i) {
        H_i <- as.matrix(fit$H[[i]])
        block_features <- colnames(dat_list[[i]])
        if (ncol(H_i) != length(block_features)) {
          stop("IntNMF: loading count does not match input features for modality ", blocks[[i]])
        }

        # IntNMF returns H as factors by features. Use retained feature names only when
        # they identify the input block completely; otherwise the input matrix defines order.
        loading_features <- colnames(H_i)
        names_reliable <- !is.null(loading_features) &&
          !anyNA(loading_features) && all(nzchar(loading_features)) &&
          !anyDuplicated(loading_features) && setequal(loading_features, block_features)
        W_i <- t(H_i)
        rownames(W_i) <- if (names_reliable) loading_features else block_features
        W_i
      })
      W_mat <- do.call(rbind, W_list)
      loading_features <- rownames(W_mat)
      if (
        nrow(W_mat) != length(input_features) || is.null(loading_features) ||
          anyNA(loading_features) || any(!nzchar(loading_features)) ||
          anyDuplicated(loading_features) || !setequal(loading_features, input_features)
      ) {
        stop("IntNMF: loading feature names must be unique and complete")
      }
      W_mat <- W_mat[input_features, , drop = FALSE]
      if (!identical(rownames(W_mat), input_features)) {
        stop("IntNMF: failed to align loading rows to the input feature order")
      }
      list(Z = Z, W = W_mat, fitted_model = fit)
    },
    # BIDIFAC: each omics layer is one column block and all samples form one row block
    # Missing modality values stay NA and are handled by impute.BIDIFAC(). Bootstrap failures
    # arose inside softSVD() during fitting (see the safeguards above), not in score extraction
    # rmt = FALSE because X is already centered and scaled; re-estimating the random-matrix
    # scaling in every replicate made fitting more fragile
    bidifac = {
      message("Starting BIDIFAC")

      # One column block per modality, one row block for all samples
      bl <- unique(layer_vec[colnames(X)])
      cols_list <- lapply(bl, function(b) {
        colnames(X)[layer_vec[colnames(X)] == b]
      })
      names(cols_list) <- bl
      feat_order <- unlist(cols_list, use.names = FALSE)
      dat <- matrix(
        vector("list", length(bl)), nrow = 1, ncol = length(bl), dimnames = list("all_samples", bl)
      )
      for (j in seq_along(bl)) {
        cols <- cols_list[[j]]
        M <- as.matrix(X[, cols, drop = FALSE])
        storage.mode(M) <- "double"
        rownames(M) <- rownames(X)
        colnames(M) <- cols
        if (all(is.na(M))) {
          stop("BIDIFAC: modality block ", bl[j], " is entirely NA in this bootstrap sample.")
        }
        if (any(is.nan(M)) || any(is.infinite(M))) {
          stop("BIDIFAC: modality block ", bl[j], " contains NaN/Inf.")
        }
        dat[1, j] <- list(M)
      }

      # Fit, imputing missing modality values when present
      sigma_one <- matrix(1, nrow = 1, ncol = length(bl))
      fit <- if (anyNA(X)) {
        impute.BIDIFAC(dat, rmt = FALSE, sigma = sigma_one, pbar = FALSE)
      } else {
        BIDIFAC(dat, rmt = FALSE, sigma = sigma_one, pbar = FALSE)
      }

      # Latent scores from the SVD of the fitted signal matrix
      Smat <- data.rearrange(fit$S, rmt = FALSE)$out
      Smat <- as.matrix(Smat)
      storage.mode(Smat) <- "double"
      rownames(Smat) <- rownames(X)
      colnames(Smat) <- feat_order
      Smat <- Smat[rownames(X), colnames(X), drop = FALSE]
      if (any(!is.finite(Smat))) {
        stop("BIDIFAC: Smat contains NA/NaN/Inf after fitting.")
      }
      sv <- svd(Smat)
      k <- min(n_factors, ncol(sv$u), ncol(sv$v), length(sv$d))
      if (k < 1L) {
        stop("BIDIFAC: no usable latent factor returned.")
      }
      Z <- sv$u[, seq_len(k), drop = FALSE] %*% diag(sv$d[seq_len(k)], nrow = k)
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(k))
      W_mat <- sv$v[, seq_len(k), drop = FALSE]
      rownames(W_mat) <- colnames(X)
      colnames(W_mat) <- colnames(Z)
      list(Z = Z, W = W_mat, fitted_model = fit)
    },
    # GFA: samples-by-features blocks with the package default options
    gfa = {
      layers <- unique(layer_vec[colnames(X)])
      Y <- lapply(layers, function(l) {
        as.matrix(X[, names(layer_vec)[layer_vec == l], drop = FALSE])
      })
      names(Y) <- layers

      fit <- GFA::gfa(Y = Y, K = n_factors, opts = GFA::getDefaultOpts())
      Z <- as.matrix(fit$X)
      rownames(Z) <- rownames(X)
      W <- as.matrix(fit$W)
      if (is.null(rownames(W))) rownames(W) <- colnames(X)
      W <- W[colnames(X), , drop = FALSE]
      list(Z = Z, W = W, fitted_model = fit)
    },
    # SGCCA (mixOmics): fully connected design, all variables kept in every component,
    # scaling disabled; block loadings are zero-padded to all features
    sgcca = {
      message("Starting SGCCA (mixOmics)")
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      Xb <- lapply(blocks, function(cols) {
        as.matrix(X[, cols, drop = FALSE])
      })

      J <- length(Xb)
      design <- 1 - diag(J)
      keepX <- lapply(Xb, function(M) rep(ncol(M), n_factors))
      sgcca_args <- list(
        X = Xb, design = design, ncomp = n_factors, keepX = keepX, scale = FALSE,
        all.outputs = FALSE
      )
      if ("scheme" %in% names(formals(mixOmics::wrapper.sgcca))) sgcca_args$scheme <- "horst"

      fit <- do.call(mixOmics::wrapper.sgcca, sgcca_args)
      Z <- do.call(cbind, fit$variates)
      rownames(Z) <- rownames(X)
      allf <- colnames(X)
      W <- do.call(cbind, lapply(seq_along(fit$loadings), function(i) {
        L <- as.matrix(fit$loadings[[i]])
        Wi <- matrix(0, nrow = length(allf), ncol = ncol(L), dimnames = list(allf, colnames(L)))
        Wi[rownames(L), ] <- L
        Wi
      }))
      list(Z = Z, W = W, fitted_model = fit)
    },
    # GPCCA: features-by-samples blocks; missing values are handled by the package
    gpcca = {
      message("Starting GPCCA")
      blocks <- unique(layer_vec[colnames(X)])
      X_list <- lapply(blocks, function(b) {
        cols <- colnames(X)[layer_vec[colnames(X)] == b]
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(X_list) <- blocks

      fit <- GPCCA::GPCCAmodel(
        X.list = X_list, d = n_factors, missVal = any(sapply(X_list, anyNA)),
        plot_RMSE = FALSE, verbose = 0
      )

      # Samples-by-factors scores and loadings in the original feature order
      Z <- as.matrix(t(fit$Z))
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      W <- as.matrix(fit$W)
      feat_order <- unlist(lapply(X_list, rownames), use.names = FALSE)
      if (nrow(W) != length(feat_order)) {
        stop()
      }
      rownames(W) <- feat_order
      colnames(W) <- colnames(Z)
      W <- W[colnames(X), , drop = FALSE]
      list(Z = Z, W = W, fitted_model = fit)
    },
    # Block PLS (mixOmics): the first block is used as the response (indY = 1)
    block_pls = {
      blocks <- unique(layer_vec[colnames(X)])
      Xb <- lapply(blocks, function(b) {
        as.matrix(X[, layer_vec[colnames(X)] == b, drop = FALSE])
      })
      names(Xb) <- blocks

      fit <- mixOmics::block.pls(X = Xb, indY = 1, ncomp = n_factors, scale = FALSE)

      Z <- do.call(cbind, fit$variates)
      rownames(Z) <- rownames(X)
      allf <- colnames(X)
      W <- do.call(cbind, lapply(names(fit$loadings), function(b) {
        L <- as.matrix(fit$loadings[[b]])
        Wb <- matrix(
          0, length(allf), ncol(L), dimnames = list(allf, paste0(b, "_comp", seq_len(ncol(L))))
        )
        Wb[rownames(L), ] <- L
        Wb
      }))
      list(Z = Z, W = W, fitted_model = fit)
    },
    # O2PLS: two complete blocks, joint rank capped by the data dimensions, no orthogonal
    # components; the X-joint scores and loadings are returned
    o2pls = {
      message("Starting O2PLS")
      bl <- unique(layer_vec[colnames(X)])
      if (length(bl) < 2) stop()
      X1 <- as.matrix(X[, layer_vec[colnames(X)] == bl[1], drop = FALSE])
      X2 <- as.matrix(X[, layer_vec[colnames(X)] == bl[2], drop = FALSE])

      k_use <- min(n_factors, nrow(X1) - 1L, ncol(X1), ncol(X2))
      if (k_use < 1L) stop("O2PLS: not enough rank for >=1 joint component")
      fit <- OmicsPLS::o2m(X1, X2, k_use, 0, 0)
      Z <- OmicsPLS::scores(fit, which_part = "Xjoint")
      if (is.null(Z) || ncol(Z) < 1) stop("scores() returned empty")
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      W <- OmicsPLS::loadings(fit, loading_name = "Xjoint")
      W <- as.matrix(W)
      rownames(W) <- colnames(X1)
      colnames(W) <- colnames(Z)
      list(Z = Z, W = W, fitted_model = fit)
    },
    # CIA: co-inertia of per-block PCAs on two complete blocks. CIA has no single loading
    # matrix, so W stacks the variable coordinates of both blocks in feature order. This is an
    # approximation: CIA links two PCA spaces rather than factorizing X
    cia = {
      message("Starting CIA")
      bl <- unique(layer_vec[colnames(X)])
      if (length(bl) < 2) stop("CIA needs >=2 blocks")
      X1 <- as.matrix(X[, layer_vec[colnames(X)] == bl[1], drop = FALSE])
      X2 <- as.matrix(X[, layer_vec[colnames(X)] == bl[2], drop = FALSE])

      # Per-block PCA followed by co-inertia analysis
      d1 <- ade4::dudi.pca(X1, center = FALSE, scale = FALSE, scannf = FALSE, nf = n_factors)
      d2 <- ade4::dudi.pca(X2, center = FALSE, scale = FALSE, scannf = FALSE, nf = n_factors)
      co <- ade4::coinertia(d1, d2, scannf = FALSE, nf = n_factors)

      Z <- as.matrix(co$lX)
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))

      # Variable coordinates of both blocks as approximate loadings
      W1 <- as.matrix(co$co[, seq_len(ncol(Z)), drop = FALSE])
      W2 <- as.matrix(co$li[, seq_len(ncol(Z)), drop = FALSE])
      if (is.null(rownames(W1))) rownames(W1) <- colnames(X1)
      if (is.null(rownames(W2))) rownames(W2) <- colnames(X2)
      colnames(W1) <- colnames(Z)
      colnames(W2) <- colnames(Z)
      W_mat <- matrix(0, nrow = ncol(X), ncol = ncol(Z), dimnames = list(colnames(X), colnames(Z)))
      W_mat[rownames(W1), ] <- W1
      W_mat[rownames(W2), ] <- W2
      list(Z = Z, W = W_mat, fitted_model = co)
    },
    # NEMO builds an integrated sample affinity graph from partially overlapping modalities
    # Latent factors come from a spectral embedding of that graph; there are no feature loadings
    nemo = {
      message("Starting NEMO")
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      omics_list <- lapply(names(blocks), function(b) {
        cols <- blocks[[b]]
        Xb <- as.matrix(X[, cols, drop = FALSE])
        # A sample belongs to a modality if at least one of its values is observed
        observed <- rowSums(!is.na(Xb)) > 0
        X_obs <- Xb[observed, , drop = FALSE]
        if (nrow(X_obs) < 2L) {
          stop(paste0("NEMO: modality ", b, " contains fewer than two observed samples."))
        }
        keep_features <- colSums(!is.na(X_obs)) > 0
        X_obs <- X_obs[, keep_features, drop = FALSE]
        if (ncol(X_obs) < 1L) {
          stop(paste0("NEMO: modality ", b, " contains no usable features."))
        }
        # NEMO cannot take cell-level NAs, so only the remaining within-modality cells are
        # mean-imputed; entirely missing modalities are not imputed
        for (j in seq_len(ncol(X_obs))) {
          if (anyNA(X_obs[, j])) {
            mu <- mean(X_obs[, j], na.rm = TRUE)
            X_obs[is.na(X_obs[, j]), j] <- mu
          }
        }
        t(X_obs)
      })
      names(omics_list) <- names(blocks)

      # Build the integrated affinity graph
      k_requested <- get0("nemo_k", ifnotfound = 10L, inherits = TRUE)
      k_max_available <- min(vapply(omics_list, ncol, integer(1))) - 1L
      if (k_max_available < 1L) {
        stop("NEMO: insufficient observed samples.")
      }
      k_use <- min(as.integer(k_requested), k_max_available)
      affinity <- NEMO::nemo.affinity.graph(omics_list, k = k_use)
      # Spectral embedding of the normalized affinity, skipping the leading eigenvector
      d <- rowSums(affinity)
      d[d < 1e-12] <- 1e-12
      D_inv_sqrt <- diag(1 / sqrt(d))
      L_sym <- D_inv_sqrt %*% affinity %*% D_inv_sqrt
      ev <- eigen(L_sym, symmetric = TRUE)
      k_end <- min(n_factors + 1L, ncol(ev$vectors))
      if (k_end < 2L) {
        stop("NEMO: not enough samples to construct latent factors.")
      }
      Z <- ev$vectors[, 2:k_end, drop = FALSE]
      rownames(Z) <- rownames(affinity)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))

      missing_from_graph <- setdiff(rownames(X), rownames(Z))
      if (length(missing_from_graph) > 0L) {
        stop("NEMO: some samples are absent from every usable omics block.")
      }
      Z <- Z[rownames(X), , drop = FALSE]
      list(Z = Z, W = NULL, fitted_model = list(affinity_graph = affinity, k = k_use))
    },
    # BSF (BSFP): features-by-samples blocks. BSFP can call BIDIFAC internally for
    # initialization and rank selection, which caused BSF bootstrap failures. Supplying the ranks
    # directly (joint = n_factors, individual = 0) and nninit = FALSE removes that dependency and
    # matches the fixed factor count used for every method
    bsf = {
      message("Starting BSF")
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      data_list <- lapply(blocks, function(cols) {
        M <- t(as.matrix(X[, cols, drop = FALSE]))
        storage.mode(M) <- "double"
        if (all(is.na(M))) {
          stop("BSF: one modality block is entirely NA in this bootstrap sample.")
        }
        if (any(is.nan(M)) || any(is.infinite(M))) {
          stop("BSF: modality block contains NaN/Inf.")
        }
        M
      })

      # Fixed ranks: n_factors joint components and no individual components
      ranks_use <- c(n_factors, rep(0L, length(data_list)))
      fit <- BSFP::bsfp(
        data = data_list, Y = NULL, nninit = FALSE, ranks = ranks_use, nsample = 1000,
        progress = FALSE, save_structures = FALSE, save_loadings_scores = FALSE,
        save_predictive_model = FALSE, save_imputations = FALSE, save_init = FALSE,
        save_last_sample = TRUE
      )

      # Scores and loadings from the last posterior sample
      V_last <- as.matrix(fit$last.iter$V.iter[[1, 1]])
      U_last <- fit$last.iter$U.iter
      storage.mode(V_last) <- "double"
      if (any(!is.finite(V_last))) {
        stop("BSF: returned non-finite sample scores.")
      }
      k <- min(n_factors, ncol(V_last))
      if (k < 1L) {
        stop("BSF: no usable factor returned.")
      }
      Z <- V_last[, seq_len(k), drop = FALSE]
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(k))

      all_feats <- colnames(X)
      W_mat <- matrix(
        0, nrow = length(all_feats), ncol = ncol(Z), dimnames = list(all_feats, colnames(Z))
      )
      for (b in seq_along(blocks)) {
        feats_b <- blocks[[b]]
        U_b <- as.matrix(U_last[[b, 1]])
        storage.mode(U_b) <- "double"
        U_b <- U_b[, seq_len(ncol(Z)), drop = FALSE]
        if (nrow(U_b) != length(feats_b)) {
          stop("BSF: loading dimension mismatch for block ", names(blocks)[b])
        }
        W_mat[feats_b, ] <- U_b
      }
      if (any(!is.finite(W_mat))) {
        stop("BSF: returned non-finite loadings.")
      }
      list(Z = Z, W = W_mat, fitted_model = fit)
    },
    # Fit Gaussian GIPCA with joint factors only
    # Exclude the intercept from patient scores and feature loadings
    gipca = {
      message("Starting GIPCA")
      load_gipca_once()
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      data_list <- lapply(blocks, function(cols) {
        as.matrix(X[, cols, drop = FALSE])
      })

      fit <- EPCAJIVEMissbio(
        data = data_list, rankj = n_factors, ranka = rep(0L, length(data_list)),
        D = vapply(data_list, ncol, integer(1)), family = rep("gaussian", length(data_list)),
        tol = 0.1, max.iter = 500, n.size = NULL, lambda = 0
      )

      Z <- as.matrix(fit$U.joint[, -1, drop = FALSE])
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      W_mat <- t(as.matrix(fit$V.joint[-1, , drop = FALSE]))
      rownames(W_mat) <- unlist(blocks, use.names = FALSE)
      colnames(W_mat) <- colnames(Z)
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(Z = Z, W = W_mat, fitted_model = fit)
    },
    # tICA is fit to a block x sample x sample tensor whose slices are within-block sample
    # correlation matrices. Its features are sample positions, so DoTICA has no loadings on the
    # columns of X; W is approximated afterwards by regressing X on the sample scores
    tica = {
      message("Starting tICA")
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      B <- length(blocks)
      n <- nrow(X)
      data_a <- array(0, dim = c(B, n, n))

      # One within-block sample-correlation slice per modality
      for (b in seq_len(B)) {
        cols_b <- blocks[[b]]
        Xb <- as.matrix(X[, cols_b, drop = FALSE])
        # Undefined correlations (e.g. constant profiles) are set to 0
        Cb <- suppressWarnings(stats::cor(t(Xb), use = "everything"))
        Cb[!is.finite(Cb)] <- 0
        diag(Cb) <- 1
        Cb <- (Cb + t(Cb)) / 2
        data_a[b, , ] <- Cb
      }

      dim_vec <- rep(n_factors, B)
      fit <- tensorICA::DoTICA(Data = data_a, dim = dim_vec, method = "FOBI")
      # Orient every projected slice as samples x components before combining
      proj_list <- lapply(fit$TICA.projS, function(M) {
        M <- as.matrix(M)
        if (ncol(M) == n) {
          t(M)
        } else if (nrow(M) == n) {
          M
        } else {
          stop()
        }
      })

      # Keep n_factors scores and regress X on them for approximate loadings
      Z_full <- do.call(cbind, proj_list)
      Z <- Z_full[, seq_len(min(n_factors, ncol(Z_full))), drop = FALSE]
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))

      Zc <- scale(Z, center = TRUE, scale = FALSE)
      Xc <- scale(as.matrix(X), center = TRUE, scale = FALSE)
      Bhat <- qr.solve(Zc, Xc)
      W_mat <- t(Bhat)
      rownames(W_mat) <- colnames(X)
      colnames(W_mat) <- colnames(Z)
      list(Z = Z, W = W_mat, fitted_model = fit)
    }
  )

  # Optionally save the fitted model
  Z <- out$Z
  W <- out$W
  model_file <- NA_character_
  if (save_models) {
    fit_obj <- if ("fitted_model" %in% names(out)) out$fitted_model else NULL
    save_obj <- list(
      model = model_name, style = toupper(style), boot = boot, imp = imp, n_factors = n_factors,
      Z = Z, W = W, fitted_model = fit_obj, rows = rownames(X), cols = colnames(X)
    )
    model_file <- save_boot_model(
      obj = save_obj, model_name = model_name, style = style, boot = boot, imp = imp,
      base_dir = model_dir
    )
  }

  # k-means seeds are tied to the bootstrap index: 100 + boot for diagnosis, 200 + boot for
  # missingness categories
  cl_diag <- score_clusters(
    Z, y_diag, X_original = X, W_loadings = W, seed = 100 + as.integer(boot)
  )

  R2_val <- if (!is.null(mmse_vec)) score_regression(Z, mmse_vec) else NA_real_

  result <- list(
    model = model_name, model_file = model_file, rand_diag = cl_diag["rand"],
    ari_diag = cl_diag["ari"], chi_diag = cl_diag["chi"], db_diag = cl_diag["db"],
    pve_diag = cl_diag["pve"], R2 = R2_val
  )

  if (!is.null(y_miss)) {
    cl_miss <- score_clusters(
      Z, y_miss, X_original = X, W_loadings = W, seed = 200 + as.integer(boot)
    )
    result$rand_miss <- cl_miss["rand"]
    result$ari_miss <- cl_miss["ari"]
    result$chi_miss <- cl_miss["chi"]
    result$db_miss <- cl_miss["db"]
    result$pve_miss <- cl_miss["pve"]
  }
  return(result)
}

# EVALUATION METRICS

# k-means (k = number of label levels, 25 starts) on the latent scores, compared with the
# labels by Rand and adjusted Rand indices, plus internal cluster criteria. PVE is the
# variance of Z %*% t(W) relative to the observed variance of X
score_clusters <- function(
  latent_matrix, true_labels, X_original = NULL, W_loadings = NULL, seed = 100
) {
  if (is.null(rownames(latent_matrix))) {
    rownames(latent_matrix) <- names(true_labels) %||% seq_len(nrow(latent_matrix))
  }

  # Align scores and labels
  ids <- intersect(rownames(latent_matrix), names(true_labels))
  latent_matrix <- latent_matrix[ids, , drop = FALSE]
  true_labels <- true_labels[ids]
  stopifnot(nrow(latent_matrix) == length(true_labels))
  k <- nlevels(droplevels(true_labels))

  # Cluster the latent scores and compare with the labels
  set.seed(seed)
  pred <- kmeans(latent_matrix, centers = k, nstart = 25)$cluster
  rand <- clusterCrit::extCriteria(as.integer(true_labels), pred, "Rand")$rand
  ari <- mclust::adjustedRandIndex(true_labels, pred)
  int <- clusterCrit::intCriteria(
    as.matrix(latent_matrix), pred, c("Calinski_Harabasz", "Davies_Bouldin")
  )

  # Variance of the reconstruction relative to the observed variance
  pve <- NA_real_
  if (!is.null(X_original) && !is.null(W_loadings)) {
    feats <- intersect(colnames(X_original), rownames(W_loadings))
    X_use <- as.matrix(X_original[rownames(latent_matrix), feats, drop = FALSE])
    W_use <- W_loadings[feats, , drop = FALSE]
    Xhat <- latent_matrix %*% t(W_use)
    mu_orig <- colMeans(X_use, na.rm = TRUE)
    Xc <- sweep(X_use, 2, mu_orig, FUN = "-")
    mu_hat <- colMeans(Xhat)
    Xhatc <- sweep(Xhat, 2, mu_hat, FUN = "-")
    var_X <- sum(matrixStats::colVars(Xc, na.rm = TRUE))
    var_Xhat <- sum(matrixStats::colVars(Xhatc, na.rm = TRUE))
    pve <- var_Xhat / var_X
  }

  c(rand = rand, ari = ari, chi = int$calinski_harabasz, db = int$davies_bouldin, pve = pve)
}

# In-sample R^2 of the clinical score regressed on all latent factors
score_regression <- function(latent_matrix, mmse_vec) {
  ids <- intersect(rownames(latent_matrix), names(mmse_vec))
  df <- data.frame(MMSE = mmse_vec[ids], latent_matrix[ids, , drop = FALSE])
  fit <- lm(MMSE ~ ., data = df)
  summary(fit)$r.squared
}

# FEATURE SCREENING

# Sure independence screening for AD versus all other diagnoses (binomial, BIC tuning)
sel_SIS <- function(X_complete, label_vec) {
  keep <- !is.na(label_vec)
  ybin <- as.integer(label_vec[keep] == "AD")
  sis <- SIS::SIS(
    x = as.matrix(X_complete[keep, , drop = FALSE]), y = ybin, family = "binomial", tune = "bic",
    nsis = ncol(X_complete)
  )
  colnames(X_complete)[sis$sis.ix0]
}

# MODEL RUNNER

# Fit each requested method in turn. A failure is recorded as a row with status "failed"
# and NA metrics rather than stopping the remaining methods
run_complete <- function(
  X, layer_vec, diag_labels, miss_labels, mmse_vec, models, n_factors = 3, save_models = FALSE,
  model_dir = CLEAN_MODEL_DIR, style = "NA", boot = NA_integer_, imp = NULL
) {
  dplyr::bind_rows(lapply(models, function(m) {
    started <- proc.time()[["elapsed"]]
    tryCatch(
      {
        # Fit and score one model
        out <- run_model(
          model_name = m, X = X, layer_vec = layer_vec[colnames(X)],
          y_diag = diag_labels[rownames(X)], y_miss = miss_labels[rownames(X)],
          mmse_vec = if (!is.null(mmse_vec)) mmse_vec[rownames(X)] else NULL,
          n_factors = n_factors, save_models = save_models, model_dir = model_dir,
          style = style, boot = boot, imp = imp
        )
        # Fill any metric columns the model did not return
        result <- tibble::as_tibble(out)
        if (!"model_file" %in% names(result)) result$model_file <- NA_character_
        for (metric in metric_columns) {
          if (!metric %in% names(result)) result[[metric]] <- NA_real_
        }

        result %>%
          dplyr::mutate(
            status = "success", error_message = NA_character_,
            n_attempted = 1L, n_successful = 1L, n_failed = 0L,
            elapsed_seconds = proc.time()[["elapsed"]] - started
          )
      },
      error = function(e) {
        message(
          "MODEL FAILED", " | model=", m, " | style=", style, " | boot=", boot,
          if (!is.null(imp)) paste0(" | imp=", imp) else "",
          " | ", conditionMessage(e)
        )
        tibble::tibble(
          model = m, model_file = NA_character_,
          rand_diag = NA_real_, ari_diag = NA_real_, chi_diag = NA_real_, db_diag = NA_real_,
          pve_diag = NA_real_, rand_miss = NA_real_, ari_miss = NA_real_, chi_miss = NA_real_,
          db_miss = NA_real_, pve_miss = NA_real_, R2 = NA_real_,
          status = "failed", error_message = conditionMessage(e),
          n_attempted = 1L, n_successful = 0L, n_failed = 1L,
          elapsed_seconds = proc.time()[["elapsed"]] - started
        )
      }
    )
  }))
}
