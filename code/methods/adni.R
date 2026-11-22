# ADNI statistical integration methods and bootstrap helpers.
# Preserved source expressions; see docs/reproducibility.md.
# Extracted from final.R:74
metric_columns <- c(
  "rand_diag",
  "ari_diag",
  "chi_diag",
  "db_diag",
  "pve_diag",
  "rand_miss",
  "ari_miss",
  "chi_miss",
  "db_miss",
  "pve_miss",
  "R2"
)

# Extracted from final.R:88
mean_or_na <- function(x) {
  if (!length(x) || all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

# Extracted from final.R:96
summarise_mi_results <- function(raw) {

  out <- raw %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(metric_columns), mean_or_na),

      model_file = {
        files <- unique(model_file[ !is.na(model_file) & nzchar(model_file) ])

        if (length(files)) {
          paste(files, collapse = ";")
        } else {
          NA_character_
        }
      },

      n_attempted = sum(n_attempted, na.rm = TRUE),
      n_successful = sum(n_successful, na.rm = TRUE),
      n_failed = sum(n_failed, na.rm = TRUE),

      error_message = {
        messages <- unique(error_message[ !is.na(error_message) & nzchar(error_message) ])

        if (length(messages)) {
          paste(messages, collapse = " | ")
        } else {
          NA_character_
        }
      },

      elapsed_seconds = sum(elapsed_seconds, na.rm = TRUE),

      .groups = "drop"
    )

  out %>%
    dplyr::mutate(
      status = dplyr::case_when(
        n_failed == 0L ~ "success",
        n_successful > 0L ~ "partial",
        TRUE ~ "failed"
      )
    )
}

# Extracted from final.R:164
bootstrap_pipe_parallel <- function(X_stand, layer_vec, diag_labels, miss_labels, mmse_vec = NULL,
                                    style = c("CC","CC_SF","SI","SI_SF","MI","MI_SF","NA"),
                                    B = 30, M = 5, models, n_factors = 3, seed = 100,
                                    boot_start = 1L,
                                    mc.cores = 20, imputed_list = NULL,
                                    save_path = NULL, save_format = "both",
                                    save_models = FALSE,
                                    model_dir = CLEAN_MODEL_DIR,
                                    verbose = TRUE) {

  style <- match.arg(style)
  save_format <- match.arg(save_format)
  ids_all <- rownames(X_stand)

  # Pre-compute single imputation if needed for SI/SI_SF styles
  if (style %in% c("SI", "SI_SF") && is.null(imputed_list)) {
    if (verbose) cat("Computing single imputation before bootstrap...\n")
    set.seed(seed)
    X_si_full <- missForest::missForest(X_stand)$ximp
  } else if (style %in% c("SI", "SI_SF")) {
    X_si_full <- imputed_list[[1]]  # Use pre-computed imputation
  }

  run_one_boot <- function(b) {
    set.seed(seed + b)
    ids_b <- sample(ids_all, length(ids_all), replace = TRUE)
    X_b <- X_stand[ids_b, , drop = FALSE]
    rownames(X_b) <- paste0(ids_b, "_b", b, "_", seq_along(ids_b))

    # labels
    diag_b <- factor(diag_labels[ids_b]); names(diag_b) <- rownames(X_b)
    miss_b <- factor(miss_labels[ids_b]); names(miss_b) <- rownames(X_b)
    mmse_b <- if (!is.null(mmse_vec)) setNames(mmse_vec[ids_b], rownames(X_b)) else NULL
    is_cc <- ids_b %in% common_ids

    # feature selection, use original diagnosis labels
    fs_diag_b <- setNames(diag_labels[ids_b], rownames(X_b))

    res <- switch(style,
                  CC = {
                    X_cc_all <- X_b[is_cc, , drop = FALSE]
                    diag_cc <- diag_b[is_cc]
                    miss_cc <- miss_b[is_cc]
                    mmse_cc <- if (!is.null(mmse_b)) mmse_b[is_cc] else NULL
                    X_cc_all <- X_cc_all[, colSums(is.na(X_cc_all)) == 0, drop = FALSE]
                    run_complete(X_cc_all, layer_vec, diag_cc, miss_cc, mmse_cc, models, n_factors,
                                 save_models = save_models, model_dir = model_dir,
                                 style = style, boot = b)
                  },

                  CC_SF = {
                    X_cc <- X_b[is_cc, , drop = FALSE]
                    X_cc <- X_cc[, colSums(is.na(X_cc)) == 0, drop = FALSE]
                    sel <- sel_SIS(X_cc, fs_diag_b[rownames(X_cc)])
                    X_sf <- X_cc[, sel, drop = FALSE]
                    run_complete(X_sf, layer_vec[sel], diag_b[rownames(X_sf)],
                                 miss_b[rownames(X_sf)], if (!is.null(mmse_b)) mmse_b[rownames(X_sf)] else NULL,
                                 models, n_factors,
                                 save_models = save_models, model_dir = model_dir,
                                 style = style, boot = b)
                  },

                  SI = {
                    # FIX: Bootstrap from pre-computed imputation
                    X_si <- X_si_full[ids_b, , drop = FALSE]
                    rownames(X_si) <- rownames(X_b)
                    run_complete(X_si, layer_vec, diag_b, miss_b, mmse_b, models, n_factors,
                                 save_models = save_models, model_dir = model_dir,
                                 style = style, boot = b)
                  },

                  SI_SF = {
                    # FIX: Bootstrap from pre-computed imputation
                    X_si <- X_si_full[ids_b, , drop = FALSE]
                    rownames(X_si) <- rownames(X_b)
                    sel <- sel_SIS(X_si, fs_diag_b[rownames(X_si)])
                    X_sf <- X_si[, sel, drop = FALSE]
                    run_complete(X_sf, layer_vec[sel], diag_b[rownames(X_sf)],
                                 miss_b[rownames(X_sf)], if (!is.null(mmse_b)) mmse_b[rownames(X_sf)] else NULL,
                                 models, n_factors,
                                 save_models = save_models, model_dir = model_dir,
                                 style = style, boot = b)
                  },

                  MI = {
                    imp_list_b <- lapply(imputed_list, function(Xi) {
                      Xi_b <- Xi[ids_b, , drop = FALSE]
                      rownames(Xi_b) <- rownames(X_b)
                      Xi_b
                    })
                    raw <- run_multiple_complete(imp_list_b, layer_vec, fs_diag_labels = fs_diag_b,
                                                 cluster_diag_labels = diag_b, cluster_miss_labels = miss_b,
                                                 clin_df = data.frame(MMSE = mmse_b), models = models,
                                                 fs = FALSE, n_factors = n_factors,
                                                 save_models = save_models, model_dir = model_dir,
                                                 style = style, boot = b)
                    
                    # metric_summary <- raw %>%
                    #   dplyr::group_by(model) %>%
                    #   dplyr::summarise(dplyr::across(c(rand_diag, ari_diag, chi_diag, db_diag, pve_diag,
                    #                                    rand_miss, ari_miss, chi_miss, db_miss, pve_miss, R2),
                    #                                  mean, na.rm = TRUE), .groups = "drop")
                    # if ("model_file" %in% names(raw)) {
                    #   file_summary <- raw %>%
                    #     dplyr::group_by(model) %>%
                    #     dplyr::summarise(model_file = paste(stats::na.omit(model_file), collapse = ";"), .groups = "drop")
                    #   metric_summary <- dplyr::left_join(metric_summary, file_summary, by = "model")
                    # }
                    # metric_summary
                    
                    summarise_mi_results(raw)
                  },

                  MI_SF = {
                    imp_list_b <- lapply(imputed_list, function(Xi) {
                      Xi_b <- Xi[ids_b, , drop = FALSE]
                      rownames(Xi_b) <- rownames(X_b)
                      Xi_b
                    })
                    raw <- run_multiple_complete(imp_list_b, layer_vec, fs_diag_labels = fs_diag_b,
                                                 cluster_diag_labels = diag_b, cluster_miss_labels = miss_b,
                                                 clin_df = data.frame(MMSE = mmse_b), models = models,
                                                 fs = TRUE, n_factors = n_factors,
                                                 save_models = save_models, model_dir = model_dir,
                                                 style = style, boot = b)
                    
                    # metric_summary <- raw %>%
                    #   dplyr::group_by(model) %>%
                    #   dplyr::summarise(dplyr::across(c(rand_diag, ari_diag, chi_diag, db_diag, pve_diag,
                    #                                    rand_miss, ari_miss, chi_miss, db_miss, pve_miss, R2),
                    #                                  mean, na.rm = TRUE), .groups = "drop")
                    # if ("model_file" %in% names(raw)) {
                    #   file_summary <- raw %>%
                    #     dplyr::group_by(model) %>%
                    #     dplyr::summarise(model_file = paste(stats::na.omit(model_file), collapse = ";"), .groups = "drop")
                    #   metric_summary <- dplyr::left_join(metric_summary, file_summary, by = "model")
                    # }
                    # metric_summary
                    
                    summarise_mi_results(raw)
                  },

                  {
                    run_NA(X_b, layer_vec, diag_b, miss_b, mmse_b, models, n_factors,
                           save_models = save_models, model_dir = model_dir,
                           style = style, boot = b)
                  }
    )

    # Ensure res is a data frame before adding boot column
    if (!is.data.frame(res)) {
      stop("run_one_boot iteration ", b, " did not return a data frame")
    }

    dplyr::mutate(res, boot = b)
  }

  if (verbose) cat("Starting parallel bootstrap with", mc.cores, "cores...\n")

  # Run parallel bootstrap with better error handling
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
  mc.cores = mc.cores,
  mc.preschedule = FALSE,
  mc.silent = FALSE
)

  # Check for errors
  errors <- sapply(results_list, function(x) inherits(x, "try-error"))
  if (any(errors)) {
    n_errors <- sum(errors)
    warning(n_errors, " bootstrap iterations failed")
    cat("First error:\n")
    print(results_list[[which(errors)[1]]])
    # Keep only successful results
    results_list <- results_list[!errors]
  }

  if (length(results_list) == 0) {
    stop("All bootstrap iterations failed")
  }

  final_results <- dplyr::bind_rows(results_list)

  # Save results if path provided
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

# Extracted from final.R:434
bidifac_svd_counts <- new.env(parent = emptyenv())

# Extracted from final.R:435
bidifac_svd_counts$total <- 0L

# Extracted from final.R:436
bidifac_svd_counts$fallback <- 0L

# Extracted from final.R:437
bidifac_svd_counts$total_by_context <- integer()

# Extracted from final.R:438
bidifac_svd_counts$fallback_by_context <- integer()

# Extracted from final.R:440
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

# Extracted from final.R:458
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

# Extracted from final.R:513
softSVD <- function(X, lambda) {
  svdX <- safe_svd_for_bidifac(X, context = "internal_softSVD")
  nuc <- pmax(svdX$d - lambda, 0)
  if (length(nuc) == 0 || all(nuc <= 0)) {
    out <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  } else {
    keep <- nuc > 0
    out <- svdX$u[, keep, drop = FALSE] %*%
      diag(nuc[keep], nrow = sum(keep)) %*%
      t(svdX$v[, keep, drop = FALSE])
  }
  list(out = out, nuc = sum(nuc))
}

# Extracted from final.R:536

# Extracted from final.R:553
run_model <- function(model_name, X, layer_vec, y_diag, y_miss = NULL,
                      mmse_vec = NULL, n_factors = 3,
                      save_models = FALSE, model_dir = CLEAN_MODEL_DIR,
                      style = "NA", boot = NA_integer_, imp = NULL) {
  out <- switch(
    tolower(model_name),

    # ---- MOFA2 ----
    # MOFA expects one matrix per view Each view is arranged as features by samples We therefore
    # split columns by layer We transpose each block to match MOFA We keep all requested factors
    # We return latent factors and weights
    mofa = {
      message("Starting MOFA")
      # Build view-specific input matrices
      views <- lapply(unique(layer_vec[colnames(X)]), function(lay) {
        cols <- names(layer_vec)[layer_vec == lay]
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(views) <- unique(layer_vec[colnames(X)])
      # Build and train the model
      m <- create_mofa_from_matrix(views)
      data_opts <- get_default_data_options(m)
      data_opts$scale_views <- FALSE
      data_opts$scale_groups <- FALSE
      data_opts$center_groups <- FALSE
      opts <- get_default_model_options(m)
      opts$num_factors <- n_factors
      train_opts <- get_default_training_options(m)
      m <- prepare_mofa(m, data_options = data_opts,
                        model_options = opts, training_options = train_opts)
      m <- run_mofa(m, use_basilisk = TRUE)
      # Extract factors and weights
      Z <- get_factors(m, factors = "all")$group1
      W_list <- get_weights(m, views = "all", factors = "all")
      W_mat <- do.call(rbind, W_list)
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(Z = Z, W = W_mat, fitted_model = m)
    },

    # ---- iClusterPlus ----
    # iClusterPlus expects separate omics matrices Here we assume gene and meta blocks We use
    # Gaussian likelihoods for both We use fixed lambda values as before We return mean latent
    # scores and betas
    icluster = {
      message("Starting iCluster")
      # Split data by modality
      gene_mat <- as.matrix(X[, layer_vec == "gene", drop = FALSE])
      meta_mat <- as.matrix(X[, layer_vec == "meta", drop = FALSE])
      # Fit the model
      fit <- iClusterPlus(dt1 = gene_mat, dt2 = meta_mat,
                          type = c("gaussian", "gaussian"), lambda = c(0.03, 0.03), K = n_factors
      )
      # Build latent matrix
      Z <- fit$meanZ
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      # Build aligned loading matrix
      rownames(fit$beta[[1]]) <- colnames(gene_mat)
      rownames(fit$beta[[2]]) <- colnames(meta_mat)
      W_mat <- rbind(fit$beta[[1]], fit$beta[[2]])
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(Z = Z, W = W_mat, fitted_model = fit)
    },

    # ---- FABIA ----
    # FABIA expects features by samples We therefore transpose the input matrix We keep centering
    # and scaling disabled We match the original factor count We return factor scores and loadings
    fabia = {
      message("Starting FABIA")
      # Fit the model
      f <- fabia(as.matrix(t(X)), p = n_factors, nL = n_factors, center = 0, norm = 0, scale = 0.0)
      # Extract factors and loadings.
      Z <- t(f@Z)
      rownames(Z) <- colnames(f@Z)
      W_mat <- f@L
      rownames(W_mat) <- colnames(X)
      list(Z = Z, W = W_mat, fitted_model = f)
    },

    # ---- SLIDE ----
    # SLIDE requires complete data. Following the example, we mean-impute missing values, center
    # columns, normalize each modality block by its Frobenius norm, select a structure by
    # bi-cross-validation, and fit the selected model.
    slide = {
      message("Starting SLIDE")
      # Optional tuning overrides.
      n_lambda <- get0("slide_n_lambda", ifnotfound = 50L, inherits = TRUE)
      lambda_min <- get0("slide_lambda_min", ifnotfound = 0.01, inherits = TRUE)
      n_fold <- get0("slide_n_fold", ifnotfound = 3L, inherits = TRUE)
      p_fold <- get0("slide_p_fold", ifnotfound = 3L, inherits = TRUE)
      center <- get0("slide_center", ifnotfound = TRUE, inherits = TRUE)
      k_max <- get0("slide_k_max", ifnotfound = 5000L, inherits = TRUE)
      eps <- get0("slide_eps", ifnotfound = 1e-6, inherits = TRUE)
      ratio_max <- get0("slide_ratio_max", ifnotfound = NULL, inherits = TRUE)
      # Build modality blocks and concatenate them.
      blocks <- unique(layer_vec[colnames(X)])
      cols_list <- lapply(blocks, function(b) {
        colnames(X)[layer_vec[colnames(X)] == b]
      })
      pvec <- as.integer(vapply(cols_list, length, integer(1)))
      X_slide <- do.call(cbind, lapply(cols_list, function(cols) {
        as.matrix(X[, cols, drop = FALSE])
      }))
      # Mean-impute missing values before running SLIDE.
      for (j in seq_len(ncol(X_slide))) {
        if (anyNA(X_slide[, j])) {
          mu <- mean(X_slide[, j], na.rm = TRUE)
          if (!is.finite(mu)) mu <- 0
          X_slide[is.na(X_slide[, j]), j] <- mu
        }
      }
      # Standardize as in the example: center columns and normalize each modality block by
      # Frobenius norm.
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
      if (is.null(ratio_max)) {
        ratio_max <- min((n_fold - 1) / n_fold, (p_fold - 1) / p_fold)
      }
      # Select candidate structures.
      out_struct <- SLIDE::create_structure_list(
        X = X_std,
        pvec = pvec,
        lambda_max = max(svec),
        standardized = TRUE,
        lambda_min = lambda_min,
        n_lambda = n_lambda,
        ratio_max = ratio_max
      )
      # Select the best structure by bi-cross-validation.
      out_bcv <- SLIDE::slide_BCV(
        X = X_std,
        pvec = pvec,
        structure_list = out_struct$Slist,
        n_fold = n_fold,
        p_fold = p_fold,
        k_max = k_max,
        eps = eps,
        center = center
      )
      # Fit the selected SLIDE model.
      fit <- SLIDE::slide_givenS(
        X = X_slide,
        pvec = pvec,
        S = out_bcv$structure_min,
        k_max = k_max,
        eps = eps,
        standardized = TRUE
      )
      # Extract sample scores and feature loadings.
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
        Z = Z,
        W = W_mat,
        fitted_model = list(
          standardized_X = X_std,
          structure_list = out_struct$Slist,
          selected_structure = out_bcv$structure_min,
          model = fit
        )
      )
    },

    # ---- JIVE ----
    # JIVE expects one matrix per block. Each block is features by samples. We use the given-rank
    # mode. Joint and individual ranks both use n_factors. We extract the joint scores and
    # loadings.
    jive = {
      message("Starting JIVE")
      # Build view matrices.
      layers <- unique(layer_vec[colnames(X)])
      views <- lapply(layers, function(lay) {
        cols <- which(layer_vec[colnames(X)] == lay)
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(views) <- layers
      # Fit the model.
      jv <- r.jive::jive(
        data = views,
        method = "given",
        rankJ = n_factors,
        rankA = rep(n_factors, length(views)),
        center = FALSE,
        scale = FALSE
      )
      pred <- r.jive::jive.predict(data.new = views, jive.output = jv)
      # Extract joint structure.
      S_mat <- as.matrix(pred$joint.scores)
      L_mat <- as.matrix(pred$joint.load)
      Z <- t(S_mat)
      W_mat <- L_mat
      rownames(Z) <- rownames(X)
      rownames(W_mat) <- colnames(X)
      list(Z = Z, W = W_mat, fitted_model = jv)
    },

    # ---- MCIA ----
    # MCIA expects a list of block matrices Each block is features by samples We keep the original
    # settings here We use CIA factors equal to n_factors We return shared scores and axes
    mcia = {
      message("Starting MCIA")
      # Build the view list.
      views <- lapply(unique(layer_vec[colnames(X)]), function(lay) {
        cols <- which(layer_vec[colnames(X)] == lay)
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(views) <- unique(layer_vec[colnames(X)])
      # Fit the model.
      mcia_res <- omicade4::mcia(df.list = views, cia.nf = n_factors,
                                 cia.scan = FALSE, nsc = TRUE, svd = TRUE)
      # Extract scores and loadings.
      Z <- as.matrix(mcia_res$mcoa$SynVar)
      rownames(Z) <- rownames(X)
      W_mat <- as.matrix(mcia_res$mcoa$axis)
      rownames(W_mat) <- colnames(X)
      list(Z = Z, W = W_mat, fitted_model = mcia_res)
    },

    # ---- RGCCA ----
    # RGCCA expects blocks as samples by features We connect every block to every block We keep
    # scaling disabled throughout We request n_factors components per block We concatenate block
    # score matrices We also pad block loadings to all features
    rgcca = {
      message("Starting RGCCA")
      # Build the input blocks
      blocks <- lapply(unique(layer_vec), function(lay) {
        as.matrix(X[, layer_vec == lay, drop = FALSE])
      })
      names(blocks) <- unique(layer_vec)
      # Build the dense design matrix
      J <- length(blocks)
      Cmat <- 1 - diag(J)
      fit <- RGCCA::rgcca(blocks = blocks, connection = Cmat, # Fit the model
                          ncomp = rep(n_factors, J), scale = FALSE, scale_block = FALSE,
                          quiet = TRUE, verbose = FALSE )
      # Concatenate block variates bc fit$Y gives one score matrix per block, and we want one
      # sample-by-factor matrix Z
      Z_list <- fit$Y
      Z <- do.call(cbind, Z_list)
      rownames(Z) <- rownames(X)
      # Expand block loadings so RGCCA’s separate per-block loading matrices become one full
      # loading matrix
      all_feats <- colnames(X)
      W_list <- vector("list", length(fit$a))
      for (i in seq_along(fit$a)) {
        a_i <- fit$a[[i]]
        rgcca_features <- rownames(a_i)
        W_i <- matrix(
          0, nrow = length(all_feats),
          ncol = ncol(a_i), dimnames = list(all_feats, colnames(a_i))
        )
        W_i[rgcca_features, ] <- a_i
        W_list[[i]] <- W_i
      }
      W_mat <- do.call(cbind, W_list)
      list(Z = Z, W = W_mat, fitted_model = fit)
    },

    # ---- IntNMF ----
    # IntNMF expects nonnegative raw-scale matrices Our working matrix is standardized We
    # therefore back-transform using X_all We split data into modality-specific blocks We fit the
    # shared NMF model We stack block-specific H matrices as loadings
    intnmf = {
      message("Starting IntNMF")
      # Recover approximate raw scale.
      ref <- X_all[, colnames(X), drop = FALSE]
      mu <- colMeans(ref, na.rm = TRUE)
      sigma <- apply(ref, 2, sd, na.rm = TRUE)
      X_raw <- as.matrix(X)
      for (j in seq_len(ncol(X_raw))) {
        X_raw[, j] <- X_raw[, j] * sigma[j] + mu[j]
      }
      # Split data by block
      blocks <- unique(layer_vec)
      dat_list <- lapply(blocks, function(b) {
        cols <- names(layer_vec)[layer_vec == b]
        as.matrix(X_raw[, cols, drop = FALSE])
      })
      names(dat_list) <- blocks
      fit <- nmf.mnnals(dat = dat_list, k = n_factors) # Fit the model
      # Extract scores and loadings.
      Z <- fit$W
      rownames(Z) <- rownames(X_raw)
      W_mat <- do.call(rbind, lapply(fit$H, t))
      rownames(W_mat) <- names(layer_vec)
      list(Z = Z, W = W_mat, fitted_model = fit)
    },

    # ---- BIDIFAC ----
    # BIDIFAC is fit as a linked block model. Here, each omics layer is one column block, and all
    # samples are kept in one row block. Missing modality values stay as NA inside the relevant
    # layer block.
    # impute.BIDIFAC() is then responsible for handling those missing values.
    # The main failure in bootstrap was not caused by the final score extraction. It occurred
    # during BIDIFAC fitting, where BIDIFAC repeatedly calls softSVD(). The patched softSVD()
    # above keeps the BIDIFAC update logic but gives its internal SVD call a fallback if LAPACK's
    # dgesdd routine fails.
    # rmt = FALSE is used because X is already globally centered and scaled before
    # this wrapper. This avoids re-estimating random-matrix scaling inside every bootstrap
    # replicate, which made the BIDIFAC path more fragile.
    bidifac = {
      message("Starting BIDIFAC")
      bl <- unique(layer_vec[colnames(X)])
      cols_list <- lapply(bl, function(b) {
        colnames(X)[layer_vec[colnames(X)] == b]
      })
      names(cols_list) <- bl
      feat_order <- unlist(cols_list, use.names = FALSE)
      dat <- matrix(
        vector("list", length(bl)),
        nrow = 1,
        ncol = length(bl),
        dimnames = list("all_samples", bl)
      )
      for (j in seq_along(bl)) {
        cols <- cols_list[[j]]
        M <- as.matrix(X[, cols, drop = FALSE])
        storage.mode(M) <- "double"
        rownames(M) <- rownames(X)
        colnames(M) <- cols
        if (all(is.na(M))) {
          stop("BIDIFAC: modality block ", bl[j],
               " is entirely NA in this bootstrap sample.")
        }
        if (any(is.nan(M)) || any(is.infinite(M))) {
          stop("BIDIFAC: modality block ", bl[j],
               " contains NaN/Inf.")
        }
        dat[1, j] <- list(M)
      }
      sigma_one <- matrix(1, nrow = 1, ncol = length(bl))

      fit <- if (anyNA(X)) {
        impute.BIDIFAC(dat, rmt = FALSE, sigma = sigma_one, pbar = FALSE)

      } else {
        BIDIFAC(dat, rmt = FALSE, sigma = sigma_one, pbar = FALSE)
      }
      Smat <- data.rearrange(fit$S, rmt = FALSE)$out
      Smat <- as.matrix(Smat)
      storage.mode(Smat) <- "double"
      rownames(Smat) <- rownames(X)
      colnames(Smat) <- feat_order
      Smat <- Smat[rownames(X), colnames(X), drop = FALSE]
      if (any(!is.finite(Smat))) {
        stop("BIDIFAC: Smat contains NA/NaN/Inf after fitting.")
      }
      #sv <- safe_svd_for_bidifac(Smat, context = "final_Smat")
      sv <- svd(Smat)
      k <- min(n_factors, ncol(sv$u), ncol(sv$v), length(sv$d))
      if (k < 1L) {
        stop("BIDIFAC: no usable latent factor returned.")
      }
      Z <- sv$u[, seq_len(k), drop = FALSE] %*%
        diag(sv$d[seq_len(k)], nrow = k)
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(k))
      W_mat <- sv$v[, seq_len(k), drop = FALSE]
      rownames(W_mat) <- colnames(X)
      colnames(W_mat) <- colnames(Z)
      list(Z = Z, W = W_mat, fitted_model = fit)
    },

    # ---- GFA ----
    # GFA expects a named list of blocks Each block is samples by features We keep the default GFA
    # options We request K equal to n_factors We align the loading matrix to X
    gfa = {
      # Build the input list
      layers <- unique(layer_vec[colnames(X)])
      Y <- lapply(layers, function(l) {
        as.matrix(X[, names(layer_vec)[layer_vec == l], drop = FALSE])
      })
      names(Y) <- layers
      fit <- GFA::gfa(Y = Y, K = n_factors, opts = GFA::getDefaultOpts()) # Fit the model
      # Extract scores and loadings
      Z <- as.matrix(fit$X)
      rownames(Z) <- rownames(X)
      W <- as.matrix(fit$W)
      if (is.null(rownames(W))) rownames(W) <- colnames(X)
      W <- W[colnames(X), , drop = FALSE]
      list(Z = Z, W = W, fitted_model = fit)
    },

    # ---- MSFA ----
    # MSFA is intentionally left empty here. The current file preserves that behavior.
    msfa = {
    },

    # ---- SGCCA ----
    # SGCCA expects a list of block matrices Here each block is samples by features We connect all
    # blocks to all blocks We keep all variables in each component We disable scaling to match
    # current setup We expand block loadings to all features
    sgcca = {
      message("Starting SGCCA (mixOmics)")
      # Split columns into blocks
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      Xb <- lapply(blocks, function(cols) {
        as.matrix(X[, cols, drop = FALSE])
      })
      J <- length(Xb)
      design <- 1 - diag(J)
      keepX <- lapply(Xb, function(M) rep(ncol(M), n_factors)) # Keep all variables
      # Fit the model.
      sgcca_args <- list(X = Xb, design = design, ncomp = n_factors, keepX = keepX, scale = FALSE, all.outputs = FALSE)
      if ("scheme" %in% names(formals(mixOmics::wrapper.sgcca))) sgcca_args$scheme <- "horst"
      fit <- do.call(mixOmics::wrapper.sgcca, sgcca_args)
      # Concatenate block variates
      Z <- do.call(cbind, fit$variates)
      rownames(Z) <- rownames(X)
      # Expand loadings to all features.
      allf <- colnames(X)
      W <- do.call(cbind, lapply(seq_along(fit$loadings), function(i) {
        L <- as.matrix(fit$loadings[[i]])
        Wi <- matrix(0, nrow = length(allf), ncol = ncol(L), dimnames = list(allf, colnames(L)))
        Wi[rownames(L), ] <- L
        Wi
      }))
      list(Z = Z, W = W, fitted_model = fit)
    },

    # ---- GPCCA ----
    # GPCCA expects one matrix per block Each block is features by samples The package can use
    # missing values directly We make d equal to n_factors We transpose fit$Z to samples by
    # factors We align fit$W to the original feature order
    gpcca = {
      message("Starting GPCCA")
      # Build block matrices
      blocks <- unique(layer_vec[colnames(X)])
      X_list <- lapply(blocks, function(b) {
        cols <- colnames(X)[layer_vec[colnames(X)] == b]
        t(as.matrix(X[, cols, drop = FALSE]))
      })
      names(X_list) <- blocks
      # Fit the model
      fit <- GPCCA::GPCCAmodel(X.list = X_list, d = n_factors,
                               missVal = any(sapply(X_list, anyNA)), plot_RMSE = FALSE, verbose = 0
      )
      # Extract latent scores
      Z <- as.matrix(t(fit$Z))
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      # Extract and validate loadings
      W <- as.matrix(fit$W)
      feat_order <- unlist(lapply(X_list, rownames), use.names = FALSE)
      if (nrow(W) != length(feat_order)) { stop() }
      # Align loadings to original order
      rownames(W) <- feat_order
      colnames(W) <- colnames(Z)
      W <- W[colnames(X), , drop = FALSE]
      list(Z = Z, W = W, fitted_model = fit)
    },

    # ---- BLOCK PLS ----
    # block.pls expects blocks as samples by features We split X by layer and fit jointly We use
    # the first block as Y We keep the requested component count We expand block loadings to all
    # features
    block_pls = {
      # Build the block list and fit model
      blocks <- unique(layer_vec[colnames(X)])
      Xb <- lapply(blocks, function(b) {
        as.matrix(X[, layer_vec[colnames(X)] == b, drop = FALSE])
      })
      names(Xb) <- blocks
      fit <- mixOmics::block.pls(X = Xb, indY = 1, ncomp = n_factors, scale = FALSE)
      # Concatenate block variates
      Z <- do.call(cbind, fit$variates)
      rownames(Z) <- rownames(X)
      # Expand loadings to all features
      allf <- colnames(X)
      W <- do.call(cbind, lapply(names(fit$loadings), function(b) {
        L <- as.matrix(fit$loadings[[b]])
        Wb <- matrix(0, length(allf), ncol(L),
                     dimnames = list(allf, paste0(b, "_comp", seq_len(ncol(L)))))
        Wb[rownames(L), ] <- L
        Wb
      }))
      list(Z = Z, W = W, fitted_model = fit)
    },

    # ---- O2PLS ----
    # O2PLS here is limited to two blocks Both blocks must be complete We center but do not scale
    # data We cap the joint rank to feasibility We extract X-joint scores and loadings
    o2pls = {
      message("Starting O2PLS")
      # Check the number of blocks
      bl <- unique(layer_vec[colnames(X)])
      if (length(bl) < 2) stop()
      # Build the first two blocks.
      X1 <- as.matrix(X[, layer_vec[colnames(X)] == bl[1], drop = FALSE])
      X2 <- as.matrix(X[, layer_vec[colnames(X)] == bl[2], drop = FALSE])
      #X1 <- scale(X1, center = TRUE, scale = FALSE)
      #X2 <- scale(X2, center = TRUE, scale = FALSE)
      # Keep the rank feasible
      k_use <- min(n_factors, nrow(X1) - 1L, ncol(X1), ncol(X2))
      if (k_use < 1L) stop("O2PLS: not enough rank for >=1 joint component")
      fit <- OmicsPLS::o2m(X1, X2, k_use, 0, 0) # Fit the model
      # Extract joint scores
      Z <- OmicsPLS::scores(fit, which_part = "Xjoint")
      if (is.null(Z) || ncol(Z) < 1) stop("scores() returned empty")
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      # Extract X-joint loadings
      W <- OmicsPLS::loadings(fit, loading_name = "Xjoint")
      W <- as.matrix(W)
      rownames(W) <- colnames(X1)
      colnames(W) <- colnames(Z)
      list(Z = Z, W = W, fitted_model = fit)
    },

    # ---- CIA ----
    # CIA here is limited to two blocks It requires complete cases across both We therefore subset
    # to shared complete rows We run PCA on each block first We then run co-inertia analysis This
    # method does not return loadings here, but we build W by taking the variable coordinates from
    # the co-inertia result for block 1 and block 2, then stacking them back into the original
    # feature order so W is features by factors; this is only an approximation because CIA is a
    # linked ordination of two separate PCA spaces, not a single matrix factorization of X,
    cia = {
      message("Starting CIA")
      # Check the number of blocks and build them
      bl <- unique(layer_vec[colnames(X)])
      if (length(bl) < 2) stop("CIA needs >=2 blocks")
      X1 <- as.matrix(X[, layer_vec[colnames(X)] == bl[1], drop = FALSE])
      X2 <- as.matrix(X[, layer_vec[colnames(X)] == bl[2], drop = FALSE])
      # Run PCA on each block
      d1 <- ade4::dudi.pca(X1, center = FALSE, scale = FALSE, scannf = FALSE, nf = n_factors)
      d2 <- ade4::dudi.pca(X2, center = FALSE, scale = FALSE, scannf = FALSE, nf = n_factors)
      co <- ade4::coinertia(d1, d2, scannf = FALSE, nf = n_factors) # Run co-inertia analysis
      # Extract sample scores
      Z <- as.matrix(co$lX)
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      # Extract block-specific variable coordinates
      W1 <- as.matrix(co$co[, seq_len(ncol(Z)), drop = FALSE])
      W2 <- as.matrix(co$li[, seq_len(ncol(Z)), drop = FALSE])
      # Restore feature names if needed
      if (is.null(rownames(W1))) rownames(W1) <- colnames(X1)
      if (is.null(rownames(W2))) rownames(W2) <- colnames(X2)
      colnames(W1) <- colnames(Z)
      colnames(W2) <- colnames(Z)
      # Stack loadings into original feature order
      W_mat <- matrix( 0, nrow = ncol(X), ncol = ncol(Z), dimnames = list(colnames(X), colnames(Z)))
      W_mat[rownames(W1), ] <- W1
      W_mat[rownames(W2), ] <- W2
      list(Z = Z, W = W_mat, fitted_model = co)
    },

    # ---- NEMO ----
    # NEMO handles partial sample overlap across modalities. It constructs an integrated sample
    # affinity graph and returns sample-level latent factors derived from that graph. NEMO does
    # not estimate original-feature loadings.
    nemo = {
      message("Starting NEMO")
      # Build modality-specific matrices. NEMO expects features by samples.
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      omics_list <- lapply(names(blocks), function(b) {
        cols <- blocks[[b]]
        Xb <- as.matrix(X[, cols, drop = FALSE])
        # A sample is present in this modality if at least one value is observed. Entirely missing
        # modality rows remain excluded from that modality.
        observed <- rowSums(!is.na(Xb)) > 0
        X_obs <- Xb[observed, , drop = FALSE]
        if (nrow(X_obs) < 2L) {
          stop(paste0("NEMO: modality ", b, " contains fewer than two observed samples."))
        }
        # Remove features that are unusable within the observed samples of this modality.
        keep_features <- colSums(!is.na(X_obs)) > 0
        X_obs <- X_obs[, keep_features, drop = FALSE]
        if (ncol(X_obs) < 1L) {
          stop(paste0("NEMO: modality ", b, " contains no usable features."))
        }
        # NEMO cannot take remaining cell-level NAs inside an observed modality. Impute only these
        # within-modality missing cells. This does not impute an entirely missing modality for a
        # participant.
        for (j in seq_len(ncol(X_obs))) {
          if (anyNA(X_obs[, j])) {
            mu <- mean(X_obs[, j], na.rm = TRUE)
            X_obs[is.na(X_obs[, j]), j] <- mu
          }
        }
        # NEMO expects features by samples.
        t(X_obs)
      })
      names(omics_list) <- names(blocks)
      # Match the example's neighborhood parameter, while keeping it feasible.
      k_requested <- get0("nemo_k", ifnotfound = 10L, inherits = TRUE)
      k_max_available <- min(vapply(omics_list, ncol, integer(1))) - 1L
      if (k_max_available < 1L) {
        stop("NEMO: insufficient observed samples.")
      }
      k_use <- min(as.integer(k_requested), k_max_available)
      # Construct the integrated NEMO affinity graph.
      affinity <- NEMO::nemo.affinity.graph(omics_list, k = k_use)
      # Derive latent sample factors from the affinity graph, following the spectral embedding
      # used in the example.
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
      # Reorder the NEMO factors to match the rows of X.
      missing_from_graph <- setdiff(rownames(X), rownames(Z))
      if (length(missing_from_graph) > 0L) {
        stop("NEMO: some samples are absent from every usable omics block.")
      }
      Z <- Z[rownames(X), , drop = FALSE]
      list(Z = Z, W = NULL, fitted_model = list(affinity_graph = affinity, k = k_use))
    },

    # ---- BSF ----
    # BSF expects one matrix per source, with features as rows and samples as columns. Therefore,
    # each modality block is transposed before fitting. The bootstrap failure for BSF was
    # connected to BIDIFAC. BSFP can use BIDIFAC internally for initialization and rank selection.
    # That means BSF can fail even when the visible model name is "bsf", because BIDIFAC is being
    # called inside the BSF fitting path. This wrapper makes BSF independent of BIDIFAC by doing
    # two things: 1. ranks are supplied directly, 2. nninit = FALSE disables BIDIFAC-based
    # initialization. This is consistent with the benchmark setup because the wrapper already
    # fixes n_factors for all methods. So BSF is evaluated with the same requested factor count
    # instead of using BIDIFAC to choose ranks separately.
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

      ranks_use <- c(n_factors, rep(0L, length(data_list)))
      fit <- BSFP::bsfp(
        data = data_list,
        Y = NULL,
        nninit = FALSE,
        ranks = ranks_use,
        nsample = 1000,
        progress = FALSE,
        save_structures = FALSE,
        save_loadings_scores = FALSE,
        save_predictive_model = FALSE,
        save_imputations = FALSE,
        save_init = FALSE,
        save_last_sample = TRUE
      )
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
        0,
        nrow = length(all_feats),
        ncol = ncol(Z),
        dimnames = list(all_feats, colnames(Z))
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

    # ---- GIPCA ----
    # GIPCA is loaded from external source files We use the joint-only model here Individual ranks
    # are all fixed at zero Family is Gaussian for every block We extract the joint score matrix
    # We align joint loadings back to X
    gipca = {
      message("Starting GIPCA")
      load_gipca_once() # Load the source functions
      # Build the data list.
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      data_list <- lapply(blocks, function(cols) {
        as.matrix(X[, cols, drop = FALSE])
      })
      # Fit the joint-only model
      fit <- EPCAJIVEMissbio(
        data = data_list, rankj = n_factors, ranka = rep(0L, length(data_list)),
        D = vapply(data_list, ncol, integer(1)), family = rep("gaussian", length(data_list)),
        tol = 0.1, max.iter = 500, n.size = NULL, lambda = 0
      )
      # Extract joint scores
      Z <- as.matrix(fit$U.joint[, -1, drop = FALSE])
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      # Extract and align joint loadings
      W_mat <- t(as.matrix(fit$V.joint[-1, , drop = FALSE]))
      rownames(W_mat) <- unlist(blocks, use.names = FALSE)
      colnames(W_mat) <- colnames(Z)
      W_mat <- W_mat[colnames(X), , drop = FALSE]
      list(Z = Z, W = W_mat, fitted_model = fit)
    },

    # ---- tICA ----
    # tICA is fit to a block-by-sample-by-sample tensor. Each block slice is a sample correlation
    # matrix.
    # Entry (i, j) is the correlation between samples i and j within one block.
    # This uses both omics blocks jointly in one tensor model. Because the tensor features are
    # sample positions, not original variables, DoTICA does not return native loadings on the
    # original X columns here. We therefore build an approximate feature loading matrix afterward
    # by regressing the original feature matrix on the final sample scores.
    tica = {
      message("Starting tICA")
      # Build block metadata This step records how the columns of X are grouped into omics blocks
      # We need that grouping because tICA is applied to one tensor slice per block, not to the
      # full feature matrix all at once The object `blocks` stores the feature names for each
      # modality B gives the number of modality blocks, n gives the number of samples
      blocks <- split(colnames(X), layer_vec[colnames(X)])
      B <- length(blocks)
      n <- nrow(X)
      data_a <- array(0, dim = c(B, n, n)) # Allocate the sample-correlation tensor
      # Build one sample-by-sample correlation slice per block so each modality contributes its
      # own similarity matrix to the tensor used by tICA.
      for (b in seq_len(B)) {
        cols_b <- blocks[[b]]
        Xb <- as.matrix(X[, cols_b, drop = FALSE])
        # Correlate sample profiles within this block.
        Cb <- suppressWarnings(stats::cor(t(Xb), use = "everything"))
        # Replace undefined correlations deterministically.
        Cb[!is.finite(Cb)] <- 0
        diag(Cb) <- 1
        Cb <- (Cb + t(Cb)) / 2
        data_a[b, , ] <- Cb
      }
      # Fit the model
      dim_vec <- rep(n_factors, B)
      fit <- tensorICA::DoTICA(Data = data_a, dim = dim_vec, method = "FOBI")
      # Convert each projected slice to sample-by-component form so all block-specific outputs
      # share one consistent orientation before they are combined.
      proj_list <- lapply(fit$TICA.projS, function(M) {
        M <- as.matrix(M)
        if (ncol(M) == n) { t(M) }
        else if (nrow(M) == n) { M }
        else { stop() }
      })
      # Combine all block-specific projections.
      Z_full <- do.call(cbind, proj_list)
      Z <- Z_full[, seq_len(min(n_factors, ncol(Z_full))), drop = FALSE]
      rownames(Z) <- rownames(X)
      colnames(Z) <- paste0("Factor", seq_len(ncol(Z)))
      # Build approximate feature loadings from X ~ Z.
      Zc <- scale(Z, center = TRUE, scale = FALSE)
      Xc <- scale(as.matrix(X), center = TRUE, scale = FALSE)
      Bhat <- qr.solve(Zc, Xc)
      W_mat <- t(Bhat)
      rownames(W_mat) <- colnames(X)
      colnames(W_mat) <- colnames(Z)
      list(Z = Z, W = W_mat, fitted_model = fit)
    }
  )
  # Pull core model outputs.
  Z <- out$Z
  W <- out$W
  # Save the fitted object if requested.
  model_file <- NA_character_
  if (save_models) {
    fit_obj <- if ("fitted_model" %in% names(out)) out$fitted_model else NULL
    save_obj <- list(
      model = model_name,
      style = toupper(style),
      boot = boot,
      imp = imp,
      n_factors = n_factors,
      Z = Z,
      W = W,
      fitted_model = fit_obj,
      rows = rownames(X),
      cols = colnames(X)
    )
    model_file <- save_boot_model(
      obj = save_obj,
      model_name = model_name,
      style = style,
      boot = boot,
      imp = imp,
      base_dir = model_dir
    )
  }
  # Score diagnosis clustering.
  cl_diag <- score_clusters(Z, y_diag, X_original = X, W_loadings = W, seed = 100 + as.integer(boot))
  # Score MMSE regression if available.
  R2_val <- if (!is.null(mmse_vec)) score_regression(Z, mmse_vec) else NA_real_
  # Build the main result list.
  result <- list(
    model = model_name,
    model_file = model_file,
    rand_diag = cl_diag["rand"],
    ari_diag = cl_diag["ari"],
    chi_diag = cl_diag["chi"],
    db_diag = cl_diag["db"],
    pve_diag = cl_diag["pve"],
    R2 = R2_val
  )
  # Score missingness clustering
  if (!is.null(y_miss)) {
    cl_miss <- score_clusters(Z, y_miss, X_original = X, W_loadings = W, seed = 200 + as.integer(boot))
    result$rand_miss <- cl_miss["rand"]
    result$ari_miss <- cl_miss["ari"]
    result$chi_miss <- cl_miss["chi"]
    result$db_miss <- cl_miss["db"]
    result$pve_miss <- cl_miss["pve"]
  }
  return(result)
}

# Extracted from final.R:1501
score_clusters <- function(latent_matrix, true_labels, X_original = NULL, W_loadings = NULL, seed = 100) {
  if (is.null(rownames(latent_matrix))) rownames(latent_matrix) <- names(true_labels) %||% seq_len(nrow(latent_matrix))
  ids <- intersect(rownames(latent_matrix), names(true_labels))
  latent_matrix <- latent_matrix[ids, , drop = FALSE]
  true_labels <- true_labels[ids]
  stopifnot(nrow(latent_matrix) == length(true_labels))
  k <- nlevels(droplevels(true_labels))
  set.seed(seed)
  pred <- kmeans(latent_matrix, centers = k, nstart = 25)$cluster
  rand <- clusterCrit::extCriteria(as.integer(true_labels), pred, "Rand")$rand
  ari <- mclust::adjustedRandIndex(true_labels, pred)
  int <- clusterCrit::intCriteria(as.matrix(latent_matrix), pred, c("Calinski_Harabasz", "Davies_Bouldin"))
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

# Extracted from final.R:1533
score_regression <- function(latent_matrix, mmse_vec) {
  ids <- intersect(rownames(latent_matrix), names(mmse_vec))
  df <- data.frame(MMSE = mmse_vec[ids], latent_matrix[ids, , drop = FALSE])
  fit <- lm(MMSE ~ ., data = df)
  summary(fit)$r.squared
}

# Extracted from final.R:1543
sel_SIS <- function(X_complete, label_vec) {
  keep <- !is.na(label_vec)
  ybin <- as.integer(label_vec[keep] == "AD")
  sis <- SIS::SIS(x = as.matrix(X_complete[keep, , drop = FALSE]), y = ybin, family = "binomial", tune = "bic", nsis = ncol(X_complete))
  colnames(X_complete)[sis$sis.ix0]
}

# Extracted from final.R:1568
run_complete <- function(
    X,
    layer_vec,
    diag_labels,
    miss_labels,
    mmse_vec,
    models,
    n_factors = 3,
    save_models = FALSE,
    model_dir = CLEAN_MODEL_DIR,
    style = "NA",
    boot = NA_integer_,
    imp = NULL
) {

  dplyr::bind_rows(
    lapply(
      models,
      function(m) {

        started <- proc.time()[["elapsed"]]

        tryCatch(
          {
            out <- run_model(
              model_name = m,
              X = X,
              layer_vec = layer_vec[colnames(X)],
              y_diag = diag_labels[rownames(X)],
              y_miss = miss_labels[rownames(X)],
              mmse_vec = if (!is.null(mmse_vec)) {
                mmse_vec[rownames(X)]
              } else {
                NULL
              },
              n_factors = n_factors,
              save_models = save_models,
              model_dir = model_dir,
              style = style,
              boot = boot,
              imp = imp
            )

            result <- tibble::as_tibble(out)

            if (!"model_file" %in% names(result)) {
              result$model_file <- NA_character_
            }

            for (metric in metric_columns) {
              if (!metric %in% names(result)) {
                result[[metric]] <- NA_real_
              }
            }

            result %>%
              dplyr::mutate(
                status = "success",
                error_message = NA_character_,
                n_attempted = 1L,
                n_successful = 1L,
                n_failed = 0L,
                elapsed_seconds =
                  proc.time()[["elapsed"]] - started
              )
          },

          error = function(e) {

            message(
              "MODEL FAILED",
              " | model=", m,
              " | style=", style,
              " | boot=", boot,
              if (!is.null(imp)) {
                paste0(" | imp=", imp)
              } else {
                ""
              },
              " | ", conditionMessage(e)
            )

            tibble::tibble(
              model = m,
              model_file = NA_character_,
              rand_diag = NA_real_,
              ari_diag = NA_real_,
              chi_diag = NA_real_,
              db_diag = NA_real_,
              pve_diag = NA_real_,
              rand_miss = NA_real_,
              ari_miss = NA_real_,
              chi_miss = NA_real_,
              db_miss = NA_real_,
              pve_miss = NA_real_,
              R2 = NA_real_,
              status = "failed",
              error_message = conditionMessage(e),
              n_attempted = 1L,
              n_successful = 0L,
              n_failed = 1L,
              elapsed_seconds =
                proc.time()[["elapsed"]] - started
            )
          }
        )
      }
    )
  )
}

# Extracted from final.R:1688

# Extracted from final.R:1709

