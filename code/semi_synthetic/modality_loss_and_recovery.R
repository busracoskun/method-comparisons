#!/usr/bin/env Rscript

# =============================================================================
# SEMI-SYNTHETIC INCOMPLETE-MODALITY EXPERIMENT: ADNI + TCGA
# =============================================================================
#
# This script preserves the settings used by the older final LSF workflow:
#   - 600 experiment tasks
#   - 100 repetitions
#   - MCAR, MAR, and MNAR
#   - 20% and 70% target missingness
#   - 3 latent factors
#   - beta = 2
#   - missForest ntree = 100 and maxiter = 10
#   - methods supplied through SEMISYNTH_METHODS
#   - SAVE_SEMISYNTH_DATA = 0 by default
#   - RETRY_FAILED = 0 by default
#
# Task mapping:
#     1-100   = MCAR, 20%, repetitions 1-100
#   101-200   = MAR,  20%, repetitions 1-100
#   201-300   = MNAR, 20%, repetitions 1-100
#   301-400   = MCAR, 70%, repetitions 1-100
#   401-500   = MAR,  70%, repetitions 1-100
#   501-600   = MNAR, 70%, repetitions 1-100
#
# Every task runs the selected condition on:
#   - every ADNI modality
#   - every available TCGA modality
#
# Dataset preparation is loaded from the local datasets/ modules.
#
# Modes:
#   SEMISYNTH_MODE=ground_truth  Fit complete-data reference models once.
#   SEMISYNTH_MODE=run           Run one SEMISYNTH_TASK_ID from 1 through 600.
#   SEMISYNTH_MODE=combine       Combine completed task outputs.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# -----------------------------------------------------------------------------
# PATHS AND OLDER-LSF SETTINGS
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv(
  "PROJECT_ROOT",
  unset = Sys.getenv("CLEAN_ROOT")
)

PROJECT_DIR <- Sys.getenv("BENCH_WORK_DIR")

OUT_DIR <- Sys.getenv(
  "SEMISYNTH_OUT_DIR",
  unset = file.path(PROJECT_DIR, "results_semi_synthetic_final_v1")
)

MODE <- tolower(Sys.getenv("SEMISYNTH_MODE", unset = "run"))

SEMISYNTH_METHODS <- Sys.getenv(
  "SEMISYNTH_METHODS",
  unset = paste(
    c(
      "icluster", "fabia", "mcia", "intnmf", "gfa", "sgcca",
      "o2pls", "cia", "tica", "slide", "jive"
    ),
    collapse = ","
  )
)

SEMISYNTH_REPETITIONS <- as.integer(
  Sys.getenv("SEMISYNTH_REPETITIONS", unset = "100")
)

SEMISYNTH_N_FACTORS <- as.integer(
  Sys.getenv("SEMISYNTH_N_FACTORS", unset = "3")
)

SEMISYNTH_BETA <- as.numeric(
  Sys.getenv("SEMISYNTH_BETA", unset = "2")
)

MISSFOREST_NTREE <- as.integer(
  Sys.getenv("MISSFOREST_NTREE", unset = "100")
)

MISSFOREST_MAXITER <- as.integer(
  Sys.getenv("MISSFOREST_MAXITER", unset = "10")
)

SAVE_SEMISYNTH_DATA <- identical(
  Sys.getenv("SAVE_SEMISYNTH_DATA", unset = "0"),
  "1"
)

RETRY_FAILED <- identical(
  Sys.getenv("RETRY_FAILED", unset = "0"),
  "1"
)

GROUND_TRUTH_FILE <- file.path(
  OUT_DIR,
  "ground_truth",
  "ground_truth_fits.rds"
)

GROUND_TRUTH_STATUS_FILE <- file.path(
  OUT_DIR,
  "ground_truth",
  "ground_truth_status.csv"
)

TASK_DIR <- file.path(OUT_DIR, "tasks")
LOCK_DIR <- file.path(OUT_DIR, "task_locks")
DATA_DIR <- file.path(OUT_DIR, "simulated_data")
COMBINED_DIR <- file.path(OUT_DIR, "combined")

for (directory in c(
  OUT_DIR,
  dirname(GROUND_TRUTH_FILE),
  TASK_DIR,
  LOCK_DIR,
  DATA_DIR,
  COMBINED_DIR
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------------------------------------------------------
# VALIDATION AND GENERAL HELPERS
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

parse_csv <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character())
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  unique(tolower(values[nzchar(values)]))
}

METHODS <- parse_csv(SEMISYNTH_METHODS)

if (!length(METHODS)) {
  stop("SEMISYNTH_METHODS did not contain any model names.")
}
if (is.na(SEMISYNTH_REPETITIONS) || SEMISYNTH_REPETITIONS != 100L) {
  stop("This final design requires SEMISYNTH_REPETITIONS=100.")
}
if (is.na(SEMISYNTH_N_FACTORS) || SEMISYNTH_N_FACTORS < 1L) {
  stop("SEMISYNTH_N_FACTORS must be at least 1.")
}
if (!is.finite(SEMISYNTH_BETA) || SEMISYNTH_BETA < 0) {
  stop("SEMISYNTH_BETA must be nonnegative.")
}
if (is.na(MISSFOREST_NTREE) || MISSFOREST_NTREE < 1L) {
  stop("MISSFOREST_NTREE must be at least 1.")
}
if (is.na(MISSFOREST_MAXITER) || MISSFOREST_MAXITER < 1L) {
  stop("MISSFOREST_MAXITER must be at least 1.")
}

required_packages <- c("dplyr", "missForest", "T4transport", "FactoMineR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

remove_triple_quote_blocks <- function(lines) {
  marker <- which(trimws(lines) == "'''")

  if (!length(marker)) return(lines)
  if (length(marker) %% 2L != 0L) {
    stop("An unmatched ''' marker was found in a setup file.")
  }

  keep <- rep(TRUE, length(lines))
  for (i in seq(1L, length(marker), by = 2L)) {
    keep[marker[i]:marker[i + 1L]] <- FALSE
  }
  lines[keep]
}

rbind_fill <- function(objects) {
  objects <- objects[!vapply(objects, is.null, logical(1))]
  if (!length(objects)) return(data.frame())
  dplyr::bind_rows(objects)
}

safe_sd <- function(x) {
  value <- suppressWarnings(stats::sd(x, na.rm = TRUE))
  if (is.finite(value)) value else NA_real_
}

# -----------------------------------------------------------------------------
# SOURCE THE MODEL AND DATA SETUP WITHOUT RUNNING BOOTSTRAPS
# -----------------------------------------------------------------------------

source_dataset_setup <- function(dataset) {
  load_dataset(dataset, observed=FALSE)
}

make_latent_runner <- function(run_model_function) {
  body_parts <- as.list(body(run_model_function))

  w_assignment <- which(vapply(
    body_parts,
    function(x) identical(x, quote(W <- out$W)),
    logical(1)
  ))

  if (length(w_assignment) != 1L) {
    stop(
      "Could not construct a latent-output runner from run_model(). ",
      "The expected 'W <- out$W' line was not found exactly once."
    )
  }

  latent_runner <- run_model_function
  body(latent_runner) <- as.call(c(
    list(as.name("{")),
    body_parts[2:w_assignment],
    list(quote(return(list(model = model_name, Z = Z, W = W))))
  ))

  latent_runner
}

prepare_dataset <- function(dataset) {
  env <- source_dataset_setup(dataset)
  X_stand <- as.matrix(env$X_stand)
  layer_vec <- env$layer_vec[colnames(X_stand)]

  truth_ids <- intersect(env$common_ids, rownames(X_stand))
  X_truth <- X_stand[truth_ids, , drop = FALSE]

  # Match the complete-case logic of the existing pipelines: retain only
  # columns that are fully observed among participants present in every block.
  X_truth <- X_truth[
    ,
    colSums(is.na(X_truth)) == 0L,
    drop = FALSE
  ]

  feature_sds <- vapply(seq_len(ncol(X_truth)), function(j) safe_sd(X_truth[, j]), numeric(1))
  keep <- is.finite(feature_sds) & feature_sds > 0
  X_truth <- X_truth[, keep, drop = FALSE]
  layer_vec <- layer_vec[colnames(X_truth)]

  if (nrow(X_truth) < 5L) {
    stop(dataset, " has fewer than five usable complete-case patients.")
  }
  if (ncol(X_truth) < SEMISYNTH_N_FACTORS) {
    stop(dataset, " has too few complete-case features.")
  }
  if (anyNA(X_truth) || any(!is.finite(X_truth))) {
    stop(dataset, " complete-case truth contains NA/NaN/Inf.")
  }

  diag_labels <- factor(env$diag_labels[rownames(X_truth)])
  names(diag_labels) <- rownames(X_truth)

  clinical_vec <- NULL
  if (exists("mmse_vec", envir = env, inherits = FALSE)) {
    clinical_vec <- env$mmse_vec[rownames(X_truth)]
    names(clinical_vec) <- rownames(X_truth)
  }

  missing_labels <- factor(rep("CC", nrow(X_truth)), levels = "CC")
  names(missing_labels) <- rownames(X_truth)

  modalities <- unique(unname(layer_vec))

  unsupported <- character()
  if (dataset == "tcga" && length(modalities) > 2L) {
    unsupported <- intersect(METHODS, c("o2pls", "cia"))
  }

  list(
    dataset = dataset,
    dataset_label = if (dataset == "adni") "ADNI" else "TCGA_GBMLGG",
    env = env,
    latent_runner = make_latent_runner(env$run_model),
    X_truth = X_truth,
    layer_vec = layer_vec,
    diag_labels = diag_labels,
    clinical_vec = clinical_vec,
    missing_labels = missing_labels,
    modalities = modalities,
    unsupported_methods = unsupported
  )
}

# -----------------------------------------------------------------------------
# MISSINGNESS GENERATION
# -----------------------------------------------------------------------------

make_missing <- function(
    X,
    block_cols,
    mechanism = c("MCAR", "MAR", "MNAR"),
    target_rate = 0.3,
    observed_score = NULL,
    latent_score = NULL,
    beta = 2,
    seed = 1
) {
  mechanism <- match.arg(mechanism)
  X <- as.matrix(X)
  n <- nrow(X)

  if (is.character(block_cols)) {
    block_cols <- match(block_cols, colnames(X))
  }
  block_cols <- block_cols[!is.na(block_cols)]

  if (!length(block_cols)) {
    stop("No valid block columns were supplied.")
  }

  mech <- switch(
    mechanism,
    MCAR = rep(0, n),
    MAR = {
      if (is.null(observed_score)) stop("MAR requires observed_score.")
      beta * as.numeric(scale(observed_score))
    },
    MNAR = {
      if (is.null(latent_score)) stop("MNAR requires latent_score.")
      beta * as.numeric(scale(latent_score))
    }
  )
  mech[!is.finite(mech)] <- 0

  alpha <- stats::uniroot(
    function(a) mean(stats::plogis(a + mech)) - target_rate,
    interval = c(-50, 50)
  )$root

  probability <- stats::plogis(alpha + mech)

  set.seed(seed)
  missing <- stats::rbinom(n, size = 1, prob = probability) == 1
  names(missing) <- rownames(X)
  names(probability) <- rownames(X)

  X_miss <- X
  X_miss[missing, block_cols] <- NA_real_

  list(
    X_miss = X_miss,
    missing = missing,
    probability = probability,
    mechanism = mechanism,
    target_rate = target_rate,
    realized_rate = mean(missing)
  )
}

block_pc1_score <- function(X_block) {
  X_block <- as.matrix(X_block)

  score <- if (ncol(X_block) == 1L) {
    as.numeric(X_block[, 1L])
  } else {
    tryCatch(
      stats::prcomp(X_block, center = TRUE, scale. = FALSE, rank. = 1L)$x[, 1L],
      error = function(e) rowMeans(X_block)
    )
  }

  average_profile <- rowMeans(X_block)
  direction_cor <- suppressWarnings(stats::cor(score, average_profile))
  if (is.finite(direction_cor) && direction_cor < 0) score <- -score

  score <- as.numeric(scale(score))
  score[!is.finite(score)] <- 0
  score
}

make_mar_score <- function(X, layer_vec, dropped_modality) {
  observed_modalities <- setdiff(unique(unname(layer_vec)), dropped_modality)
  if (!length(observed_modalities)) {
    stop("No observed modality remains for MAR.")
  }

  scores <- vapply(
    observed_modalities,
    function(modality) {
      cols <- names(layer_vec)[layer_vec == modality]
      block_pc1_score(X[, cols, drop = FALSE])
    },
    numeric(nrow(X))
  )

  if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1L)
  score <- rowMeans(scores)
  names(score) <- rownames(X)
  score
}

make_mnar_score <- function(X, layer_vec, dropped_modality) {
  cols <- names(layer_vec)[layer_vec == dropped_modality]
  score <- block_pc1_score(X[, cols, drop = FALSE])
  names(score) <- rownames(X)
  score
}

# -----------------------------------------------------------------------------
# LATENT METRICS
# -----------------------------------------------------------------------------

latent_metrics <- function(Z0, Z1, W0 = NULL, W1 = NULL, block = NULL, gw = TRUE, top_k = 50L) {
  Z0 <- as.matrix(Z0)
  Z1 <- as.matrix(Z1)

  if (!is.null(rownames(Z0)) && !is.null(rownames(Z1))) {
    ids <- intersect(rownames(Z0), rownames(Z1))
    Z0 <- Z0[ids, , drop = FALSE]
    Z1 <- Z1[ids, , drop = FALSE]
  }

  k_z <- min(ncol(Z0), ncol(Z1))
  if (k_z < 1L) stop("No shared latent factors were returned.")

  Z0 <- Z0[, seq_len(k_z), drop = FALSE]
  Z1 <- Z1[, seq_len(k_z), drop = FALSE]

  relerr <- function(A, B) {
    denominator <- norm(A, type = "F")
    if (!is.finite(denominator) || denominator <= 0) {
      NA_real_
    } else {
      norm(A - B, type = "F") / denominator
    }
  }

  output <- c(
    rv = FactoMineR::coeffRV(Z0, Z1)$rv,
    gw = NA_real_,
    recon_relerr_full = NA_real_,
    recon_cor = NA_real_,
    rv_W = NA_real_,
    wgene_rankcor = NA_real_,
    wgene_jaccard = NA_real_
  )

  if (isTRUE(gw)) {
    output["gw"] <- T4transport::gwdist(stats::dist(Z0), stats::dist(Z1))$distance
  }

  if (!is.null(W0) && !is.null(W1)) {
    W0 <- as.matrix(W0)
    W1 <- as.matrix(W1)

    features <- if (!is.null(rownames(W0)) && !is.null(rownames(W1))) {
      intersect(rownames(W0), rownames(W1))
    } else {
      seq_len(min(nrow(W0), nrow(W1)))
    }

    if (length(features) >= 2L) {
      W0 <- W0[features, , drop = FALSE]
      W1 <- W1[features, , drop = FALSE]
      k_w <- min(ncol(W0), ncol(W1), ncol(Z0), ncol(Z1))

      if (k_w >= 1L) {
        W0 <- W0[, seq_len(k_w), drop = FALSE]
        W1 <- W1[, seq_len(k_w), drop = FALSE]

        Xhat0 <- Z0[, seq_len(k_w), drop = FALSE] %*% t(W0)
        Xhat1 <- Z1[, seq_len(k_w), drop = FALSE] %*% t(W1)
        colnames(Xhat0) <- colnames(Xhat1) <- as.character(features)

        output["recon_relerr_full"] <- relerr(Xhat0, Xhat1)
        output["recon_cor"] <- suppressWarnings(stats::cor(as.numeric(Xhat0), as.numeric(Xhat1)))
        output["rv_W"] <- FactoMineR::coeffRV(W0, W1)$rv

        importance0 <- sqrt(rowSums(W0^2))
        importance1 <- sqrt(rowSums(W1^2))
        output["wgene_rankcor"] <- suppressWarnings(
          stats::cor(importance0, importance1, method = "spearman")
        )

        if (!is.null(rownames(W0))) {
          K <- min(as.integer(top_k), length(features))
          names(importance0) <- rownames(W0)
          names(importance1) <- rownames(W1)
          top0 <- names(sort(importance0, decreasing = TRUE))[seq_len(K)]
          top1 <- names(sort(importance1, decreasing = TRUE))[seq_len(K)]
          output["wgene_jaccard"] <-
            length(intersect(top0, top1)) / length(union(top0, top1))
        }

        if (!is.null(block)) {
          block_use <- block[as.character(features)]
          for (modality in unique(block_use[!is.na(block_use)])) {
            columns <- which(block_use == modality)
            output[paste0("recon_relerr_", modality)] <- relerr(
              Xhat0[, columns, drop = FALSE],
              Xhat1[, columns, drop = FALSE]
            )
          }
        }
      }
    }
  }

  output
}

# -----------------------------------------------------------------------------
# MODEL FITTING
# -----------------------------------------------------------------------------

fit_latent <- function(setup, model, X, seed, run_id, style) {
  started <- proc.time()[["elapsed"]]
  set.seed(seed)

  tryCatch(
    {
      result <- setup$latent_runner(
        model_name = model,
        X = X,
        layer_vec = setup$layer_vec[colnames(X)],
        y_diag = setup$diag_labels[rownames(X)],
        y_miss = setup$missing_labels[rownames(X)],
        mmse_vec = if (!is.null(setup$clinical_vec)) {
          setup$clinical_vec[rownames(X)]
        } else {
          NULL
        },
        n_factors = SEMISYNTH_N_FACTORS,
        save_models = FALSE,
        style = style,
        boot = run_id,
        imp = NULL
      )

      Z <- as.matrix(result$Z)
      if (is.null(rownames(Z))) rownames(Z) <- rownames(X)
      Z <- Z[rownames(X), , drop = FALSE]

      if (!ncol(Z) || any(!is.finite(Z))) {
        stop("The fitted latent matrix is empty or non-finite.")
      }

      W <- result$W
      if (!is.null(W)) {
        W <- as.matrix(W)
        if (any(!is.finite(W))) {
          stop("The fitted loading matrix contains NA/NaN/Inf.")
        }
      }

      list(
        status = "success",
        error_message = NA_character_,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        Z = Z,
        W = W
      )
    },
    error = function(e) {
      list(
        status = "failed",
        error_message = conditionMessage(e),
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        Z = NULL,
        W = NULL
      )
    }
  )
}

compare_fit <- function(baseline, fitted, layer_vec) {
  if (!identical(baseline$status, "success")) {
    return(list(
      status = "baseline_failed",
      error_message = baseline$error_message
    ))
  }

  if (!identical(fitted$status, "success")) {
    return(list(
      status = "fit_failed",
      error_message = fitted$error_message
    ))
  }

  tryCatch(
    {
      values <- latent_metrics(
        Z0 = baseline$Z,
        Z1 = fitted$Z,
        W0 = baseline$W,
        W1 = fitted$W,
        block = layer_vec,
        gw = TRUE,
        top_k = 50L
      )

      c(
        list(status = "success", error_message = NA_character_),
        as.list(values)
      )
    },
    error = function(e) {
      list(
        status = "metric_failed",
        error_message = conditionMessage(e)
      )
    }
  )
}

# -----------------------------------------------------------------------------
# GROUND-TRUTH MODE
# -----------------------------------------------------------------------------

build_ground_truth <- function() {
  cat("============================================================\n")
  cat("BUILDING ADNI + TCGA GROUND-TRUTH FITS\n")
  cat("Methods:", paste(METHODS, collapse = ", "), "\n")
  cat("Latent factors:", SEMISYNTH_N_FACTORS, "\n")
  cat("Output:", GROUND_TRUTH_FILE, "\n")
  cat("============================================================\n")

  ground_truth <- list(
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    methods = METHODS,
    n_factors = SEMISYNTH_N_FACTORS,
    datasets = list()
  )

  status_rows <- list()
  status_index <- 0L

  for (dataset in strsplit(Sys.getenv("BENCH_DATASET"), ",", fixed=TRUE)[[1]]) {
    setup <- prepare_dataset(dataset)
    fits <- setNames(vector("list", length(METHODS)), METHODS)

    cat(
      setup$dataset_label,
      "complete-case patients:", nrow(setup$X_truth),
      "features:", ncol(setup$X_truth),
      "modalities:", paste(setup$modalities, collapse = ", "),
      "\n"
    )

    for (method_index in seq_along(METHODS)) {
      model <- METHODS[method_index]

      if (model %in% setup$unsupported_methods) {
        fit <- list(
          status = "unsupported",
          error_message = paste0(
            model,
            " is a two-block wrapper and is not run on TCGA with ",
            length(setup$modalities),
            " modalities."
          ),
          elapsed_seconds = 0,
          Z = NULL,
          W = NULL
        )
      } else {
        cat("Ground truth:", setup$dataset_label, model, "\n")
        fit <- fit_latent(
          setup = setup,
          model = model,
          X = setup$X_truth,
          seed = 10000L + method_index,
          run_id = 0L,
          style = "SIM_GROUND_TRUTH"
        )
      }

      fits[[model]] <- fit
      status_index <- status_index + 1L
      status_rows[[status_index]] <- data.frame(
        dataset = setup$dataset_label,
        model = model,
        status = fit$status,
        error_message = fit$error_message,
        elapsed_seconds = fit$elapsed_seconds,
        n_samples = nrow(setup$X_truth),
        n_features = ncol(setup$X_truth),
        modalities = paste(setup$modalities, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }

    ground_truth$datasets[[dataset]] <- list(
      dataset_label = setup$dataset_label,
      sample_ids = rownames(setup$X_truth),
      feature_names = colnames(setup$X_truth),
      layer_vec = setup$layer_vec,
      modalities = setup$modalities,
      unsupported_methods = setup$unsupported_methods,
      fits = fits
    )
  }

  saveRDS(ground_truth, GROUND_TRUTH_FILE)
  write.csv(dplyr::bind_rows(status_rows), GROUND_TRUTH_STATUS_FILE, row.names = FALSE)

  cat("Saved:", GROUND_TRUTH_FILE, "\n")
  cat("Saved:", GROUND_TRUTH_STATUS_FILE, "\n")
  invisible(ground_truth)
}

# -----------------------------------------------------------------------------
# TASK MAPPING
# -----------------------------------------------------------------------------

condition_table <- data.frame(
  condition_id = 1:6,
  mechanism = c("MCAR", "MAR", "MNAR", "MCAR", "MAR", "MNAR"),
  target_rate = c(0.2, 0.2, 0.2, 0.7, 0.7, 0.7),
  stringsAsFactors = FALSE
)

map_task <- function(task_id) {
  total_tasks <- nrow(condition_table) * SEMISYNTH_REPETITIONS

  if (is.na(task_id) || task_id < 1L || task_id > total_tasks) {
    stop("SEMISYNTH_TASK_ID must be between 1 and ", total_tasks, ".")
  }

  condition_index <- ((task_id - 1L) %/% SEMISYNTH_REPETITIONS) + 1L
  repetition <- ((task_id - 1L) %% SEMISYNTH_REPETITIONS) + 1L

  list(
    task_id = task_id,
    repetition = repetition,
    condition_id = condition_table$condition_id[condition_index],
    mechanism = condition_table$mechanism[condition_index],
    target_rate = condition_table$target_rate[condition_index]
  )
}

# -----------------------------------------------------------------------------
# RUN ONE OF THE 600 TASKS
# -----------------------------------------------------------------------------

run_task <- function(task_id) {
  task <- map_task(task_id)
  task_tag <- sprintf("%04d", task_id)

  result_rds <- file.path(TASK_DIR, paste0("task_", task_tag, "_results.rds"))
  result_csv <- file.path(TASK_DIR, paste0("task_", task_tag, "_results.csv"))
  failure_csv <- file.path(TASK_DIR, paste0("task_", task_tag, "_failures.csv"))
  metadata_rds <- file.path(TASK_DIR, paste0("task_", task_tag, "_metadata.rds"))
  checkpoint_rds <- file.path(TASK_DIR, paste0("task_", task_tag, "_checkpoint.rds"))

  if (file.exists(result_rds) && !RETRY_FAILED) {
    cat("Task", task_id, "already has a completed result and RETRY_FAILED=0.\n")
    return(invisible(readRDS(result_rds)))
  }

  if (!file.exists(GROUND_TRUTH_FILE)) {
    stop(
      "Ground-truth cache does not exist: ",
      GROUND_TRUTH_FILE,
      "\nRun SEMISYNTH_MODE=ground_truth first."
    )
  }

  ground_truth <- readRDS(GROUND_TRUTH_FILE)

  if (!identical(ground_truth$methods, METHODS)) {
    stop(
      "The methods in the ground-truth cache do not match SEMISYNTH_METHODS. ",
      "Rebuild the ground-truth cache."
    )
  }
  if (!identical(as.integer(ground_truth$n_factors), SEMISYNTH_N_FACTORS)) {
    stop(
      "The ground-truth factor count does not match SEMISYNTH_N_FACTORS. ",
      "Rebuild the ground-truth cache."
    )
  }

  cat("============================================================\n")
  cat("SEMI-SYNTHETIC TASK", task_id, "OF 600\n")
  cat("Condition:", task$mechanism, task$target_rate, "\n")
  cat("Repetition:", task$repetition, "of", SEMISYNTH_REPETITIONS, "\n")
  cat("Methods:", paste(METHODS, collapse = ", "), "\n")
  cat("missForest ntree:", MISSFOREST_NTREE, "\n")
  cat("missForest maxiter:", MISSFOREST_MAXITER, "\n")
  cat("Save simulated data:", SAVE_SEMISYNTH_DATA, "\n")
  cat("============================================================\n")

  # Resume from a previously saved checkpoint, if one exists.
  checkpoint_df <- NULL

  if (file.exists(checkpoint_rds)) {
    checkpoint_df <- readRDS(checkpoint_rds)

    cat(
      "Resuming task", task_id,
      "from checkpoint with", nrow(checkpoint_df),
      "completed rows.\n"
    )
  }

  rows <- if (!is.null(checkpoint_df) && nrow(checkpoint_df) > 0L) {
    lapply(seq_len(nrow(checkpoint_df)), function(i) checkpoint_df[i, , drop = FALSE])
  } else {
    list()
  }

  row_index <- length(rows)

  for (dataset_index in seq_along(c("adni", "tcga"))) {
    dataset <- c("adni", "tcga")[dataset_index]
    if (!dataset %in% strsplit(Sys.getenv("BENCH_DATASET"), ",", fixed=TRUE)[[1]]) next
    setup <- prepare_dataset(dataset)
    baseline_dataset <- ground_truth$datasets[[dataset]]

    if (is.null(baseline_dataset)) {
      stop("Ground-truth cache is missing dataset: ", dataset)
    }

    if (!identical(rownames(setup$X_truth), baseline_dataset$sample_ids)) {
      stop(
        setup$dataset_label,
        " sample IDs changed after the ground-truth cache was created."
      )
    }
    if (!identical(colnames(setup$X_truth), baseline_dataset$feature_names)) {
      stop(
        setup$dataset_label,
        " feature columns changed after the ground-truth cache was created."
      )
    }

    for (modality_index in seq_along(setup$modalities)) {
      modality <- setup$modalities[modality_index]
      block_cols <- names(setup$layer_vec)[setup$layer_vec == modality]

            # Do not repeat a dataset/modality block that is already complete
      # in the checkpoint.
      if (!is.null(checkpoint_df) && nrow(checkpoint_df) > 0L) {

        completed_models <- unique(
          checkpoint_df$model[
            checkpoint_df$dataset == setup$dataset_label &
            checkpoint_df$dropped_modality == modality
          ]
        )

        if (all(METHODS %in% completed_models)) {
          cat(
            "Checkpoint: skipping completed block:",
            setup$dataset_label, "/", modality, "\n"
          )
          next
        }
      }

      scenario_seed <- as.integer(
        1000000L +
          task_id * 100L +
          (if (identical(Sys.getenv("TCGA_INTNMF_REPAIR"), "1")) 1L else dataset_index) * 10L +
          modality_index
      )

      mar_score <- if (task$mechanism == "MAR") {
        make_mar_score(setup$X_truth, setup$layer_vec, modality)
      } else {
        NULL
      }

      mnar_score <- if (task$mechanism == "MNAR") {
        make_mnar_score(setup$X_truth, setup$layer_vec, modality)
      } else {
        NULL
      }

      simulated <- make_missing(
        X = setup$X_truth,
        block_cols = block_cols,
        mechanism = task$mechanism,
        target_rate = task$target_rate,
        observed_score = mar_score,
        latent_score = mnar_score,
        beta = SEMISYNTH_BETA,
        seed = scenario_seed
      )

      imputation_started <- proc.time()[["elapsed"]]
      imputation <- tryCatch(
        {
          set.seed(scenario_seed + 500000L)
          fitted <- missForest::missForest(
            xmis = simulated$X_miss,
            maxiter = MISSFOREST_MAXITER,
            ntree = MISSFOREST_NTREE,
            verbose = FALSE
          )

          list(
            status = "success",
            error_message = NA_character_,
            elapsed_seconds = proc.time()[["elapsed"]] - imputation_started,
            X = as.matrix(fitted$ximp),
            OOBerror = fitted$OOBerror
          )
        },
        error = function(e) {
          list(
            status = "failed",
            error_message = conditionMessage(e),
            elapsed_seconds = proc.time()[["elapsed"]] - imputation_started,
            X = NULL,
            OOBerror = NA
          )
        }
      )

      if (SAVE_SEMISYNTH_DATA) {
        data_file <- file.path(
          DATA_DIR,
          paste0(
            "task_", task_tag, "_",
            dataset, "_",
            gsub("[^A-Za-z0-9_-]+", "_", modality),
            ".rds"
          )
        )
        saveRDS(
          list(
            task = task,
            dataset = setup$dataset_label,
            modality = modality,
            simulated = simulated,
            imputation = imputation
          ),
          data_file
        )
      }

      for (method_index in seq_along(METHODS)) {
        model <- METHODS[method_index]
        baseline <- baseline_dataset$fits[[model]]

        if (is.null(baseline)) {
          baseline <- list(
            status = "baseline_failed",
            error_message = "The ground-truth cache contains no fit for this model.",
            elapsed_seconds = NA_real_,
            Z = NULL,
            W = NULL
          )
        }

        fitted <- if (model %in% setup$unsupported_methods) {
          list(
            status = "unsupported",
            error_message = paste0(
              model,
              " is a two-block wrapper and is not run on TCGA with ",
              length(setup$modalities),
              " modalities."
            ),
            elapsed_seconds = 0,
            Z = NULL,
            W = NULL
          )
        } else if (!identical(imputation$status, "success")) {
          list(
            status = "imputation_failed",
            error_message = imputation$error_message,
            elapsed_seconds = 0,
            Z = NULL,
            W = NULL
          )
        } else if (!identical(baseline$status, "success")) {
          list(
            status = "not_run",
            error_message = "Skipped because the complete-data baseline failed.",
            elapsed_seconds = 0,
            Z = NULL,
            W = NULL
          )
        } else {
          fit_latent(
            setup = setup,
            model = model,
            X = imputation$X,
            seed = scenario_seed + method_index,
            run_id = task_id,
            style = "SIM_MISSFOREST"
          )
        }

        comparison <- compare_fit(
          baseline = baseline,
          fitted = fitted,
          layer_vec = setup$layer_vec
        )

        row <- list(
          task_id = task_id,
          condition_id = task$condition_id,
          repetition = task$repetition,
          dataset = setup$dataset_label,
          dropped_modality = modality,
          mechanism = task$mechanism,
          target_rate = task$target_rate,
          realized_rate = simulated$realized_rate,
          n_samples = nrow(setup$X_truth),
          n_missing_samples = sum(simulated$missing),
          n_features = ncol(setup$X_truth),
          model = model,
          n_factors = SEMISYNTH_N_FACTORS,
          beta = SEMISYNTH_BETA,
          imputation = "missForest",
          missforest_ntree = MISSFOREST_NTREE,
          missforest_maxiter = MISSFOREST_MAXITER,
          imputation_status = imputation$status,
          imputation_error = imputation$error_message,
          imputation_elapsed_seconds = imputation$elapsed_seconds,
          fit_status = fitted$status,
          fit_error = fitted$error_message,
          fit_elapsed_seconds = fitted$elapsed_seconds,
          status = comparison$status,
          error_message = comparison$error_message,
          n_attempted = 1L,
          n_successful = as.integer(comparison$status == "success"),
          n_failed = as.integer(comparison$status != "success"),
          seed = scenario_seed
        )

        metric_names <- setdiff(
          names(comparison),
          c("status", "error_message")
        )
        for (metric_name in metric_names) {
          row[[metric_name]] <- comparison[[metric_name]]
        }

        row_index <- row_index + 1L
        rows[[row_index]] <- tibble::as_tibble(row)
      }

      saveRDS(dplyr::bind_rows(rows), checkpoint_rds)
    }
  }

  results <- dplyr::bind_rows(rows)

  saveRDS(results, result_rds)
  write.csv(results, result_csv, row.names = FALSE)
  write.csv(
    results[results$status != "success", , drop = FALSE],
    failure_csv,
    row.names = FALSE
  )

  saveRDS(
    list(
      task = task,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      settings = list(
        methods = METHODS,
        repetitions = SEMISYNTH_REPETITIONS,
        n_factors = SEMISYNTH_N_FACTORS,
        beta = SEMISYNTH_BETA,
        missforest_ntree = MISSFOREST_NTREE,
        missforest_maxiter = MISSFOREST_MAXITER,
        save_semisynth_data = SAVE_SEMISYNTH_DATA,
        retry_failed = RETRY_FAILED
      )
    ),
    metadata_rds
  )

  if (file.exists(checkpoint_rds)) unlink(checkpoint_rds)

  cat("============================================================\n")
  cat("TASK", task_id, "COMPLETE\n")
  cat("Results RDS:", result_rds, "\n")
  cat("Results CSV:", result_csv, "\n")
  cat("Failures CSV:", failure_csv, "\n")
  cat("Successful rows:", sum(results$status == "success"), "\n")
  cat("Failed rows:", sum(results$status != "success"), "\n")
  cat("============================================================\n")

  invisible(results)
}

# -----------------------------------------------------------------------------
# COMBINE COMPLETED TASKS
# -----------------------------------------------------------------------------

combine_tasks <- function() {
  files <- list.files(
    TASK_DIR,
    pattern = "^task_[0-9]{4}_results\\.rds$",
    full.names = TRUE
  )

  if (!length(files)) {
    stop("No completed task result RDS files were found in: ", TASK_DIR)
  }

  results <- dplyr::bind_rows(lapply(files, readRDS))

  combined_rds <- file.path(COMBINED_DIR, "semi_synthetic_ADNI_TCGA_all_tasks.rds")
  combined_csv <- file.path(COMBINED_DIR, "semi_synthetic_ADNI_TCGA_all_tasks.csv")
  counts_csv <- file.path(COMBINED_DIR, "semi_synthetic_ADNI_TCGA_run_counts.csv")
  failures_csv <- file.path(COMBINED_DIR, "semi_synthetic_ADNI_TCGA_failures.csv")
  coverage_csv <- file.path(COMBINED_DIR, "semi_synthetic_ADNI_TCGA_task_coverage.csv")

  saveRDS(results, combined_rds)
  write.csv(results, combined_csv, row.names = FALSE)

  counts <- results %>%
    dplyr::group_by(dataset, dropped_modality, mechanism, target_rate, model) %>%
    dplyr::summarise(
      n_attempted = dplyr::n(),
      n_successful = sum(status == "success"),
      n_failed = sum(status != "success"),
      success_rate = n_successful / n_attempted,
      .groups = "drop"
    )
  write.csv(counts, counts_csv, row.names = FALSE)

  failures <- results[results$status != "success", , drop = FALSE]
  write.csv(failures, failures_csv, row.names = FALSE)

  expected_tasks <- seq_len(nrow(condition_table) * SEMISYNTH_REPETITIONS)
  completed_tasks <- sort(unique(results$task_id))
  coverage <- data.frame(task_id = expected_tasks, completed = expected_tasks %in% completed_tasks)
  write.csv(coverage, coverage_csv, row.names = FALSE)

  cat("============================================================\n")
  cat("COMBINED SEMI-SYNTHETIC RESULTS\n")
  cat("Completed tasks:", length(completed_tasks), "of", length(expected_tasks), "\n")
  cat("Combined RDS:", combined_rds, "\n")
  cat("Combined CSV:", combined_csv, "\n")
  cat("Run counts:", counts_csv, "\n")
  cat("Failures:", failures_csv, "\n")
  cat("Task coverage:", coverage_csv, "\n")
  cat("============================================================\n")

  invisible(results)
}

# -----------------------------------------------------------------------------
# ENTRY POINT
# -----------------------------------------------------------------------------

