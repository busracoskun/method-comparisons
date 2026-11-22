# Evaluate saved ADNI full-cohort model fits.
# ==============================================================================
# build_full_run_metrics_from_FullRuns.R
# UPDATED VERSION
#
# PURPOSE
# -------
# Re-score the already-saved ADNI full-run model objects under:
#   <configured output>/adni/{fits,fits_python}
#
# This does NOT refit models.
#
# IMPORTANT FIXES IN THIS VERSION
# --------------------------------
# 1. CC_SF missingness metrics:
#    Complete-case data can contain only one missingness label ("CC").
#    In that case, clustering-vs-missingness metrics are undefined, so they
#    are recorded as NA instead of causing the entire saved model to fail.
#
# 2. MOFA/CMF feature snapshots:
#    Diagnosis clustering and MMSE R2 depend only on saved latent scores Z.
#    They are therefore still scored even if the saved feature set no longer
#    exactly matches the current analysis_variables.rds feature set.
#
#    Reconstruction PVE is computed only when the saved feature set can be
#    matched exactly and safely. Otherwise PVE is NA, with a note explaining why.
#
# 3. Missing/non-overlapping labels:
#    Rows with missing labels are dropped only for that particular clustering
#    metric rather than causing the whole model to fail.
#
# 4. Duplicated latent-component names:
#    Some saved Block PLS / RGCCA / sGCCA objects contain repeated names
#    such as comp1. These are made unique only for the MMSE regression;
#    the latent-score values themselves are unchanged.
#
# OUTPUTS
# -------
# <configured output>/adni/full_rescored/full_run_metrics_from_FullRuns_raw.csv
# <configured output>/adni/full_rescored/full_run_metrics_from_FullRuns.csv
# <configured output>/adni/full_rescored/full_run_metric_errors.csv
# <configured output>/adni/full_rescored/full_run_metric_notes.csv
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------------------------

source(file.path(Sys.getenv("CLEAN_ROOT"), "code","metrics", "runtime.R"))

required_packages <- c(
  "dplyr",
  "clusterCrit",
  "mclust",
  "matrixStats"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv("BENCH_WORK_DIR")
FULL_DIR <- c(file.path(PROJECT_DIR, "fits"), file.path(PROJECT_DIR, "fits_python"))
FULL_DIR <- FULL_DIR[dir.exists(FULL_DIR)]
OUT_DIR <- file.path(PROJECT_DIR, "full_rescored")

ANALYSIS_VARS_FILE <- file.path(
  PROJECT_DIR,
  "variables",
  "analysis_variables.rds"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!length(FULL_DIR)) {
  stop("FullRuns directory not found: ", FULL_DIR)
}

if (!file.exists(ANALYSIS_VARS_FILE)) {
  stop("analysis_variables.rds not found: ", ANALYSIS_VARS_FILE)
}

# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

compute_pve <- function(latent_matrix, X_original, W_loadings) {

  if (is.null(X_original) || is.null(W_loadings)) {
    return(NA_real_)
  }

  latent_matrix <- as.matrix(latent_matrix)
  X_original <- as.matrix(X_original)
  W_loadings <- as.matrix(W_loadings)

  if (
    is.null(rownames(latent_matrix)) ||
    is.null(rownames(X_original)) ||
    is.null(rownames(W_loadings))
  ) {
    return(NA_real_)
  }

  ids <- intersect(rownames(latent_matrix), rownames(X_original))

  feats <- intersect(colnames(X_original), rownames(W_loadings))

  if (!length(ids) || !length(feats)) {
    return(NA_real_)
  }

  Z <- latent_matrix[
    ids,
    ,
    drop = FALSE
  ]

  X_use <- X_original[
    ids,
    feats,
    drop = FALSE
  ]

  W_use <- W_loadings[
    feats,
    ,
    drop = FALSE
  ]

  if (any(!is.finite(Z)) || any(!is.finite(W_use))) {
    return(NA_real_)
  }

  Xhat <- Z %*% t(W_use)

  mu_orig <- colMeans(X_use, na.rm = TRUE)

  Xc <- sweep(
    X_use,
    2,
    mu_orig,
    FUN = "-"
  )

  mu_hat <- colMeans(Xhat, na.rm = TRUE)

  Xhatc <- sweep(
    Xhat,
    2,
    mu_hat,
    FUN = "-"
  )

  var_X <- sum(matrixStats::colVars(Xc, na.rm = TRUE))

  var_Xhat <- sum(matrixStats::colVars(Xhatc, na.rm = TRUE))

  if (!is.finite(var_X) || var_X <= 0 || !is.finite(var_Xhat)) {
    return(NA_real_)
  }

  var_Xhat / var_X
}


score_clusters <- function(latent_matrix, true_labels, pve_value = NA_real_, seed = 100L) {

  latent_matrix <- as.matrix(latent_matrix)

  if (is.null(rownames(latent_matrix))) {
    stop("Latent matrix has no row names.")
  }

  if (is.null(names(true_labels))) {
    stop("True-label vector has no sample names.")
  }

  ids <- intersect(rownames(latent_matrix), names(true_labels))

  if (!length(ids)) {
    stop("No overlapping sample IDs for clustering.")
  }

  Z <- latent_matrix[
    ids,
    ,
    drop = FALSE
  ]

  y <- true_labels[ids]

  keep <- !is.na(y)

  Z <- Z[
    keep,
    ,
    drop = FALSE
  ]

  y <- y[keep]

  if (!nrow(Z)) {
    stop("No non-missing labels available for clustering.")
  }

  if (any(!is.finite(Z))) {
    stop("Latent matrix contains non-finite values.")
  }

  y <- droplevels(factor(y))

  k <- nlevels(y)

  # IMPORTANT:
  # This is expected for CC_SF missingness labels, because all complete-case
  # subjects can have the same missingness class ("CC").
  # The clustering metrics are undefined, but the saved model itself is valid.
  if (k < 2L) {
    return(c(rand = NA_real_, ari = NA_real_, chi = NA_real_, db = NA_real_, pve = pve_value))
  }

  if (nrow(Z) <= k) {
    return(c(rand = NA_real_, ari = NA_real_, chi = NA_real_, db = NA_real_, pve = pve_value))
  }

  set.seed(seed)

  pred <- kmeans(Z, centers = k, nstart = 25)$cluster

  rand <- clusterCrit::extCriteria(
    as.integer(y),
    pred,
    "Rand"
  )$rand

  ari <- mclust::adjustedRandIndex(y, pred)

  internal <- clusterCrit::intCriteria(
    as.matrix(Z),
    pred,
    c(
      "Calinski_Harabasz",
      "Davies_Bouldin"
    )
  )

  c(
    rand = rand,
    ari = ari,
    chi = internal$calinski_harabasz,
    db = internal$davies_bouldin,
    pve = pve_value
  )
}


score_regression <- function(latent_matrix, mmse_vec) {

  latent_matrix <- as.matrix(latent_matrix)

  if (is.null(rownames(latent_matrix)) || is.null(names(mmse_vec))) {
    return(NA_real_)
  }

  ids <- intersect(rownames(latent_matrix), names(mmse_vec))

  if (!length(ids)) {
    return(NA_real_)
  }

  Z <- latent_matrix[
    ids,
    ,
    drop = FALSE
  ]

  # Some saved methods (notably Block PLS / RGCCA / sGCCA) can store
  # duplicated latent-component names such as "comp1" more than once.
  # lm(MMSE ~ .) requires unique predictor names. Renaming the columns is
  # harmless because the numerical latent scores are unchanged.
  if (is.null(colnames(Z))) {
    colnames(Z) <- paste0(
      "latent_",
      seq_len(ncol(Z))
    )
  } else {
    colnames(Z) <- make.unique(make.names(colnames(Z), unique = FALSE))
  }

  df <- data.frame(MMSE = mmse_vec[ids], Z, check.names = FALSE)

  df <- df[
    complete.cases(df),
    ,
    drop = FALSE
  ]

  if (nrow(df) <= ncol(Z) + 1L) {
    return(NA_real_)
  }

  fit <- lm(MMSE ~ ., data = df)

  summary(fit)$r.squared
}

# ------------------------------------------------------------------------------
# 3. LOAD ANALYSIS VARIABLES
# ------------------------------------------------------------------------------

av <- readRDS(ANALYSIS_VARS_FILE)

required_av <- c(
  "X_na",
  "X_cc_sf",
  "diag_labels",
  "missingness_label",
  "mmse_vec",
  "imp_list_file"
)

missing_av <- setdiff(required_av, names(av))

if (length(missing_av)) {
  stop(
    "analysis_variables.rds is missing: ",
    paste(missing_av, collapse = ", ")
  )
}

if (!file.exists(av$imp_list_file)) {
  stop(
    "Saved MI imputation file not found: ",
    av$imp_list_file
  )
}

imp_list <- readRDS(av$imp_list_file)

# ------------------------------------------------------------------------------
# 4. FIND SAVED FULL-RUN RDS OBJECTS
# ------------------------------------------------------------------------------

files <- list.files(
  FULL_DIR,
  pattern = "_(NA|CC_SF|MI_SF)_0(_imp[0-9]+)?\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(files)) {
  stop(
    "No full-run .rds model objects found under: ",
    FULL_DIR
  )
}

cat(
  "Found ",
  length(files),
  " full-run RDS files.\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 5. FILE-NAME PARSER
# ------------------------------------------------------------------------------

parse_filename <- function(path) {

  f <- basename(path)

  m <- regexec(
    "^(.+)_(NA|CC_SF|MI_SF)_0(?:_imp([0-9]+))?\\.rds$",
    f,
    perl = TRUE
  )

  z <- regmatches(f, m)[[1]]

  if (!length(z)) {
    stop(
      "Could not parse full-run filename: ",
      f
    )
  }

  list(
    model = tolower(z[2]),
    style = toupper(z[3]),
    imp = if (length(z) >= 4L && nzchar(z[4])) {
      as.integer(z[4])
    } else {
      NA_integer_
    }
  )
}

# ------------------------------------------------------------------------------
# 6. SCORE ONE SAVED FULL RUN
# ------------------------------------------------------------------------------

score_one_file <- function(path) {

  meta <- parse_filename(path)

  obj <- readRDS(path)

  if (!is.list(obj)) {
    stop("Saved object is not a list.")
  }

  if (is.null(obj$Z)) {
    stop("Saved object does not contain Z.")
  }

  model <- tolower(as.character(obj$model %||% meta$model))

  style <- toupper(as.character(obj$style %||% meta$style))

  imp <- obj$imp %||% meta$imp

  if (length(imp) == 0L || is.null(imp) || is.na(imp)) {
    imp <- NA_integer_
  } else {
    imp <- as.integer(imp)
  }

  Z <- as.matrix(obj$Z)

  if (is.null(rownames(Z))) {

    if (!is.null(obj$rows) && length(obj$rows) == nrow(Z)) {
      rownames(Z) <- as.character(obj$rows)
    } else {
      stop(
        "Z has no row names and saved rows cannot restore them."
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Choose the data object for PVE only.
  #
  # IMPORTANT:
  # Diagnosis ARI / missingness ARI / R2 below do NOT depend on this feature
  # matrix. Therefore a feature mismatch must never discard those metrics.
  # ---------------------------------------------------------------------------

  X_source <- switch(
    style,

    "NA" = av$X_na,

    "CC_SF" = av$X_cc_sf,

    "MI_SF" = {
      if (is.na(imp) || imp < 1L || imp > length(imp_list)) {
        stop(
          "MI_SF file has invalid imputation index."
        )
      }

      imp_list[[imp]]
    },

    stop(
      "Unsupported style in FullRuns: ",
      style
    )
  )

  W <- obj$W

  if (!is.null(W)) {

    W <- as.matrix(W)

    if (is.null(rownames(W)) && !is.null(obj$cols) && length(obj$cols) == nrow(W)) {
      rownames(W) <- as.character(obj$cols)
    }
  }

  pve_value <- NA_real_
  pve_status <- "unavailable"

  saved_cols <- obj$cols

  if (is.null(saved_cols) || !length(saved_cols)) {
    if (!is.null(W) && !is.null(rownames(W))) {
      saved_cols <- rownames(W)
    }
  }

  rows_match <- all(rownames(Z) %in% rownames(X_source))

  cols_match <- (
    !is.null(saved_cols) &&
    length(saved_cols) > 0L &&
    all(saved_cols %in% colnames(X_source))
  )

  W_match <- (
    !is.null(W) &&
    !is.null(rownames(W)) &&
    !is.null(saved_cols) &&
    length(saved_cols) > 0L &&
    all(saved_cols %in% rownames(W))
  )

  if (rows_match && cols_match && W_match) {

    X_for_pve <- X_source[
      rownames(Z),
      saved_cols,
      drop = FALSE
    ]

    W_for_pve <- W[
      saved_cols,
      ,
      drop = FALSE
    ]

    pve_value <- compute_pve(latent_matrix = Z, X_original = X_for_pve, W_loadings = W_for_pve)

    if (is.finite(pve_value)) {
      pve_status <- "success"
    } else {
      pve_status <- "calculation_failed"
    }

  } else {

    reasons <- character()

    if (!rows_match) {
      reasons <- c(
        reasons,
        "saved sample IDs do not exactly match current analysis data"
      )
    }

    if (!cols_match) {
      reasons <- c(
        reasons,
        "saved feature set does not exactly match current analysis data"
      )
    }

    if (!W_match) {
      reasons <- c(
        reasons,
        "saved W/loadings cannot be aligned exactly to saved feature set"
      )
    }

    pve_status <- paste(
      reasons,
      collapse = "; "
    )
  }

  # ---------------------------------------------------------------------------
  # Z-BASED METRICS
  # ---------------------------------------------------------------------------

  diag_score <- score_clusters(
    latent_matrix = Z,
    true_labels = av$diag_labels,
    pve_value = pve_value,
    seed = 100L
  )

  miss_score <- score_clusters(
    latent_matrix = Z,
    true_labels = av$missingness_label,
    pve_value = pve_value,
    seed = 200L
  )

  R2_value <- score_regression(latent_matrix = Z, mmse_vec = av$mmse_vec)

  # Record how many label groups were actually available.
  diag_ids <- intersect(rownames(Z), names(av$diag_labels))

  miss_ids <- intersect(rownames(Z), names(av$missingness_label))

  n_diag_levels <- length(unique(na.omit(av$diag_labels[diag_ids])))

  n_missingness_levels <- length(unique(na.omit(av$missingness_label[miss_ids])))

  tibble::tibble(
    model = model,
    style = style,
    imp = imp,
    source_file = path,

    n_diag_levels = n_diag_levels,
    n_missingness_levels = n_missingness_levels,

    pve_status = pve_status,

    rand_diag = unname(
      diag_score["rand"]
    ),

    ari_diag = unname(
      diag_score["ari"]
    ),

    chi_diag = unname(
      diag_score["chi"]
    ),

    db_diag = unname(
      diag_score["db"]
    ),

    pve_diag = unname(
      diag_score["pve"]
    ),

    R2 = R2_value,

    rand_miss = unname(
      miss_score["rand"]
    ),

    ari_miss = unname(
      miss_score["ari"]
    ),

    chi_miss = unname(
      miss_score["chi"]
    ),

    db_miss = unname(
      miss_score["db"]
    ),

    pve_miss = unname(
      miss_score["pve"]
    )
  )
}

# ------------------------------------------------------------------------------
# 7. SCORE ALL SAVED OBJECTS SAFELY
# ------------------------------------------------------------------------------

results <- list()
errors <- list()

for (i in seq_along(files)) {

  f <- files[i]

  cat(
    sprintf(
      "[%03d/%03d] %s\n",
      i,
      length(files),
      f
    )
  )

  out <- tryCatch(score_one_file(f), error = function(e) e)

  if (inherits(out, "error")) {

    errors[[length(errors) + 1L]] <- data.frame(
      source_file = f,
      error_message = conditionMessage(out),
      stringsAsFactors = FALSE
    )

    cat(
      "  ERROR: ",
      conditionMessage(out),
      "\n",
      sep = ""
    )

  } else {

    results[[length(results) + 1L]] <- out

    note_text <- character()

    if (out$n_missingness_levels < 2L) {
      note_text <- c(
        note_text,
        "missingness clustering undefined (<2 label levels)"
      )
    }

    if (
      out$pve_status != "success"
    ) {
      note_text <- c(
        note_text,
        paste0(
          "PVE unavailable: ",
          out$pve_status
        )
      )
    }

    cat(
      "  OK: ",
      out$model,
      " / ",
      out$style,
      if (!is.na(out$imp)) {
        paste0(" / imp ", out$imp)
      } else {
        ""
      },
      if (length(note_text)) {
        paste0(
          " [",
          paste(
            note_text,
            collapse = " | "
          ),
          "]"
        )
      } else {
        ""
      },
      "\n",
      sep = ""
    )
  }

  rm(out)
  gc(verbose = FALSE)
}

raw_metrics <- bind_rows(results)

error_table <- if (length(errors)) {
  bind_rows(errors)
} else {
  data.frame(source_file = character(), error_message = character())
}

if (!nrow(raw_metrics)) {
  stop(
    "No full-run metrics could be reconstructed."
  )
}

# ------------------------------------------------------------------------------
# 8. AGGREGATE TO ONE LOGICAL FULL RUN PER MODEL / STYLE
#
# MI_SF is averaged across saved imputations, matching the benchmark workflow.
# NA and CC_SF normally have one saved object each.
# ------------------------------------------------------------------------------

metric_cols <- c(
  "rand_diag",
  "ari_diag",
  "chi_diag",
  "db_diag",
  "pve_diag",
  "R2",
  "rand_miss",
  "ari_miss",
  "chi_miss",
  "db_miss",
  "pve_miss"
)

full_summary <- raw_metrics %>%
  group_by(model, style) %>%
  summarise(
    n_saved_runs = n(),

    n_pve_available = sum(is.finite(pve_diag)),

    n_missingness_metric_available = sum(is.finite(ari_miss)),

    across(all_of(metric_cols), safe_mean),

    source_files = paste(
      basename(source_file),
      collapse = ";"
    ),

    .groups = "drop"
  ) %>%
  arrange(style, model)

# ------------------------------------------------------------------------------
# 9. NOTES TABLE
# ------------------------------------------------------------------------------

notes_table <- raw_metrics %>%
  transmute(
    model,
    style,
    imp,
    source_file,
    n_diag_levels,
    n_missingness_levels,
    pve_status,
    missingness_metric_note = ifelse(
      n_missingness_levels < 2L,
      "Undefined because fewer than two missingness label levels are present",
      ""
    )
  ) %>%
  filter(
    pve_status != "success" |
      missingness_metric_note != ""
  )

# ------------------------------------------------------------------------------
# 10. WRITE
# ------------------------------------------------------------------------------

raw_file <- file.path(
  OUT_DIR,
  "full_run_metrics_from_FullRuns_raw.csv"
)

summary_file <- file.path(
  OUT_DIR,
  "full_run_metrics_from_FullRuns.csv"
)

error_file <- file.path(
  OUT_DIR,
  "full_run_metric_errors.csv"
)

notes_file <- file.path(
  OUT_DIR,
  "full_run_metric_notes.csv"
)

write.csv(
  raw_metrics,
  raw_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  full_summary,
  summary_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  error_table,
  error_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  notes_table,
  notes_file,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------------------------
# 11. REPORT
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FULL-RUN METRIC RECONSTRUCTION COMPLETE\n")
cat("============================================================\n")

cat(
  "Successfully scored saved objects: ",
  nrow(raw_metrics),
  "\n",
  sep = ""
)

cat(
  "True unrecoverable errors: ",
  nrow(error_table),
  "\n",
  sep = ""
)

cat(
  "Objects with undefined missingness clustering: ",
  sum(raw_metrics$n_missingness_levels < 2L),
  "\n",
  sep = ""
)

cat(
  "Objects without comparable PVE: ",
  sum(!is.finite(raw_metrics$pve_diag)),
  "\n",
  sep = ""
)

cat("\nLogical full runs by model/style:\n")

print(
  full_summary %>%
    select(model, style, n_saved_runs, n_pve_available, n_missingness_metric_available),
  n = Inf
)

cat("\nCreated:\n")
cat("  ", raw_file, "\n", sep = "")
cat("  ", summary_file, "\n", sep = "")
cat("  ", error_file, "\n", sep = "")
cat("  ", notes_file, "\n", sep = "")
cat("============================================================\n")