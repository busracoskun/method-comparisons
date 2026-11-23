# Check the 1000 bootstrap inputs per dataset: distinct patients, entirely missing rows and
# columns, and per-modality support and zero-variance features. See README.md

diag_out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
for (part in c("private", "public")) {
  dir.create(file.path(diag_out, part), recursive = TRUE, showWarnings = FALSE)
}

args <- commandArgs(TRUE)
root <- if (length(args)) {
  args[1]
} else {
  Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
}

p <- Sys.getenv("BENCH_ARCHIVE_DIR")
options(stringsAsFactors = FALSE)

all <- list()
for (ds in c("adni", "tcga")) {
  v <- readRDS(file.path(
    p, "variables",
    if (ds == "adni") "analysis_variables.rds" else "tcga_bootstrap_analysis_variables.rds"
  ))

  X <- v$X_stand
  lv <- v$layer_vec[colnames(X)]
  blocks <- split(seq_len(ncol(X)), lv)
  cat(ds, dim(X), "\n")

  # Input properties of each bootstrap replicate
  rows <- lapply(1:1000, function(b) {
    set.seed(100 + b)
    ix <- sample(rownames(X), nrow(X), replace = TRUE)
    Y <- X[ix, , drop = FALSE]
    z <- data.frame(
      dataset = ds, boot = b, n = nrow(Y), p = ncol(Y), unique_samples = length(unique(ix)),
      missing_fraction = mean(is.na(Y)), all_missing_rows = sum(rowSums(!is.na(Y)) == 0),
      all_missing_columns = sum(colSums(!is.na(Y)) == 0),
      nonfinite_nonmissing = sum(!is.finite(Y) & !is.na(Y))
    )

    # Per-modality support in the replicate
    for (nm in names(blocks)) {
      A <- Y[, blocks[[nm]], drop = FALSE]
      obs <- rowSums(!is.na(A)) > 0
      vari <- apply(A, 2, var, na.rm = TRUE)
      z[[paste0(nm, "_observed_rows")]] <- sum(obs)
      z[[paste0(nm, "_unique_observed_rows")]] <- length(unique(ix[obs]))
      z[[paste0(nm, "_zero_variance_columns")]] <- sum(vari == 0, na.rm = TRUE)
      z[[paste0(nm, "_min_feature_observations")]] <- min(colSums(!is.na(A)))
    }
    z
  })
  all[[ds]] <- do.call(rbind, rows)
}

# Write per-replicate checks and the overview
write.csv(do.call(rbind, lapply(names(all), function(ds) {
  x <- all[[ds]]
  write.csv(
    x, file.path(root, "private", paste0(ds, "_bootstrap_input_checks.csv")), row.names = FALSE
  )
  data.frame(
    dataset = ds, n = nrow(x), minimum_unique = min(x$unique_samples),
    maximum_unique = max(x$unique_samples)
  )
})), file.path(root, "public", "input_overview.csv"), row.names = FALSE)
