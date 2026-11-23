# Export the latent scores of saved full-cohort GPCCA/NA, MCIA/CC_SF, and GFA/MI_SF fits for
# Figure S3 (usage: <archive_root> <cache_dir>). No model fitting or data preparation

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
cache <- args[2]

entries <- list()
exclusions <- list()
for (dataset in c("adni", "tcga")) {
  vars <- readRDS(file.path(
    root, "results", dataset, "variables",
    if (dataset == "adni") "analysis_variables.rds" else "tcga_bootstrap_analysis_variables.rds"
  ))

  # Availability pattern and diagnosis of each patient
  X <- vars$X_stand
  lv <- vars$layer_vec[colnames(X)]
  modalities <- unique(unname(lv))
  available <- sapply(modalities, function(m) rowSums(!is.na(X[, lv == m, drop = FALSE])) > 0)
  colnames(available) <- modalities
  patterns <- apply(
    available, 1,
    function(x) if (all(x)) "CC" else paste0("missing_", paste(names(x)[!x], collapse = "+"))
  )
  diagnosis <- toupper(as.character(vars$diag_labels[rownames(X)]))
  names(diagnosis) <- names(patterns) <- rownames(X)
  valid <- if (dataset == "adni") c("CN", "MCI", "AD") else c("GBM", "LGG")

  # Export the scores of each saved fit for eligible patients
  for (spec in list(c("gpcca", "NA"), c("mcia", "CC_SF"), c("gfa", "MI_SF"))) {
    model <- spec[1]
    strategy <- spec[2]
    fitfile <- file.path(
      root, "results", dataset, "fits", model,
      paste0(model, "_", strategy, "_0", if (strategy == "MI_SF") "_imp1" else "", ".rds")
    )

    fit <- readRDS(fitfile)
    Z <- as.matrix(fit$Z)
    ids <- rownames(Z)
    if (is.null(ids)) ids <- fit$rows
    stopifnot(
      !anyDuplicated(ids), length(ids) == nrow(Z), all(ids %in% rownames(X)), all(is.finite(Z)),
      ncol(Z) > 1
    )

    keep <- !is.na(diagnosis[ids]) & diagnosis[ids] %in% valid
    exclusions[[length(exclusions) + 1]] <- data.frame(
      dataset = dataset, model = model, strategy = strategy, source_n = length(ids),
      eligible_n = sum(keep), excluded_n = sum(!keep),
      excluded_ids = paste(ids[!keep], collapse = ";")
    )
    ids <- ids[keep]
    Z <- Z[keep, , drop = FALSE]
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))

    folder <- file.path(cache, dataset)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    name <- paste0(model, "_", strategy, "_full.csv")
    write.csv(
      data.frame(
        patient_id = ids, diagnosis = diagnosis[ids], availability = patterns[ids], Z,
        check.names = FALSE
      ),
      file.path(folder, name), row.names = FALSE
    )

    entries[[length(entries) + 1]] <- list(
      dataset = dataset, model = model, strategy = strategy,
      imputation = if (strategy == "MI_SF") 1L else NULL, file = paste0(dataset, "/", name),
      source = sub(paste0(root, "/"), "", fitfile, fixed = TRUE)
    )
  }
}

# Write the export index and coverage
jsonlite::write_json(
  entries, file.path(cache, "observed.json"), pretty = TRUE, auto_unbox = TRUE, null = "null"
)
write.csv(
  do.call(rbind, exclusions), file.path(cache, "observed_export_coverage.csv"), row.names = FALSE
)
