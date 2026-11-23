#!/usr/bin/env Rscript
# Reproduce one selected tICA scenario for Figure 7 (usage: <dataset> <modality>; run from the
# archive root). Function definitions are parsed from the historical analysis scripts, which
# are never sourced or executed. The refit is Procrustes-aligned to the reference scores and
# projected onto a seeded reference UMAP
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
root <- normalizePath(".")
out <- file.path(root, "figures_revised")
.libPaths(c(file.path(Sys.getenv("HOME"), "R", "4.5"), .libPaths()))
# Fail before imputation if the fitting environment cannot load its dependencies
stopifnot(
  requireNamespace("tensorICA", quietly = TRUE),
  requireNamespace("missForest", quietly = TRUE),
  requireNamespace("FactoMineR", quietly = TRUE)
)

# Selected repetitions for this dataset and modality
ds <- args[1]
mod <- args[2]
selection <- read.csv(file.path(out, "tables", "semisynthetic_UMAP_selection.csv"))
selection <- selection[selection$plot_dataset == ds & selection$dropped_modality == mod, ]
stopifnot(nrow(selection) == 2)
# Evaluate only named function definitions, avoiding bootstrap execution,
# package installation, data preparation, or writes to original directories
e <- new.env(parent = .GlobalEnv)

# Load selected functions without executing source scripts
assign_functions <- function(path, env) {
  txt <- readLines(path, warn = FALSE)
  # final.R contains quoted historical blocks marked by triple single quotes
  markers <- which(trimws(txt) == "'''")
  if (length(markers)) {
    for (i in seq(1, length(markers), by = 2)) txt[markers[i]:markers[i + 1]] <- ""
  }
  for (expr in parse(text = txt)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.call(expr[[3]]) &&
      identical(expr[[3]][[1]], as.name("function"))) eval(expr, env)
  }
}

assign_functions(file.path(root, if (ds == "ADNI") "final.R" else "final_tcga.R"), e)
assign_functions(
  file.path(
    root, if (ds == "ADNI") "semi_synthetic_adni_recover_431_440.R" else "semi_synthetic_tcga_only.R"
  ),
  e
)

e$SEMISYNTH_N_FACTORS <- 3L

# Truth matrix and reference tICA fit
vars <- readRDS(file.path(
  root, "variables",
  if (ds == "ADNI") "analysis_variables.rds" else "tcga_bootstrap_analysis_variables.rds"
))
archive <- if (ds == "ADNI") {
  "results_semi_synthetic_adni/source_archives/results_semi_synthetic_final_v1"
} else {
  "results_semi_synthetic_tcga/source_archives/results_semi_synthetic_TCGA_only_v1"
}
ground <- readRDS(file.path(root, archive, "ground_truth/ground_truth_fits.rds"))
baseline <- ground$datasets[[tolower(ds)]]
stopifnot("tica" %in% ground$methods, baseline$fits$tica$status == "success")
X <- vars$X_stand[baseline$sample_ids, baseline$feature_names, drop = FALSE]
lv <- baseline$layer_vec[colnames(X)]
stopifnot(
  !anyNA(X), all(is.finite(X)), !anyNA(lv), identical(colnames(X), baseline$feature_names),
  identical(rownames(X), baseline$sample_ids)
)

setup <- list(
  latent_runner = e$make_latent_runner(e$run_model), layer_vec = lv, diag_labels = vars$diag_labels,
  missing_labels = setNames(rep("CC", nrow(X)), rownames(X)), clinical_vec = vars$mmse_vec
)

Z0 <- baseline$fits$tica$Z[rownames(X), , drop = FALSE]
stopifnot(all(is.finite(Z0)), identical(rownames(Z0), rownames(X)))

# Reference UMAP of the complete-data scores
prefix <- paste(ds, mod, sep = "_")
ref_file <- file.path(out, "cache", paste0(prefix, "_reference_UMAP.rds"))
set.seed(20261018)
reference <- uwot::umap(
  Z0, n_neighbors = 15, min_dist = .15, metric = "euclidean", n_epochs = 500, scale = FALSE,
  n_threads = 1, n_sgd_threads = 1, ret_model = TRUE, verbose = FALSE
)
saveRDS(reference, ref_file)
write.csv(
  data.frame(patient_id = rownames(Z0), Z0),
  file.path(out, "cache", paste0(prefix, "_reference_Z.csv")), row.names = FALSE
)
write.csv(
  data.frame(
    patient_id = rownames(Z0), UMAP1 = reference$embedding[, 1], UMAP2 = reference$embedding[, 2]
  ),
  file.path(out, "cache", paste0(prefix, "_reference_UMAP.csv")), row.names = FALSE
)

records <- list()

# Refit tICA for each selected rate and project it onto the reference UMAP
for (i in seq_len(nrow(selection))) {
  r <- selection[i, ]
  rate <- r$target_rate
  tag <- paste0(prefix, "_", round(rate * 100))
  cache <- file.path(out, "cache", paste0(tag, "_fit.rds"))

  if (file.exists(cache)) {
    refit <- readRDS(cache)
  } else {
    seed <- as.integer(r$seed)
    imputed_file <- file.path(out, "cache", paste0(tag, "_imputed.rds"))
    if (file.exists(imputed_file)) {
      imputed <- readRDS(imputed_file)
    } else {
      simulated <- e$make_missing(
        X = X, block_cols = names(lv)[lv == mod], mechanism = "MCAR", target_rate = rate,
        observed_score = NULL, latent_score = NULL, beta = 2, seed = seed
      )
      cat("Imputing", tag, "seed", seed, "\n")
      flush.console()
      set.seed(seed + 500000L)
      imp <- missForest::missForest(
        xmis = simulated$X_miss, maxiter = 10, ntree = 100, verbose = FALSE
      )
      imputed <- list(X = as.matrix(imp$ximp), realized_rate = simulated$realized_rate)
      saveRDS(imputed, imputed_file)
    }

    fit <- e$fit_latent(
      setup = setup, model = "tica", X = imputed$X,
      seed = seed + match("tica", ground$methods), run_id = as.integer(r$task_id), style = "SIM"
    )
    if (fit$status != "success") stop("Fit failed: ", fit$error_message)
    refit <- list(
      Z = fit$Z, W = fit$W, seed = seed, realized_rate = imputed$realized_rate,
      reference_file = file.path(archive, "ground_truth/ground_truth_fits.rds"),
      features = colnames(X), ids = rownames(X)
    )
    saveRDS(refit, cache)
  }

  Z <- refit$Z[rownames(Z0), , drop = FALSE]
  stopifnot(identical(rownames(Z), rownames(Z0)), all(is.finite(Z)))
  rv <- FactoMineR::coeffRV(Z0, Z)$rv
  # Orthogonal Procrustes rotation plus one global scale removes arbitrary factor signs,
  # rotations, and amplitude while preserving within-fit distances up to scale
  A <- scale(Z, center = TRUE, scale = FALSE)
  B <- scale(Z0, center = TRUE, scale = FALSE)
  s <- svd(crossprod(A, B))
  rotation <- s$u %*% t(s$v)
  multiplier <- sum(s$d) / sum(A * A)
  aligned <- sweep(multiplier * A %*% rotation, 2, colMeans(Z0), "+")

  set.seed(20261018)
  projected <- uwot::umap_transform(
    aligned, reference, n_threads = 1, n_sgd_threads = 1, verbose = FALSE
  )
  write.csv(
    data.frame(patient_id = rownames(Z), Z), file.path(out, "cache", paste0(tag, "_Z.csv")),
    row.names = FALSE
  )
  write.csv(
    data.frame(patient_id = rownames(Z), UMAP1 = projected[, 1], UMAP2 = projected[, 2]),
    file.path(out, "cache", paste0(tag, "_UMAP.csv")), row.names = FALSE
  )

  records[[i]] <- data.frame(
    dataset = ds, modality = mod, model = "tica", mechanism = "MCAR", repetition = r$repetition,
    task_id = r$task_id, target_rate = rate, realized_rate = refit$realized_rate, seed = r$seed,
    saved_rv = r$rv, visualization_rv = rv, rv_difference = rv - r$rv, n = nrow(Z),
    features = ncol(X), global_scale = multiplier,
    reference = file.path(archive, "ground_truth/ground_truth_fits.rds")
  )
  cat(tag, "RV saved", r$rv, "reproduced", rv, "\n")
  flush.console()
}

write.csv(
  do.call(rbind, records), file.path(out, "tables", paste0(prefix, "_visualization_manifest.csv")),
  row.names = FALSE
)
