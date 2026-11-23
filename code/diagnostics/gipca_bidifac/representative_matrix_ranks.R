# Column rank of each centered complete-case block in representative bootstrap samples
# (no fitting, no feature removal). See README.md

diag_out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
for (part in c("private", "public")) {
  dir.create(file.path(diag_out, part), recursive = TRUE, showWarnings = FALSE)
}

root <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
p <- Sys.getenv("BENCH_ARCHIVE_DIR")

out <- list()
for (ds in c("adni", "tcga")) {
  v <- readRDS(file.path(
    p, "variables",
    if (ds == "adni") "analysis_variables.rds" else "tcga_bootstrap_analysis_variables.rds"
  ))

  X <- v$X_stand
  lv <- v$layer_vec[colnames(X)]

  # Ranks of the complete-case block in each representative replicate
  for (b in if (ds == "adni") c(1, 3, 8, 865) else c(1, 4, 6)) {
    set.seed(100 + b)
    ix <- sample(seq_len(nrow(X)), nrow(X), TRUE)
    for (bl in unique(lv)) {
      A <- X[ix, lv == bl, drop = FALSE]
      A <- A[complete.cases(A), , drop = FALSE]
      out[[length(out) + 1]] <- data.frame(
        dataset = ds, bootstrap = b, modality = bl, complete_rows = nrow(A), features = ncol(A),
        centered_rank = if (nrow(A)) {
          qr(scale(A, scale = FALSE), tol = 1e-10, LAPACK = FALSE)$rank
        } else {
          NA
        }
      )
    }
  }
}

write.csv(
  do.call(rbind, out), file.path(root, "private/representative_ranks.csv"), row.names = FALSE
)
