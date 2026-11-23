# Per-modality dimensions, zero-variance and duplicated features, covariance rank, and
# invertibility of the prepared complete-case inputs (X_cc, X_cc_sf) for TCGA and ADNI
# (usage: <prepared_root>)

# Paths
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gpcca")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rows <- list()
for (ds in c("tcga", "adni")) {
  v0 <- readRDS(file.path(
    root, ds, "variables",
    if (ds == "tcga") "tcga_bootstrap_analysis_variables.rds" else "analysis_variables.rds"
  ))

  for (name in c("X_cc", "X_cc_sf")) {
    X <- v0[[name]]
    lv <- v0$layer_vec[colnames(X)]
    for (mod in unique(lv)) {
      Z <- X[, lv == mod, drop = FALSE]

      # Covariance spectrum and invertibility
      C <- stats::cov(Z)
      eigenvalues <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
      inv <- tryCatch(
        {
          Matrix::solve(Matrix::Matrix(C, sparse = TRUE), tol = 0)
          "success"
        },
        error = function(e) conditionMessage(e)
      )

      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, input = name, modality = mod, patients = nrow(Z), features = ncol(Z),
        unique_patients = length(unique(rownames(Z))), zero_variance = sum(apply(Z, 2, sd) == 0),
        duplicate_columns = sum(duplicated(t(Z))),
        rank_1e10 = sum(eigenvalues > max(eigenvalues) * 1e-10),
        smallest_eigenvalue = min(eigenvalues), largest_eigenvalue = max(eigenvalues),
        matrix_inverse = inv
      )
    }
  }
}

write.csv(
  do.call(rbind, rows), file.path(out, "full_input_covariance_diagnostics.csv"), row.names = FALSE
)
print(do.call(rbind, rows))
