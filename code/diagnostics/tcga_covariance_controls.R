# Covariance controls for the GPCCA failures (usage: <prepared_root>): for the full-cohort CC_SF
# input and bootstrap replicates 1, 7, and 72, compare the rank, eigenvalues, and Matrix::solve
# of each per-view covariance as is, after removing duplicated features, and with a diagonal
# ridge. Bootstrap inputs come from tcga_screening.R's screening replay

# Paths and prepared variables
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gpcca")
v <- readRDS(file.path(root, "tcga/variables/tcga_bootstrap_analysis_variables.rds"))

rows <- list()

# Full-cohort input (boot 0) and selected bootstrap replicates
for (boot in c(0L, 1L, 7L, 72L)) {
  if (boot == 0L) {
    X <- v$X_cc_sf
    unique_n <- nrow(X)
  } else {
    r <- readRDS(file.path(out, "screening_replay", sprintf("boot_%04d.rds", boot)))
    X <- v$X_stand[r$cc_ids, r$selected, drop = FALSE]
    unique_n <- length(unique(r$cc_ids))
  }

  for (mod in unique(v$layer_vec[colnames(X)])) {
    for (policy in c("original", "deduplicated", "initial_ridge")) {
      Z <- X[, v$layer_vec[colnames(X)] == mod, drop = FALSE]
      if (policy == "deduplicated") Z <- Z[, !duplicated(t(Z)), drop = FALSE]

      # Covariance under each policy and whether Matrix::solve inverts it
      C <- cov(Z)
      if (policy == "initial_ridge") C <- C + diag(diag(C), nrow = nrow(C))
      vals <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
      psi <- Matrix::bdiag(list(C))
      residual <- NA_real_
      message <- tryCatch(
        {
          inv <- Matrix::solve(psi, tol = 0)
          residual <- norm(as.matrix(psi %*% inv) - diag(ncol(C)), "F") / sqrt(ncol(C))
          "returned_inverse"
        },
        error = function(e) conditionMessage(e)
      )

      rows[[length(rows) + 1]] <- data.frame(
        boot = boot, policy = policy, modality = mod, patients = nrow(Z),
        unique_patients = unique_n, features = ncol(Z),
        max_centered_rank = min(ncol(Z), unique_n - 1L),
        covariance_rank_relative_1e10 = sum(vals > max(vals) * 1e-10),
        smallest_eigenvalue = min(vals), largest_eigenvalue = max(vals),
        positive_definite_relative_1e10 = min(vals) > max(vals) * 1e-10, Matrix_solve = message,
        inverse_relative_residual = residual
      )
    }
  }
}

write.csv(do.call(rbind, rows), file.path(out, "covariance_controls.csv"), row.names = FALSE)
cat("Wrote", length(rows), "original/deduplicated/ridge per-view covariance controls\n")
