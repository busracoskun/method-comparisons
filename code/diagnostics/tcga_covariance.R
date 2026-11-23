# Per-modality covariance diagnostics of the TCGA complete-case inputs (X_cc, X_cc_sf): zero
# variance, duplicated columns, near-collinear pairs, and covariance rank. Also checks the
# duplicated methylation pair (RNU5D, RNU5E) in the raw prepared file. Read-only: does not call
# GPCCA, SIS, imputation, or bootstrap code

# Paths and prepared variables
root <- Sys.getenv("BENCH_PREPARED_DIR", "private_outputs")
out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gpcca")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
v <- readRDS(file.path(root, "tcga/variables/tcga_bootstrap_analysis_variables.rds"))

# Covariance summaries for each input and modality
rows <- list()
pairs <- list()
for (input in c("X_cc", "X_cc_sf")) {
  X <- v[[input]]
  lv <- v$layer_vec[colnames(X)]
  for (mod in unique(lv)) {
    Z <- X[, lv == mod, drop = FALSE]

    sdv <- apply(Z, 2, sd)
    C <- cov(Z)
    ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    R <- cor(Z)

    # Pairwise correlations and exact duplicate pairs
    ind <- which(upper.tri(R), arr.ind = TRUE)
    cors <- R[ind]
    dup <- vapply(
      seq_len(nrow(ind)), function(k) identical(unname(Z[, ind[k, 1]]), unname(Z[, ind[k, 2]])),
      logical(1)
    )
    strongest <- which.max(ifelse(dup, -Inf, abs(cors)))

    rows[[length(rows) + 1]] <- data.frame(
      input = input, modality = mod, n = nrow(Z), features = ncol(Z), zero_variance = sum(sdv == 0),
      minimum_sd = min(sdv), duplicate_columns = sum(duplicated(t(Z))),
      exact_duplicate_pairs = sum(dup), abs_r_ge_095 = sum(abs(cors) >= .95),
      abs_r_ge_099 = sum(abs(cors) >= .99), abs_r_ge_0999 = sum(abs(cors) >= .999),
      nonduplicate_abs_r_ge_099 = sum(abs(cors) >= .99 & !dup),
      nonduplicate_abs_r_ge_0999 = sum(abs(cors) >= .999 & !dup),
      max_nonduplicate_abs_r = max(abs(cors[!dup])),
      max_nonduplicate_pair1 = colnames(Z)[ind[strongest, 1]],
      max_nonduplicate_pair2 = colnames(Z)[ind[strongest, 2]],
      covariance_rank_relative_1e10 = sum(ev > max(ev) * 1e-10), smallest_eigenvalue = min(ev),
      largest_eigenvalue = max(ev), stringsAsFactors = FALSE
    )

    # Highly correlated methylation pairs
    if (mod == "methylation") {
      keep <- which(abs(cors) >= .95)
      pairs[[length(pairs) + 1]] <- data.frame(
        input = input, feature1 = colnames(Z)[ind[keep, 1]], feature2 = colnames(Z)[ind[keep, 2]],
        pearson_r = cors[keep], exact_duplicate = dup[keep]
      )
    }
  }
}
write.csv(do.call(rbind, rows), file.path(out, "input_rechecks.csv"), row.names = FALSE)
write.csv(
  do.call(rbind, pairs), file.path(out, "methylation_high_correlation_pairs.csv"), row.names = FALSE
)

# Check the duplicated methylation pair in the raw prepared file
methylation <- file.path(Sys.getenv("TCGA_INPUT_DIR", "data/tcga"), "imeth_fs.csv")
if (!file.exists(methylation)) methylation <- paste0(methylation, ".gz")
raw <- read.csv(methylation, check.names = FALSE)
stopifnot(all(c("RNU5D", "RNU5E") %in% names(raw)))
write.csv(
  data.frame(
    input = "prepared methylation before analysis scaling/SIS", patients = nrow(raw),
    RNU5D_RNU5E_exactly_identical = identical(raw$RNU5D, raw$RNU5E),
    both_selected = all(c("RNU5D", "RNU5E") %in% colnames(v$X_cc_sf))
  ),
  file.path(out, "raw_duplicate_recheck.csv"), row.names = FALSE
)
