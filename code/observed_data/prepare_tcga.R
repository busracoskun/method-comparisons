# Build the TCGA observed-data inputs for the NA, CC_SF, and MI_SF strategies
# Sourced after load_prepared_tcga.R. Cached feature sets and imputations are reused only
# when they match the current matrix; otherwise they are recomputed

all_na_cols <- which(colSums(!is.na(X_stand)) == 0L)
if (length(all_na_cols)) {
  X_stand <- X_stand[, -all_na_cols, drop = FALSE]
  X_all <- X_all[, colnames(X_stand), drop = FALSE]
  layer_vec <- layer_vec[colnames(X_stand)]
}

var_dir <- VARIABLE_DIR
dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)

# NA: incomplete matrix passed directly to methods that support missing values
X_na <- X_stand
lv_na <- layer_vec[colnames(X_na)]

# CC_SF: patients observed in every modality, followed by SIS feature screening
ids_cc <- intersect(common_ids, rownames(X_na))
X_cc <- X_na[ids_cc, , drop = FALSE]
X_cc <- X_cc[, colSums(is.na(X_cc)) == 0L, drop = FALSE]
lv_cc <- layer_vec[colnames(X_cc)]

if (nrow(X_cc) < 2L || ncol(X_cc) < 2L) {
  stop("TCGA complete-case matrix is too small for the benchmark.")
}

# Reuse the cached SIS selection when it matches the data, otherwise recompute it
feats_cc_file <- cache_file("tcga_feats_cc_full.rds")
recompute_feats <- TRUE
if (file.exists(feats_cc_file)) {
  candidate_feats <- readRDS(feats_cc_file)
  if (length(candidate_feats) && all(candidate_feats %in% colnames(X_cc))) {
    feats_cc <- candidate_feats
    recompute_feats <- FALSE
    cat("Loading saved TCGA feats_cc from:", feats_cc_file, "\n")
  } else {
    cat("Saved TCGA feats_cc is incompatible with current data; recomputing.\n")
  }
}
if (recompute_feats) {
  cat("Computing TCGA feats_cc and saving to:", feats_cc_file, "\n")
  feats_cc <- sel_SIS(X_cc, diag_labels[rownames(X_cc)])
  saveRDS(feats_cc, file.path(var_dir, basename(feats_cc_file)))
}

X_cc_sf <- X_cc[, feats_cc, drop = FALSE]
lv_cc_sf <- lv_cc[feats_cc]

# MI_SF: M missForest imputations, each seeded with imp_seed + m
M <- 5L
imp_seed <- 10L
imp_list_file <- cache_file(paste0("tcga_imp_list_M", M, "_seed", imp_seed, ".rds"))

# A cached imputation list must match X_na exactly and be complete
valid_imp_list <- function(x) {
  is.list(x) && length(x) == M &&
    all(vapply(x, function(Xi) {
      is.matrix(Xi) && identical(dim(Xi), dim(X_na)) && identical(rownames(Xi), rownames(X_na)) &&
        identical(colnames(Xi), colnames(X_na)) &&
        !anyNA(Xi)
    }, logical(1)))
}

# Reuse valid cached imputations or compute them
if (file.exists(imp_list_file)) {
  cat("Loading saved TCGA imp_list from:", imp_list_file, "\n")
  imp_list <- readRDS(imp_list_file)
  if (!valid_imp_list(imp_list)) {
    cat("Saved TCGA imp_list is incompatible with current data; recomputing.\n")
    imp_list <- NULL
  }
} else {
  imp_list <- NULL
}

if (is.null(imp_list)) {
  cat("Computing TCGA imp_list and saving to:", imp_list_file, "\n")
  imp_list <- lapply(seq_len(M), function(m) {
    set.seed(imp_seed + m)
    missForest::missForest(X_na, verbose = FALSE)$ximp
  })
  saveRDS(imp_list, file.path(var_dir, basename(imp_list_file)))
}

# Save the analysis variables used by all methods
tcga_runtime_vars_file <- file.path(var_dir, "tcga_bootstrap_analysis_variables.rds")
tcga_runtime_vars <- list(
  source_input_dir = TCGA_INPUT_DIR,
  input_files = c(
    mrna = TCGA_MRNA_FILE, methylation = TCGA_METHYLATION_FILE, rppa = TCGA_RPPA_FILE,
    mirna = if (file.exists(TCGA_MIRNA_FILE)) TCGA_MIRNA_FILE else NA_character_,
    clinical = TCGA_CLINICAL_FILE
  ),
  X_stand = X_stand, X_na = X_na, X_cc = X_cc, X_cc_sf = X_cc_sf, X_all = X_all,
  layer_vec = layer_vec, feature_name_map = feature_name_map,
  lv_na = lv_na, lv_cc = lv_cc, lv_cc_sf = lv_cc_sf,
  diag_labels = diag_labels, missingness_label = missingness_label, mmse_vec = mmse_vec,
  regression_name = regression_name, common_ids = common_ids, ids_by_modality = ids_by_modality,
  feats_cc = feats_cc, imp_list_file = imp_list_file, M = M, imp_seed = imp_seed
)
saveRDS(tcga_runtime_vars, tcga_runtime_vars_file)
cat("Saved TCGA bootstrap variables to:", tcga_runtime_vars_file, "\n")
