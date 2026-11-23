# Build the ADNI observed-data inputs for the NA, CC_SF, and MI_SF strategies
# Sourced after load_prepared_adni.R. Feature sets and imputations are cached in
# VARIABLE_DIR (or read from CACHE_INPUT_DIR) so that all methods share them

all_na_cols <- which(colSums(!is.na(X_stand)) == 0)
if (length(all_na_cols)) {
  X_stand <- X_stand[, -all_na_cols, drop = FALSE]
  layer_vec <- layer_vec[colnames(X_stand)]
}

var_dir <- VARIABLE_DIR
dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)

# NA: incomplete matrix passed directly to methods that support missing values
X_na <- X_stand
lv_na <- layer_vec[colnames(X_na)]

# CC_SF: complete cases followed by SIS feature screening
ids_cc <- intersect(intersect(ids_gene, ids_meta), rownames(X_na))
X_cc <- X_na[ids_cc, , drop = FALSE]
X_cc <- X_cc[, colSums(is.na(X_cc)) == 0, drop = FALSE]
lv_cc <- layer_vec[colnames(X_cc)]

# Reuse the cached SIS selection or compute it
feats_cc_file <- cache_file("feats_cc_full.rds")
if (file.exists(feats_cc_file)) {
  cat("Loading saved feats_cc from:", feats_cc_file, "\n")
  feats_cc <- readRDS(feats_cc_file)
} else {
  cat("Computing feats_cc and saving to:", feats_cc_file, "\n")
  feats_cc <- sel_SIS(X_cc, diag_labels[rownames(X_cc)])
  saveRDS(feats_cc, file.path(var_dir, basename(feats_cc_file)))
}

X_cc_sf <- X_cc[, feats_cc, drop = FALSE]
lv_cc_sf <- lv_cc[feats_cc]

# MI_SF: M missForest imputations, each seeded with imp_seed + m
M <- 5L
imp_seed <- 10L

# Reuse the cached imputations or compute them
imp_list_file <- cache_file(paste0("imp_list_M", M, "_seed", imp_seed, ".rds"))
if (file.exists(imp_list_file)) {
  cat("Loading saved imp_list from:", imp_list_file, "\n")
  imp_list <- readRDS(imp_list_file)
} else {
  cat("Computing imp_list and saving to:", imp_list_file, "\n")
  set.seed(imp_seed)
  imp_list <- lapply(seq_len(M), function(m) {
    set.seed(imp_seed + m)
    missForest::missForest(X_na, verbose = FALSE)$ximp
  })
  saveRDS(imp_list, file.path(var_dir, basename(imp_list_file)))
}

# Save the analysis variables used by all methods
analysis_vars_file <- file.path(var_dir, "analysis_variables.rds")
analysis_vars <- list(
  X_stand = X_stand, X_na = X_na, X_cc = X_cc, X_cc_sf = X_cc_sf, layer_vec = layer_vec,
  lv_na = lv_na, lv_cc = lv_cc, lv_cc_sf = lv_cc_sf, diag_labels = diag_labels,
  missingness_label = missingness_label, mmse_vec = mmse_vec, common_ids = common_ids,
  feats_cc = feats_cc, imp_list_file = imp_list_file, M = M, imp_seed = imp_seed
)
saveRDS(analysis_vars, analysis_vars_file)
cat("Saved analysis variables to:", analysis_vars_file, "\n")
