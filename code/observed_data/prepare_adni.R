# Preserved source expressions; see docs/reproducibility.md.
# Extracted from final.R:1739
all_na_cols <- which(colSums(!is.na(X_stand)) == 0)

# Extracted from final.R:1740
if (length(all_na_cols)) {
  X_stand <- X_stand[, -all_na_cols, drop = FALSE]
  layer_vec <- layer_vec[colnames(X_stand)]
}

# Extracted from final.R:1745
var_dir <- VARIABLE_DIR

# Extracted from final.R:1746
dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)

# Extracted from final.R:1749
X_na <- X_stand

# Extracted from final.R:1750
lv_na <- layer_vec[colnames(X_na)]

# Extracted from final.R:1753
ids_cc <- intersect(intersect(ids_gene, ids_meta), rownames(X_na))

# Extracted from final.R:1754
X_cc <- X_na[ids_cc, , drop = FALSE]

# Extracted from final.R:1755
X_cc <- X_cc[, colSums(is.na(X_cc)) == 0, drop = FALSE]

# Extracted from final.R:1756
lv_cc <- layer_vec[colnames(X_cc)]

# Extracted from final.R:1760
feats_cc_file <- cache_file("feats_cc_full.rds")

# Extracted from final.R:1762
if (file.exists(feats_cc_file)) {
  cat("Loading saved feats_cc from:", feats_cc_file, "\n")
  feats_cc <- readRDS(feats_cc_file)
} else {
  cat("Computing feats_cc and saving to:", feats_cc_file, "\n")
  feats_cc <- sel_SIS(X_cc, diag_labels[rownames(X_cc)])
  saveRDS(feats_cc, file.path(var_dir, basename(feats_cc_file)))
}

# Extracted from final.R:1771
X_cc_sf <- X_cc[, feats_cc, drop = FALSE]

# Extracted from final.R:1772
lv_cc_sf <- lv_cc[feats_cc]

# Extracted from final.R:1776
M <- 5L

# Extracted from final.R:1777
imp_seed <- 10L

# Extracted from final.R:1778
imp_list_file <- cache_file(paste0("imp_list_M", M, "_seed", imp_seed, ".rds"))

# Extracted from final.R:1780
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

# Extracted from final.R:1794
analysis_vars_file <- file.path(var_dir, "analysis_variables.rds")

# Extracted from final.R:1796
analysis_vars <- list(
  X_stand = X_stand,
  X_na = X_na,
  X_cc = X_cc,
  X_cc_sf = X_cc_sf,
  layer_vec = layer_vec,
  lv_na = lv_na,
  lv_cc = lv_cc,
  lv_cc_sf = lv_cc_sf,
  diag_labels = diag_labels,
  missingness_label = missingness_label,
  mmse_vec = mmse_vec,
  common_ids = common_ids,
  feats_cc = feats_cc,
  imp_list_file = imp_list_file,
  M = M,
  imp_seed = imp_seed
)

# Extracted from final.R:1815
saveRDS(analysis_vars, analysis_vars_file)

# Extracted from final.R:1816
cat("Saved analysis variables to:", analysis_vars_file, "\n")

