# Preserved source expressions; see docs/reproducibility.md.
# Extracted from final_tcga.R:2000
all_na_cols <- which(colSums(!is.na(X_stand)) == 0L)

# Extracted from final_tcga.R:2001
if (length(all_na_cols)) {
  X_stand <- X_stand[, -all_na_cols, drop = FALSE]
  X_all <- X_all[, colnames(X_stand), drop = FALSE]
  layer_vec <- layer_vec[colnames(X_stand)]
}

# Extracted from final_tcga.R:2007
var_dir <- VARIABLE_DIR

# Extracted from final_tcga.R:2008
dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)

# Extracted from final_tcga.R:2011
X_na <- X_stand

# Extracted from final_tcga.R:2012
lv_na <- layer_vec[colnames(X_na)]

# Extracted from final_tcga.R:2015
ids_cc <- intersect(common_ids, rownames(X_na))

# Extracted from final_tcga.R:2016
X_cc <- X_na[ids_cc, , drop = FALSE]

# Extracted from final_tcga.R:2017
X_cc <- X_cc[, colSums(is.na(X_cc)) == 0L, drop = FALSE]

# Extracted from final_tcga.R:2018
lv_cc <- layer_vec[colnames(X_cc)]

# Extracted from final_tcga.R:2020
if (nrow(X_cc) < 2L || ncol(X_cc) < 2L) {
  stop("TCGA complete-case matrix is too small for the benchmark.")
}

# Extracted from final_tcga.R:2028
feats_cc_file <- cache_file("tcga_feats_cc_full.rds")

# Extracted from final_tcga.R:2029
recompute_feats <- TRUE

# Extracted from final_tcga.R:2031
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

# Extracted from final_tcga.R:2042
if (recompute_feats) {
  cat("Computing TCGA feats_cc and saving to:", feats_cc_file, "\n")
  feats_cc <- sel_SIS(X_cc, diag_labels[rownames(X_cc)])
  saveRDS(feats_cc, file.path(var_dir, basename(feats_cc_file)))
}

# Extracted from final_tcga.R:2048
X_cc_sf <- X_cc[, feats_cc, drop = FALSE]

# Extracted from final_tcga.R:2049
lv_cc_sf <- lv_cc[feats_cc]

# Extracted from final_tcga.R:2055
M <- 5L

# Extracted from final_tcga.R:2056
imp_seed <- 10L

# Extracted from final_tcga.R:2057
imp_list_file <- cache_file(paste0("tcga_imp_list_M", M, "_seed", imp_seed, ".rds"))

# Extracted from final_tcga.R:2062
valid_imp_list <- function(x) {
  is.list(x) &&
    length(x) == M &&
    all(vapply(x, function(Xi) {
      is.matrix(Xi) &&
        identical(dim(Xi), dim(X_na)) &&
        identical(rownames(Xi), rownames(X_na)) &&
        identical(colnames(Xi), colnames(X_na)) &&
        !anyNA(Xi)
    }, logical(1)))
}

# Extracted from final_tcga.R:2074
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

# Extracted from final_tcga.R:2086
if (is.null(imp_list)) {
  cat("Computing TCGA imp_list and saving to:", imp_list_file, "\n")
  imp_list <- lapply(seq_len(M), function(m) {
    set.seed(imp_seed + m)
    missForest::missForest(X_na, verbose = FALSE)$ximp
  })
  saveRDS(imp_list, file.path(var_dir, basename(imp_list_file)))
}

# Extracted from final_tcga.R:2096
tcga_runtime_vars_file <- file.path(var_dir, "tcga_bootstrap_analysis_variables.rds")

# Extracted from final_tcga.R:2098
tcga_runtime_vars <- list(
  source_input_dir = TCGA_INPUT_DIR,
  input_files = c(
    mrna = TCGA_MRNA_FILE,
    methylation = TCGA_METHYLATION_FILE,
    rppa = TCGA_RPPA_FILE,
    mirna = if (file.exists(TCGA_MIRNA_FILE)) TCGA_MIRNA_FILE else NA_character_,
    clinical = TCGA_CLINICAL_FILE
  ),
  X_stand = X_stand,
  X_na = X_na,
  X_cc = X_cc,
  X_cc_sf = X_cc_sf,
  X_all = X_all,
  layer_vec = layer_vec,
  feature_name_map = feature_name_map,
  lv_na = lv_na,
  lv_cc = lv_cc,
  lv_cc_sf = lv_cc_sf,
  diag_labels = diag_labels,
  missingness_label = missingness_label,
  mmse_vec = mmse_vec,
  regression_name = regression_name,
  common_ids = common_ids,
  ids_by_modality = ids_by_modality,
  feats_cc = feats_cc,
  imp_list_file = imp_list_file,
  M = M,
  imp_seed = imp_seed
)

# Extracted from final_tcga.R:2129
saveRDS(tcga_runtime_vars, tcga_runtime_vars_file)

# Extracted from final_tcga.R:2130
cat("Saved TCGA bootstrap variables to:", tcga_runtime_vars_file, "\n")

