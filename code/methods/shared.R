# Model-running helpers shared by the ADNI and TCGA method files

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Historical GIPCA source was downloaded from mutable GitHub master
# Supply the exact source snapshot to reproduce this method; never substitute it
load_gipca_once <- local({
  loaded <- FALSE
  function() {
    if (loaded) {
      return(invisible(TRUE))
    }
    folder <- Sys.getenv(
      "GIPCA_SOURCE_DIR", file.path(CLEAN_ROOT, "external_sources", "gipca", "Functions")
    )

    files <- list.files(folder, pattern = "\\.R$", full.names = TRUE)
    if (!length(files)) {
      stop(
        "GIPCA source snapshot missing. See docs/reproducibility.md; supply the documented source snapshot."
      )
    }
    for (f in files) source(f, local = .GlobalEnv)
    loaded <<- TRUE
    invisible(TRUE)
  }
})

# Save one fitted model under <base_dir>/<model>/
save_boot_model <- function(obj, model_name, style, boot, imp = NULL, base_dir = CLEAN_MODEL_DIR) {
  model_name <- tolower(model_name)
  style <- toupper(style)
  model_dir <- file.path(base_dir, model_name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  imp_part <- if (!is.null(imp)) paste0("_imp", imp) else ""
  model_file <- file.path(model_dir, paste0(model_name, "_", style, "_", boot, imp_part, ".rds"))
  saveRDS(obj, model_file)
  model_file
}

# Fit methods that accept incomplete inputs (the NA strategy) on the incomplete matrix
run_NA <- function(
  X_na, layer_vec, diag_labels, miss_labels, mmse_vec, models, n_factors = 3,
  save_models = FALSE, model_dir = CLEAN_MODEL_DIR, style = "NA", boot = NA_integer_
) {
  run_complete(
    X = X_na, layer_vec = layer_vec, diag_labels = diag_labels, miss_labels = miss_labels,
    mmse_vec = mmse_vec, models = models, n_factors = n_factors, save_models = save_models,
    model_dir = model_dir, style = style, boot = boot
  )
}

# Fit methods on each imputed matrix. With fs = TRUE, SIS is run on every imputation and
# only features selected in at least 4 imputations are kept (frequency threshold)
run_multiple_complete <- function(
  imp_list, layer_vec, fs_diag_labels, cluster_diag_labels, cluster_miss_labels, clin_df, models,
  fs = FALSE, n_factors = 3, save_models = FALSE, model_dir = CLEAN_MODEL_DIR,
  style = "MI", boot = NA_integer_
) {
  if (fs) {
    feat_lists <- lapply(
      imp_list, function(Xi) sel_SIS(X_complete = Xi, label_vec = fs_diag_labels[rownames(Xi)])
    )
    freq <- table(unlist(feat_lists))
    common_feats <- names(freq)[freq >= 4]
    if (!length(common_feats)) stop("frequency threshold not met")
  }

  # Fit the models on each imputed matrix
  dplyr::bind_rows(lapply(seq_along(imp_list), function(i) {
    Xi <- imp_list[[i]]
    if (fs) {
      Xi <- Xi[, common_feats, drop = FALSE]
      layer_sub <- layer_vec[common_feats]
    } else {
      layer_sub <- layer_vec[colnames(Xi)]
    }

    diag_i <- factor(cluster_diag_labels[rownames(Xi)])
    miss_i <- factor(cluster_miss_labels[rownames(Xi)])
    mmse_i <- if (!is.null(clin_df)) setNames(clin_df[rownames(Xi), "MMSE"], rownames(Xi)) else NULL
    run_complete(
      X = Xi, layer_vec = layer_sub, diag_labels = diag_i, miss_labels = miss_i, mmse_vec = mmse_i,
      models = models, n_factors = n_factors, save_models = save_models, model_dir = model_dir,
      style = style, boot = boot, imp = i
    ) %>% dplyr::mutate(imp = i)
  }))
}
