# Run the observed-data benchmark for the R methods with the original method settings
# OBSERVED_MODE selects "bootstrap" (resampled replicates), "full" (full cohort), or
# "prepare" (build the cached inputs only)

source(file.path(Sys.getenv("CLEAN_ROOT"), "code", "load_dataset.R"))
dataset <- Sys.getenv("BENCH_DATASET")
e <- load_dataset(dataset, observed = TRUE)
mode <- Sys.getenv("OBSERVED_MODE", "bootstrap")
if (mode == "prepare") quit(save = "no")

# Methods to run; only some of them accept incomplete inputs
methods <- strsplit(Sys.getenv("R_METHODS"), ",", fixed = TRUE)[[1]]
na_methods <- intersect(methods, c("rgcca", "gpcca", "block_pls", "gipca", "bsf", "nemo", "bidifac"))
# CIA and O2PLS support two blocks only
if (dataset == "tcga" && length(unique(e$layer_vec)) > 2) {
  methods <- setdiff(methods, c("cia", "o2pls"))
}

# Run settings and output folder
B <- as.integer(Sys.getenv("B_BOOT", "10"))
first <- as.integer(Sys.getenv("BOOT_START", "1"))
cores <- as.integer(Sys.getenv("BOOT_CORES", "1"))
save_models <- Sys.getenv("SAVE_BOOT_MODELS") == "1"
out <- file.path(Sys.getenv("BENCH_WORK_DIR"), if (mode == "full") "full_r" else "bootstrap_r")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Fit each strategy
rows <- list()
for (style in c("NA", "CC_SF", "MI_SF")) {
  selected <- if (style == "NA") na_methods else methods
  if (!length(selected)) next
  if (mode == "bootstrap") {
    # Methods that accept NA inputs and those that do not are run as separate groups
    groups <- if (style == "NA") {
      list(na_methods)
    } else {
      list(intersect(selected, na_methods), setdiff(selected, na_methods))
    }
    result <- dplyr::bind_rows(lapply(groups, function(group) {
      if (!length(group)) {
        return(NULL)
      }
      e$bootstrap_pipe_parallel(
        X_stand = e$X_stand, layer_vec = e$layer_vec, diag_labels = e$diag_labels,
        miss_labels = e$missingness_label, mmse_vec = e$mmse_vec, style = style,
        B = B, boot_start = first, models = group, n_factors = 3L, seed = 100L, mc.cores = cores,
        imputed_list = e$imp_list, save_models = save_models, model_dir = e$CLEAN_MODEL_DIR
      )
    }))
  } else if (mode == "full") {
    # Historical full-fit initializations were not explicitly seeded here
    if (style == "NA") {
      result <- e$run_NA(
        e$X_na, e$lv_na, e$diag_labels, e$missingness_label,
        e$mmse_vec, selected, 3L, save_models, e$CLEAN_MODEL_DIR, style, 0L
      )
    } else if (style == "CC_SF") {
      result <- e$run_complete(
        e$X_cc_sf, e$lv_cc_sf, e$diag_labels, e$missingness_label,
        e$mmse_vec, selected, 3L, save_models, e$CLEAN_MODEL_DIR, style, 0L
      )
    } else {
      raw <- e$run_multiple_complete(
        e$imp_list, e$layer_vec, e$diag_labels, e$diag_labels,
        e$missingness_label, data.frame(MMSE = e$mmse_vec), selected,
        fs = TRUE, n_factors = 3L, save_models = save_models, model_dir = e$CLEAN_MODEL_DIR,
        style = style, boot = 0L
      )
      saveRDS(raw, file.path(out, "MI_SF_raw.rds"))
      result <- e$summarise_mi_results(raw)
    }
  } else {
    stop("Invalid observed mode")
  }
  result$style <- style
  rows[[style]] <- result
}

# Save the combined results
all <- dplyr::bind_rows(rows)
tag <- if (mode == "full") "full" else sprintf("%04d-%04d", first, first + B - 1L)
saveRDS(all, file.path(out, paste0(tag, ".rds")))
write.csv(all, file.path(out, paste0(tag, ".csv")), row.names = FALSE)
