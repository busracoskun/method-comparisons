# Optional: recreate the MCIA/GFA latent scores for the selected MCAR repetitions behind
# Figures S1 and S2 (usage: <selection.csv> <cache_dir> <results_dir>; launched by
# recreate_supplementary_latents.py). Saved reference fits and scenarios are reused when
# available; the full benchmark is never run

args <- commandArgs(trailingOnly = TRUE)
source(file.path(Sys.getenv("CLEAN_ROOT"), "code/load_dataset.R"))
source(file.path(CLEAN_ROOT, "code/semi_synthetic/modality_loss_and_recovery.R"))

selection <- read.csv(args[1], stringsAsFactors = FALSE)
cache <- args[2]
results <- args[3]

# Write the latent scores of one fit with patient IDs and diagnoses
write_scores <- function(fit, setup, path) {
  if (!identical(fit$status, "success")) stop("Latent fit failed: ", fit$error_message)
  Z <- as.matrix(fit$Z)
  if (is.null(rownames(Z)) || !setequal(rownames(Z), rownames(setup$X_truth))) {
    stop("Latent fit did not preserve patient IDs")
  }
  colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  write.csv(
    data.frame(
      patient_id = rownames(Z), diagnosis = as.character(setup$diag_labels[rownames(Z)]), Z,
      check.names = FALSE
    ),
    path, row.names = FALSE
  )
}

for (dataset in unique(selection$dataset)) {
  Sys.setenv(BENCH_DATASET = dataset)
  setup <- prepare_dataset(dataset)
  folder <- file.path(cache, dataset)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

  # Saved reference fits, if available
  ground <- file.path(results, dataset, "semisynthetic_r/ground_truth/ground_truth_fits.rds")
  truth <- if (file.exists(ground)) readRDS(ground) else NULL
  if (!is.null(truth) && !identical(truth$methods, METHODS)) {
    stop(
      "Saved reference method ordering differs from default; recover original settings before proceeding"
    )
  }
  selected <- selection[selection$dataset == dataset, , drop = FALSE]

  # Reference scores for each selected method
  for (model in unique(selected$model)) {
    path <- file.path(folder, paste0(model, "_reference.csv"))
    if (!file.exists(path)) {
      fit <- if (!is.null(truth)) {
        truth$datasets[[dataset]]$fits[[model]]
      } else {
        fit_latent(
          setup, model, setup$X_truth, 10000L + match(model, METHODS), 0L, "SIM_GROUND_TRUTH"
        )
      }
      write_scores(fit, setup, path)
      saveRDS(list(
        source = if (is.null(truth)) "recreated reference" else ground,
        seed = 10000L + match(model, METHODS), methods = METHODS, n_factors = 3L
      ), paste0(path, ".provenance.rds"))
    }
  }

  # Recreate each selected scenario
  for (i in seq_len(nrow(selected))) {
    r <- selected[i, ]
    task <- map_task(r$task_id)
    if (task$mechanism != "MCAR" || task$target_rate != r$rate) {
      stop("Selection is not the expected MCAR task")
    }
    if (!r$modality %in% setup$modalities) stop("Unknown modality")
    if (!r$model %in% c("mcia", "gfa")) stop("Unexpected selected method")

    path <- file.path(folder, paste0(r$model, "_", r$modality, "_task", r$task_id, ".csv"))
    if (file.exists(path)) next
    # Same scenario-seed formula as run_task() in modality_loss_and_recovery.R
    seed <- as.integer(
      1000000L + r$task_id * 100L + match(dataset, c("adni", "tcga")) * 10L +
        match(r$modality, setup$modalities)
    )

    # Reuse the saved scenario or imputation when available, otherwise recreate it
    simulated_file <- file.path(
      results, dataset, "semisynthetic_r/simulated_data",
      paste0("task_", sprintf("%04d", r$task_id), "_", dataset, "_", r$modality, ".rds")
    )
    imp_path <- file.path(folder, paste0("imputation_", r$modality, "_task", r$task_id, ".rds"))
    if (file.exists(simulated_file)) {
      old <- readRDS(simulated_file)
      if (old$modality != r$modality || old$task$task_id != r$task_id) {
        stop("Saved scenario mismatch")
      }
      imputed <- old$imputation$X
    } else if (file.exists(imp_path)) {
      imputed <- readRDS(imp_path)
    } else {
      sim <- make_missing(
        setup$X_truth, names(setup$layer_vec)[setup$layer_vec == r$modality], "MCAR", r$rate,
        beta = SEMISYNTH_BETA, seed = seed
      )
      set.seed(seed + 500000L)
      imputed <- as.matrix(
        missForest::missForest(sim$X_miss, ntree = 100L, maxiter = 10L, verbose = FALSE)$ximp
      )
      saveRDS(imputed, imp_path)
    }

    # Refit, save the scores, and compare the RV with the recorded value
    fit <- fit_latent(
      setup, r$model, imputed, seed + match(r$model, METHODS), r$task_id, "SIM_MISSFOREST"
    )
    write_scores(fit, setup, path)
    reference <- read.csv(file.path(folder, paste0(r$model, "_reference.csv")), row.names = 1)
    zcols <- grep("^Z", names(reference), value = TRUE)
    recovered_rv <- FactoMineR::coeffRV(
      as.matrix(reference[rownames(fit$Z), zcols]), as.matrix(fit$Z)
    )$rv
    write.csv(data.frame(
      recorded_rv = r$rv, recreated_rv = recovered_rv, scenario_seed = seed,
      imputation_seed = seed + 500000L, fit_seed = seed + match(r$model, METHODS),
      reference_seed = 10000L + match(r$model, METHODS)
    ), paste0(path, ".provenance.csv"), row.names = FALSE)
    # A recreated fit whose RV differs from the recorded one is flagged, never relabelled
    if (abs(recovered_rv - r$rv) > 1e-6) {
      warning(
        "Recreated RV differs from recorded RV: ", path, "; inspect provenance before using UMAPs"
      )
    }
  }
}
