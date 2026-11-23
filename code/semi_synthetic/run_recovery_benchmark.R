# Run one stage of the R semi-synthetic modality-loss experiment (SEMISYNTH_MODE)

source(file.path(Sys.getenv("CLEAN_ROOT"), "code", "load_dataset.R"))
source(file.path(CLEAN_ROOT, "code", "semi_synthetic", "modality_loss_and_recovery.R"))

switch(MODE,
  ground_truth = build_ground_truth(),
  run = run_task(as.integer(Sys.getenv("SEMISYNTH_TASK_ID"))),
  combine = combine_tasks(),
  stop("Unknown mode")
)
