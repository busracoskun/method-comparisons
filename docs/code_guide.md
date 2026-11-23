# Code guide

Start with `run_benchmark.py`. It selects a dataset, analysis, and implementation, sets the paths, and calls the scripts in `code/`. Run `python run_benchmark.py --help` for the command-line options; [settings.md](settings.md) lists the defaults.

## Layout

```
run_benchmark.py            Entry point for all experiments
project_paths.py            Project-relative paths and configuration loading
code/
  load_dataset.R            R session setup; load_dataset() builds one dataset environment
  download_external_sources.py
                            Fetch the pinned GIPCA/BIDIFAC sources (config/external_sources.json)
  preprocessing/            Load and align the prepared ADNI and TCGA inputs
  observed_data/            Observed-data analysis: NA, CC_SF, and MI_SF inputs, bootstrap and
                            full-cohort drivers (R and Python), and re-scoring of saved ADNI fits
  semi_synthetic/           Modality-loss experiment: masks, imputation, recovery metrics, drivers
  methods/                  Method wrappers, evaluation metrics, and bootstrap code
    adni.R, tcga.R          The 18 R methods (dataset-specific versions)
    shared.R                Helpers shared by adni.R and tcga.R
    mofa_cmf.py             MOFA and CMF, and the Python bootstrap
    mofa_cmf_tcga.py        TCGA data loading, SIS, and imputation cache for mofa_cmf.py
  summaries/                Tables from saved or newly run results
  figures/                  Manuscript figures; patient_visualizations/ holds the optional
                            UMAP helpers that need the private coordinate archive
  diagnostics/              Numerical-failure, feature-screening, and runtime checks
```

## Where to start

| Task | Start here |
|---|---|
| Prepare ADNI/TCGA | `code/preprocessing/load_prepared_adni.R`, `load_prepared_tcga.R`, then `code/observed_data/prepare_adni.R`, `prepare_tcga.R` |
| Run observed-data analysis (R) | `code/observed_data/run_observed_benchmark.R` |
| Run observed-data analysis (Python) | `code/observed_data/run_adni_mofa_cmf_bootstrap.py`, `run_adni_mofa_cmf_full_cohort.py`, `run_tcga_mofa_cmf.py` |
| Fit statistical methods | `code/methods/adni.R`, `tcga.R`, `mofa_cmf.py`, `mofa_cmf_tcga.py` |
| Run modality-loss experiments | `code/semi_synthetic/run_recovery_benchmark.R`, `run_mofa_cmf_recovery.py` (both use `modality_loss_and_recovery.R`) |
| Calculate metrics | `score_clusters()` and `score_regression()` in `code/methods/`; `latent_metrics()` in `code/semi_synthetic/modality_loss_and_recovery.R` |
| Combine outputs of new runs | `code/summaries/combine_run_outputs.py` |
| Summarize saved results | `code/summaries/summarize_results.py`, `format_full_cohort_results.py` |
| Make the figures | `code/figures/render_figures.py` |
| Check numerical failures and runtime | `code/diagnostics/` |

The method files keep each model's original wrapper and settings together, so that the calculations used in the paper are unchanged; `run_complete()` selects which wrappers run. The Python MOFA/CMF code calls R (through rpy2) for k-means, SIS screening, missForest, and bootstrap resampling so that seeds and samples match the R pipeline.

For a particular figure or result, use [results_guide.md](results_guide.md). The archived numerical results are in `results/`; new outputs are written to `build/` or the Git-ignored private output directory.
