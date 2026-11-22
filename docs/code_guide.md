# Code guide

Start with `run_benchmark.py`. It selects a dataset, analysis, and implementation. The scientific calculations are in `code/`.

| Task | Start here |
|---|---|
| Prepare ADNI/TCGA | `code/preprocessing/load_prepared_adni.R`, `load_prepared_tcga.R` |
| Run observed-data analysis (R) | `code/observed_data/run_observed_benchmark.R` |
| Run observed-data analysis (Python) | `code/observed_data/run_adni_mofa_cmf_bootstrap.py`, `run_tcga_mofa_cmf.py` |
| Fit statistical methods | `code/methods/adni.R`, `tcga.R`, `python_methods.py`, `python_methods_tcga.py` |
| Run modality-loss experiments | `code/semi_synthetic/run_recovery_benchmark.R`, `run_mofa_cmf_recovery.py` |
| Calculate metrics | `code/metrics/` and the original method wrappers |
| Summarize saved results | `code/summarize_results.py` |
| Make the figures | `code/figures/render_figures.py` |
| Check numerical failures and runtime | `code/diagnostics/` |

The long method files contain the original model-specific functions and settings. They are kept together to preserve the calculations used in the paper. The wrappers select which functions run.

For a particular figure or result, use [results_guide.md](results_guide.md). The archived numerical results are in `results/`; new outputs are written to `build/` or the Git-ignored private output directory.
