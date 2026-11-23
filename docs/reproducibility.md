# Reproducibility

## Inputs and environment

The prepared ADNI and TCGA matrices are not distributed. Their expected paths and required fields are listed in the [README](../README.md). The original cohort preparation is in `code/preprocessing/`; the experiment settings are in `config/default.json` and the analysis drivers.

To configure input paths without editing source files, copy `config/host.example.json` to `config/host.local.json`, then set `adni_input_dir` and `tcga_input_dir`. The local configuration file is ignored by Git. Inputs must come from authorized data sources.

Python dependencies for plots are in `requirements.txt`; experiment dependencies are in `environment.yml`, `requirements_experiments.txt`, and `config/install.R`. Recorded package versions are in `environment/`. Some historical methods require externally sourced R implementations. Their URLs and expected hashes are in `config/external_sources.json`.

## Analyses

The main entry point is `run_benchmark.py`. Use `--check` to validate paths without fitting. Observed runs use CC_SF, MI_SF, and NA strategies; semi-synthetic runs use 600 task IDs across three missingness mechanisms, two rates, and 100 repetitions.

```bash
python run_benchmark.py adni-observed --mode bootstrap --count 1 --check
python run_benchmark.py tcga-observed --mode bootstrap --count 1 --check
python run_benchmark.py adni-semisynthetic --mode ground_truth
python run_benchmark.py adni-semisynthetic --task 1
```

Bootstrap and semi-synthetic results from the paper are already saved under `results/`. Outputs from new runs belong under `private_outputs/`.

Original historical results include a few unavailable records and fitting failures. Reported success counts refer to recorded pipeline returns; they do not certify optimization convergence. Numerical diagnostics are in `results/diagnostics/`.

## Figures

Run `python code/figures/render_figures.py` to regenerate all plots based only on aggregate numerical results. Pass `--figure 5` or `--figure s10` for an individual plot. The figure-generation scripts write to `build/figures/` by default.

### UMAP inputs

Figures 4, 7, S1, S2, and S3 use patient-level UMAP coordinates. These are not included. To regenerate them, use the original **authorized private** coordinate archive:

```bash
python code/figures/render_figures.py --figure 4 --coordinate-archive /path/to/authorized_archive
```

The script `code/figures/embeddings.py` documents the expected files, including:

- `results/historical_manuscript_figures/tables/observed_UMAP_manifest.csv` and matching `UMAP_<DATASET>_<METHOD>_<STRATEGY>.csv` files;
- `results/historical_manuscript_figures/cache/` for the selected tICA reference and refit UMAPs;
- `supplementary/tables/S1_ADNI_semisynthetic_UMAP.csv`, `S2_TCGA_semisynthetic_UMAP.csv`, and `S3_observed_UMAP.csv`.

Keep this archive outside the public repository. The final figures are supplied in `figures/`, subject to redistribution clearance.

### Historical patient-visualization helpers

`select_semisynthetic_examples.py` reproduces the historical tICA repetition-selection
calculation from the included `results/semi_synthetic/adni_semi_synthetic_results.csv` and
`tcga_semi_synthetic_results.csv` files. It writes its generated combined selection under
`build/patient_visualizations/` and does not modify the published selection records.

`export_python_scores.py` requires two private patient-level fitted models; the aggregate
metrics in `results/` are not sufficient. Pass a trusted fit directory containing
`mofa/mofa_NA_0.pkl` and `mofa/mofa_MI_SF_0_imp1.pkl`:

```bash
python code/figures/patient_visualizations/export_python_scores.py \
  --fits-dir /path/to/authorized/fits_python
```

The exported patient IDs and latent scores are written to
`private_outputs/patient_visualizations/python_scores/` by default. The required fitted models
are not included in this repository; only unpickle artifacts from a trusted source.

## Scope

The repository preserves the final reported results and the main scripts used to obtain them. It excludes temporary cluster recovery scripts, obsolete outputs, private intermediate objects, and restricted datasets. Some historic environment and random-initialization states cannot be reconstructed exactly. Reproducing a model run can therefore differ from reproducing a figure using its saved results.
