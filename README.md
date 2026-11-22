# Multi-omics integration with missing modalities

Code and saved results for **Benchmarking statistical dimension reduction frameworks for integrative multi-omics analysis with missing modalities**.

The study compares 20 methods using ADNI and TCGA. Observed-data analyses compare complete cases, imputation, and direct fitting of incomplete inputs. Semi-synthetic analyses measure recovery after removing a modality under MCAR, MAR, and MNAR.

## Find a file

| Folder | Contents |
|---|---|
| [`code/`](code/) | Data preparation, model fitting, metrics, and plotting scripts |
| [`results/`](results/) | Saved bootstrap records, recovery metrics, and summary tables |
| [`figures/`](figures/) | Main Figures 1–8 and Supplementary Figures S1–S10 |
| [`docs/results_guide.md`](docs/results_guide.md) | Figure-to-script and result-to-source index |
| [`docs/code_guide.md`](docs/code_guide.md) | Short guide to the analysis code |

## Regenerate figures and tables

With Python 3.11 or newer:

```bash
python -m pip install -r requirements.txt
python code/figures/render_figures.py
python code/summarize_results.py --verify
python code/build_supplementary_workbook.py
```

Results are written to `build/`, not over the archived results. To regenerate one figure, use `--figure 2` or `--figure s10`.

Figures 4, 7, S1, S2, and S3 need a **private, authorized coordinate archive**. Their published images are included, but patient-level coordinates are not. See [UMAP inputs](docs/reproducibility.md#umap-inputs).

## Rerun the experiments

The original ADNI and TCGA inputs are **not included**. Put authorized prepared inputs in `data/` or configure private paths in `config/host.local.json`.

| Expected file | Data |
|---|---|
| `data/adni/rnaseq.csv` (or `.csv.gz`) | Gene expression; `RID` and numeric features |
| `data/adni/meta.csv` | Metabolomics; `RID` and numeric features |
| `data/adni/clin.csv` | Diagnosis and MMSE by `RID` |
| `data/tcga/irsem_fs.csv` | mRNA expression |
| `data/tcga/imeth_fs.csv` | DNA methylation |
| `data/tcga/irppa_fs.csv` | RPPA |
| `data/tcga/clin_common.csv` | Disease group, Karnofsky score, and missingness category |

Each assay matrix has patients in rows and features in columns. Missing values are coded as `NA`. Use the prepared inputs described in [reproducibility](docs/reproducibility.md), not arbitrary reprocessed datasets.

Install experiment dependencies with `environment.yml` and `config/install.R`. Check inputs before running an experiment:

```bash
python run_benchmark.py adni-observed --mode bootstrap --count 1 --check
python run_benchmark.py adni-observed --mode bootstrap --count 1
python run_benchmark.py tcga-semisynthetic --mode ground_truth
python run_benchmark.py tcga-semisynthetic --task 1
```

See [reproducibility](docs/reproducibility.md) for setup, expected outputs, and limitations.

## Data and reuse

Patient-level datasets, imputations, model objects, latent scores, and bootstrap membership are excluded. The `results/` folder contains shareable run-level metrics and summary data.
