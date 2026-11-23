# GIPCA and BIDIFAC failure diagnostics

These scripts investigate the GIPCA (ADNI) and BIDIFAC (TCGA) bootstrap fitting failures reported in the manuscript. They are path-only adaptations of the recorded diagnostics: fitting settings are unchanged, and no benchmark script is run.

## Inputs

| Variable | Meaning |
|---|---|
| `BENCH_ARCHIVE_DIR` | Original benchmark directory, including `variables/`, the bootstrap result CSVs, and `vendor/bidifac/BIDIFAC.R` |
| `DIAGNOSTIC_OUTPUT_DIR` | Output folder; default `private_outputs/diagnostics/gipca_bidifac` (Git-ignored) |
| `GIPCA_SOURCE_DIR`, `BIDIFAC_FILE` | Optional locations of the external method sources |

## Scripts

| Script | Purpose |
|---|---|
| `bootstrap_input_checks.R` | Per-replicate input checks (distinct patients, missing rows/columns, per-modality support) |
| `distinct_observation_support.R` | Distinct observed patients per feature, overall and per replicate |
| `fit_status_audit.py` | Validates the saved fit records and summarizes failures; run after the two scripts above, whose private outputs it reads |
| `gipca_failure_reproduction.R` | Reproduces GIPCA failures on selected ADNI replicates |
| `gipca_aliased_regression_trace.R` | Traces aliased regressions inside GIPCA |
| `bidifac_failure_probe.R` | Locates the failing SVD call inside BIDIFAC |
| `representative_matrix_ranks.R` | Ranks of centered complete-case blocks in representative replicates |

Outputs under `private/` may contain patient-level records and stay local. Outputs under `public/` contain aggregate counts only; the released copies are `results/diagnostics/gipca_bidifac_*.csv`.
