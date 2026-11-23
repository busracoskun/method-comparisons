#!/usr/bin/env python3
"""Run one ADNI observed-data bootstrap replicate (BOOT_ID) for MOFA and CMF.

The inputs and imputations are the ones prepared by the R pipeline, and the replicate is drawn
with R's sampler, so MOFA/CMF see exactly the same bootstrap samples as the R methods.
"""

from __future__ import annotations

import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro

# SETTINGS

BENCH_DIR = Path(os.environ["BENCH_WORK_DIR"])
PYTHON_METHODS_FILE = Path(__file__).resolve().parents[2] / "code/methods" / "mofa_cmf.py"
ANALYSIS_VARS_FILE = BENCH_DIR / "variables" / "analysis_variables.rds"
IMP_LIST_FILE = BENCH_DIR / "variables" / "imp_list_M5_seed10.rds"
OUT_DIR = BENCH_DIR / "results_boot_python"
external_imp = (
    Path(os.environ.get("CACHE_INPUT_DIR") or BENCH_DIR / "variables") / "imp_list_M5_seed10.rds"
)
if external_imp.exists():
    IMP_LIST_FILE = external_imp

BOOT_ID = int(os.environ["BOOT_ID"])
MODELS = tuple(os.environ.get("PY_METHODS", "mofa,cmf").split(","))
STYLES = ("NA", "CC_SF", "MI_SF")
N_FACTORS = 3
BOOTSTRAP_SEED = 100
M_IMPUTATIONS = 5
SF_FREQ_THRESHOLD = 4

# Load the shared MOFA/CMF module; registering it in sys.modules is required for its dataclass
if not PYTHON_METHODS_FILE.exists():
    raise FileNotFoundError(f"Missing mofa_cmf.py: {PYTHON_METHODS_FILE}")
spec = importlib.util.spec_from_file_location("python_methods_existing", PYTHON_METHODS_FILE)
pm = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = pm
spec.loader.exec_module(pm)


# RDS CONVERSION HELPERS


def r_matrix_to_df(x) -> pd.DataFrame:
    arr = np.asarray(x, dtype=float)
    rn = [str(v) for v in ro.r["rownames"](x)]
    cn = [str(v) for v in ro.r["colnames"](x)]
    return pd.DataFrame(arr, index=pd.Index(rn), columns=pd.Index(cn))


def named_character_to_series(x) -> pd.Series:
    values_r = ro.r["as.character"](x)
    is_na = np.asarray(ro.r["is.na"](values_r), dtype=bool)
    values = [None if is_na[i] else str(values_r[i]) for i in range(len(values_r))]
    names = [str(v) for v in x.names]
    return pd.Series(values, index=pd.Index(names), dtype="object")


def named_numeric_to_series(x) -> pd.Series:
    values = np.asarray(x, dtype=float)
    names = [str(v) for v in x.names]
    return pd.Series(values, index=pd.Index(names), dtype=float)


# LOAD THE VARIABLES USED BY THE R ADNI BOOTSTRAPS

if not ANALYSIS_VARS_FILE.exists():
    raise FileNotFoundError(f"Missing analysis variables: {ANALYSIS_VARS_FILE}")
if not IMP_LIST_FILE.exists():
    raise FileNotFoundError(f"Missing imputation cache: {IMP_LIST_FILE}")

analysis = ro.r["readRDS"](str(ANALYSIS_VARS_FILE))
X_stand = r_matrix_to_df(analysis.rx2("X_stand"))
layer_r = analysis.rx2("layer_vec")
layer_vec = pd.Series(
    [str(v) for v in layer_r], index=pd.Index([str(v) for v in layer_r.names]), dtype="object"
)
diag_labels = named_character_to_series(analysis.rx2("diag_labels"))
miss_labels = named_character_to_series(analysis.rx2("missingness_label"))
mmse_vec = named_numeric_to_series(analysis.rx2("mmse_vec"))
common_ids = pd.Index([str(v) for v in analysis.rx2("common_ids")])

imp_r = ro.r["readRDS"](str(IMP_LIST_FILE))
imp_list = [r_matrix_to_df(imp_r[i]) for i in range(len(imp_r))]
if len(imp_list) != M_IMPUTATIONS:
    raise RuntimeError(f"Expected {M_IMPUTATIONS} imputations, found {len(imp_list)}.")

settings = pm.Settings(
    data_dir=BENCH_DIR.parent / "preprocessed_data", out_dir=OUT_DIR,
    model_dir=BENCH_DIR / "Bootstraps_python", save_models=False, n_factors=N_FACTORS,
    m_imputations=M_IMPUTATIONS, imputation_seed=10, bootstrap_B=1, bootstrap_seed=BOOTSTRAP_SEED,
    n_jobs_bootstrap=1, styles_to_run=STYLES, models=MODELS,
    sf_frequency_threshold=SF_FREQ_THRESHOLD, use_r_bootstrap_sampling=True,
    cmf_lambda=0.1, cmf_niter=30, cmf_method="als", cmf_use_item_view_info=True
)
OUT_DIR.mkdir(parents=True, exist_ok=True)


# FITTING AND SCORING

# k-means seeds follow the R pipeline: 100 + boot for diagnosis, 200 + boot for missingness
METRICS = [
    "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag", "R2",
    "rand_miss", "ari_miss", "chi_miss", "db_miss", "pve_miss",
]


def fit_and_score(model_name, X, layer, y_diag, y_miss, y_mmse, style, boot, fit_seed, imp=None):
    if model_name == "mofa":
        Z, W, _ = pm.fit_mofa(X, layer.reindex(X.columns), N_FACTORS, fit_seed)
    elif model_name == "cmf":
        Z, W, _ = pm.fit_cmf(
            X, layer.reindex(X.columns), N_FACTORS, fit_seed,
            lambda_=0.1, niter=30, method="als", use_item_view_info=True
        )
    else:
        raise ValueError(f"Unsupported model: {model_name}")

    diag_score = pm.score_clusters(
        Z, y_diag.reindex(X.index), X_original=X, W_loadings=W, seed=100 + boot
    )
    miss_score = pm.score_clusters(
        Z, y_miss.reindex(X.index), X_original=X, W_loadings=W, seed=200 + boot
    )
    return {
        "model": model_name, "model_file": np.nan,
        "rand_diag": diag_score["rand"], "ari_diag": diag_score["ari"],
        "chi_diag": diag_score["chi"], "db_diag": diag_score["db"], "pve_diag": diag_score["pve"],
        "R2": pm.score_regression(Z, y_mmse.reindex(X.index)),
        "rand_miss": miss_score["rand"], "ari_miss": miss_score["ari"],
        "chi_miss": miss_score["chi"], "db_miss": miss_score["db"], "pve_miss": miss_score["pve"],
        "style": style, "boot": boot, "imp": imp,
    }


def failure_result(model, style, boot, error_message, attempted=1):
    """Record a failed fit with NA metrics."""
    row = {"model": model, "model_file": np.nan}
    for metric in METRICS:
        row[metric] = np.nan
    row.update({
        "style": style, "boot": boot, "status": "failed", "error_message": str(error_message),
        "n_attempted": attempted, "n_successful": 0, "n_failed": attempted,
        "elapsed_seconds": np.nan,
    })
    return row


# BOOTSTRAP SAMPLE (SEED 100 + BOOT_ID, R SAMPLER)

ids_all = X_stand.index.astype(str).to_numpy()
sampling_seed = BOOTSTRAP_SEED + BOOT_ID
sampled_ids = pm.sample_ids_for_bootstrap(ids_all, seed=sampling_seed, use_r=True)
new_index = pd.Index([f"{rid}_b{BOOT_ID}_{i + 1}" for i, rid in enumerate(sampled_ids)])
X_b = X_stand.reindex(sampled_ids).copy()
X_b.index = new_index


def bootstrap_series(series: pd.Series) -> pd.Series:
    return pd.Series(series.reindex(sampled_ids).to_numpy(), index=new_index, name=series.name)


diag_b = bootstrap_series(diag_labels)
miss_b = bootstrap_series(miss_labels)
mmse_b = bootstrap_series(mmse_vec)
is_cc = np.isin(sampled_ids, common_ids.to_numpy())


# NA: FIT THE INCOMPLETE MATRIX

rows = []
for model in MODELS:
    started = time.perf_counter()
    try:
        result = fit_and_score(
            model_name=model, X=X_b, layer=layer_vec, y_diag=diag_b, y_miss=miss_b, y_mmse=mmse_b,
            style="NA", boot=BOOT_ID, fit_seed=sampling_seed
        )
        result.update({
            "status": "success", "error_message": np.nan, "n_attempted": 1, "n_successful": 1,
            "n_failed": 0, "elapsed_seconds": time.perf_counter() - started,
        })
        rows.append(result)
    except Exception as exc:
        row = failure_result(model, "NA", BOOT_ID, exc, attempted=1)
        row["elapsed_seconds"] = time.perf_counter() - started
        rows.append(row)


# CC_SF: COMPLETE CASES, SIS SCREENING WITHIN THE REPLICATE

X_cc = X_b.loc[is_cc].copy()
X_cc = X_cc.dropna(axis=1, how="any")
try:
    selected_cc = pm.sel_SIS(X_cc, diag_b.reindex(X_cc.index), seed=sampling_seed)
    X_cc_sf = X_cc.loc[:, selected_cc]
    for model in MODELS:
        started = time.perf_counter()
        try:
            result = fit_and_score(
                model_name=model, X=X_cc_sf, layer=layer_vec, y_diag=diag_b, y_miss=miss_b,
                y_mmse=mmse_b, style="CC_SF", boot=BOOT_ID, fit_seed=sampling_seed
            )
            result.update({
                "status": "success", "error_message": np.nan, "n_attempted": 1, "n_successful": 1,
                "n_failed": 0, "elapsed_seconds": time.perf_counter() - started,
            })
            rows.append(result)
        except Exception as exc:
            row = failure_result(model, "CC_SF", BOOT_ID, exc, attempted=1)
            row["elapsed_seconds"] = time.perf_counter() - started
            rows.append(row)
except Exception as exc:
    for model in MODELS:
        rows.append(failure_result(model, "CC_SF", BOOT_ID, f"SIS failed: {exc}", attempted=1))


# MI_SF: SIS PER IMPUTATION WITH A SELECTION-FREQUENCY THRESHOLD

imp_list_b = []
for Xi in imp_list:
    Xi_b = Xi.reindex(sampled_ids).copy()
    Xi_b.index = new_index
    imp_list_b.append(Xi_b)

try:
    feature_lists = [pm.sel_SIS(Xi, diag_b, seed=sampling_seed) for Xi in imp_list_b]
    frequency = {}
    for feature_list in feature_lists:
        for feature in feature_list:
            frequency[feature] = frequency.get(feature, 0) + 1
    common_features = [
        feature for feature, count in frequency.items() if count >= SF_FREQ_THRESHOLD
    ]
    if not common_features:
        raise RuntimeError("MI_SF selected zero common features.")

    for model in MODELS:
        started = time.perf_counter()
        imp_results = []
        imp_errors = []
        # Imputation i is fit with seed sampling_seed + i; metrics average the successful fits
        for imp_idx, Xi in enumerate(imp_list_b, start=1):
            Xi_sf = Xi.loc[:, common_features]
            try:
                result = fit_and_score(
                    model_name=model, X=Xi_sf, layer=layer_vec, y_diag=diag_b, y_miss=miss_b,
                    y_mmse=mmse_b, style="MI_SF", boot=BOOT_ID,
                    fit_seed=(sampling_seed + imp_idx), imp=imp_idx
                )
                imp_results.append(result)
            except Exception as exc:
                imp_errors.append(f"imp{imp_idx}: {exc}")

        n_success = len(imp_results)
        n_failed = M_IMPUTATIONS - n_success
        if n_success:
            raw = pd.DataFrame(imp_results)
            row = {"model": model, "model_file": np.nan}
            for metric in METRICS:
                row[metric] = pd.to_numeric(raw[metric], errors="coerce").mean()
            row.update({
                "style": "MI_SF", "boot": BOOT_ID,
                "status": "success" if n_failed == 0 else "partial",
                "error_message": np.nan if not imp_errors else " | ".join(imp_errors),
                "n_attempted": M_IMPUTATIONS, "n_successful": n_success, "n_failed": n_failed,
                "elapsed_seconds": time.perf_counter() - started,
            })
            rows.append(row)
        else:
            row = failure_result(
                model, "MI_SF", BOOT_ID, " | ".join(imp_errors), attempted=M_IMPUTATIONS
            )
            row["elapsed_seconds"] = time.perf_counter() - started
            rows.append(row)
except Exception as exc:
    for model in MODELS:
        rows.append(failure_result(
            model, "MI_SF", BOOT_ID, f"MI_SF SIS failed: {exc}", attempted=M_IMPUTATIONS
        ))


# SAVE

result_df = pd.DataFrame(rows)
column_order = [
    "model", "model_file", "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag", "R2",
    "rand_miss", "ari_miss", "chi_miss", "db_miss", "pve_miss", "style", "boot", "status",
    "error_message", "n_attempted", "n_successful", "n_failed", "elapsed_seconds",
]
for column in column_order:
    if column not in result_df:
        result_df[column] = np.nan
result_df = result_df[column_order]
result_df = result_df.sort_values(["boot", "style", "model"]).reset_index(drop=True)
out_file = OUT_DIR / (f"boot_ALL_styles_mofa_cmf_{BOOT_ID:04d}.csv")
result_df.to_csv(out_file, index=False)

print("============================================================")
print(f"BOOTSTRAP {BOOT_ID} COMPLETE")
print(f"Sampling seed: {sampling_seed}")
print(f"Rows: {len(result_df)}")
print(
    result_df[["model", "style", "status", "n_attempted", "n_successful", "n_failed"]].to_string(
        index=False
    )
)
print(f"Saved: {out_file}")
print("============================================================")
