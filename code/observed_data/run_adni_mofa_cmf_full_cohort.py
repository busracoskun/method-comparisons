#!/usr/bin/env python3
"""Fit MOFA and CMF once on the full ADNI cohort for the NA, CC_SF, and MI_SF strategies.

Uses the mofa_cmf.py wrappers for data loading, SIS, missForest, model fitting, and scoring.
Each fit is saved as RDS (Z, W, metadata; readable by evaluate_saved_adni_fits.R) and as PKL,
and the MOFA/CMF rows of the full-run result tables are replaced.
"""

from __future__ import annotations

import importlib.util
import os
import sys
import pickle
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np
import pandas as pd

# SETTINGS

BENCH_DIR = Path(os.environ["BENCH_WORK_DIR"])
DATA_DIR = Path(os.environ["ADNI_INPUT_DIR"])
RESULTS_FULL_DIR = BENCH_DIR / "full_python"
FULL_MODEL_DIR = BENCH_DIR / "fits_python"
PYTHON_METHODS_FILE = Path(__file__).resolve().parents[2] / "code/methods" / "mofa_cmf.py"

MODELS = tuple(os.environ.get("PY_METHODS", "mofa,cmf").split(","))
N_FACTORS = int(os.environ.get("N_FACTORS", "3"))
M_IMPUTATIONS = int(os.environ.get("M_IMPUTATIONS", "5"))
IMPUTATION_SEED = int(os.environ.get("IMPUTATION_SEED", "10"))
SF_FREQ_THRESHOLD = int(os.environ.get("SF_FREQ_THRESHOLD", "4"))

# Fixed model-fitting seeds for every (strategy, model, imputation)
SEEDS = {
    ("NA", "mofa", None): 50001, ("NA", "cmf", None): 50002, ("CC_SF", "mofa", None): 50011,
    ("CC_SF", "cmf", None): 50012
}
for imp in range(1, M_IMPUTATIONS + 1):
    SEEDS[("MI_SF", "mofa", imp)] = 50100 + imp
    SEEDS[("MI_SF", "cmf", imp)] = 50200 + imp

# Load the shared MOFA/CMF module; @dataclass needs the module registered in sys.modules
if not PYTHON_METHODS_FILE.exists():
    raise FileNotFoundError(f"Cannot find mofa_cmf.py at: {PYTHON_METHODS_FILE}")
spec = importlib.util.spec_from_file_location("python_methods_existing", PYTHON_METHODS_FILE)
pm = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = pm
spec.loader.exec_module(pm)

RESULT_COLS = [
    "model", "model_file", "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag", "R2",
    "rand_miss", "ari_miss", "chi_miss", "db_miss", "pve_miss", "style", "run",
]
RAW_RESULT_COLS = RESULT_COLS[:13] + ["imp"] + RESULT_COLS[13:]


# RESULT-TABLE HELPERS


def backup_existing_results() -> None:
    RESULTS_FULL_DIR.mkdir(parents=True, exist_ok=True)

    # Copy any existing full-run tables into a timestamped backup folder
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = RESULTS_FULL_DIR / f"backup_before_mofa_cmf_{timestamp}"
    backup_dir.mkdir(parents=True, exist_ok=True)
    targets = [
        "full_NA_all.csv", "full_NA_all.rds", "full_CC_SF_all.csv", "full_CC_SF_all.rds",
        "full_MI_SF_raw_all.csv", "full_MI_SF_raw_all.rds", "full_MI_SF_all.csv",
        "full_MI_SF_all.rds", "full_ALL_styles.csv", "full_ALL_styles.rds",
    ]
    copied = 0
    for name in targets:
        src = RESULTS_FULL_DIR / name
        if src.exists():
            shutil.copy2(src, backup_dir / name)
            copied += 1
    print(f"Backed up {copied} existing result files to: {backup_dir}")


# Existing table, or an empty one with the expected columns
def read_csv_or_empty(path: Path, cols: Sequence[str]) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame(columns=list(cols))


def save_csv_and_rds(df: pd.DataFrame, csv_path: Path, rds_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    import rpy2.robjects as ro

    # Write the same table as RDS through R
    ro.globalenv["csv_path"] = str(csv_path)
    ro.globalenv["rds_path"] = str(rds_path)
    ro.r("""
      df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
      saveRDS(df, rds_path)
    """)


def replace_models(
    existing: pd.DataFrame, new_rows: pd.DataFrame, cols: Sequence[str]
) -> pd.DataFrame:
    """Replace the MOFA/CMF rows of an existing result table, keeping all other methods."""

    # Drop old MOFA/CMF rows and append the new ones
    if existing.empty:
        out = new_rows.copy()
    else:
        existing = existing.loc[~existing["model"].astype(str).str.lower().isin(MODELS)].copy()
        out = pd.concat([existing, new_rows], ignore_index=True)

    # Expected columns first, any extra columns after
    for col in cols:
        if col not in out.columns:
            out[col] = np.nan
    extra_cols = [c for c in out.columns if c not in cols]
    return out[list(cols) + extra_cols]


# SCORING AND SAVING FITS


def score_clusters_compatible(
    Z: pd.DataFrame, labels: pd.Series, X_original: Optional[pd.DataFrame],
    W_loadings: Optional[pd.DataFrame], seed: int
) -> Dict[str, float]:
    """score_clusters(), except single-class labels get the values the saved R fits record.

    With one class (e.g. CC_SF missingness labels, all "CC"), Rand = ARI = 1, CH = NA, DB = 0,
    and only PVE is computed.
    """

    # More than one class: use the regular scoring
    labels = labels.reindex(Z.index)
    keep = labels.notna()
    n_classes = labels.loc[keep].astype(str).nunique()
    if n_classes != 1:
        return pm.score_clusters(Z, labels, X_original=X_original, W_loadings=W_loadings, seed=seed)

    # Single class: clustering metrics are fixed, only PVE is computed
    pve = np.nan
    if X_original is not None and W_loadings is not None:
        feats = X_original.columns.intersection(W_loadings.index, sort=False)
        if len(feats) > 0:
            X_use = X_original.loc[Z.index, feats].to_numpy(dtype=float)
            W_use = W_loadings.loc[feats].to_numpy(dtype=float)
            X_hat = Z.to_numpy(dtype=float) @ W_use.T
            var_x = np.nansum(np.nanvar(X_use, axis=0, ddof=1))
            var_hat = np.nansum(np.var(X_hat, axis=0, ddof=1))
            if np.isfinite(var_x) and var_x > 0:
                pve = float(var_hat / var_x)
    return {"rand": 1.0, "ari": 1.0, "chi": np.nan, "db": 0.0, "pve": pve}


def save_model_rds_and_pkl(
    model_name: str, style: str, seed: int, Z: pd.DataFrame, W: Optional[pd.DataFrame],
    fit_obj: Any, imp: Optional[int] = None
) -> str:
    """Save Z, W, and metadata as RDS (same layout as the R fits) and the raw fit as PKL."""

    # Output file names for this fit
    model_dir = FULL_MODEL_DIR / model_name.lower()
    model_dir.mkdir(parents=True, exist_ok=True)
    imp_part = f"_imp{imp}" if imp is not None else ""
    stem = f"{model_name.lower()}_{style.upper()}_0{imp_part}"
    rds_path = model_dir / f"{stem}.rds"
    pkl_path = model_dir / f"{stem}.pkl"
    W_save = W if W is not None else pd.DataFrame(index=[], columns=[])

    # Pass Z and W to R through temporary CSV files and save them as RDS
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        z_csv = td_path / "Z.csv"
        w_csv = td_path / "W.csv"
        Z.to_csv(z_csv)
        W_save.to_csv(w_csv)
        import rpy2.robjects as ro

        ro.globalenv["z_csv"] = str(z_csv)
        ro.globalenv["w_csv"] = str(w_csv)
        ro.globalenv["out_file"] = str(rds_path)
        ro.globalenv["model_name"] = model_name.lower()
        ro.globalenv["style_name"] = style.upper()
        ro.globalenv["seed_value"] = int(seed)
        ro.globalenv["n_factors_value"] = int(N_FACTORS)
        ro.globalenv["imp_value"] = ro.NA_Integer if imp is None else int(imp)
        ro.r("""
          Z <- read.csv(z_csv, row.names = 1, check.names = FALSE)
          W <- read.csv(w_csv, row.names = 1, check.names = FALSE)
          obj <- list(
            model = model_name, style = style_name, boot = 0L, imp = imp_value, seed = seed_value,
            n_factors = n_factors_value, Z = as.matrix(Z), W = as.matrix(W), fitted_model = NULL,
            rows = rownames(Z), cols = rownames(W),
            note = paste(
              "Python full-run output. Raw fitted object is saved separately",
              "as PKL when pickleable."
            )
          )
          saveRDS(obj, out_file)
        """)

    # Pickle the full fit, or just Z, W, and metadata if the fit cannot be pickled
    try:
        with open(pkl_path, "wb") as f:
            pickle.dump(
                {
                    "model": model_name.lower(), "style": style.upper(), "boot": 0, "imp": imp,
                    "seed": seed, "n_factors": N_FACTORS, "Z": Z, "W": W, "fitted_model": fit_obj,
                    "rows": list(Z.index), "cols": list(W.index) if W is not None else [],
                },
                f,
            )
    except Exception:
        with open(pkl_path, "wb") as f:
            pickle.dump(
                {
                    "model": model_name.lower(), "style": style.upper(), "boot": 0, "imp": imp,
                    "seed": seed, "n_factors": N_FACTORS, "Z": Z, "W": W, "fitted_model": None,
                    "pickle_note": "Raw fitted object could not be pickled.",
                    "rows": list(Z.index), "cols": list(W.index) if W is not None else [],
                },
                f,
            )
    return str(rds_path)


def fit_one_model(
    model_name: str, X: pd.DataFrame, layer_vec: pd.Series, diag_labels: pd.Series,
    miss_labels: pd.Series, mmse_vec: pd.Series, style: str, seed: int, imp: Optional[int] = None
) -> Dict[str, Any]:
    """Fit, save, and score one full-cohort model; k-means uses the model seed."""
    print(f"Running {model_name} | style={style} | imp={imp} | seed={seed}", flush=True)

    # Fit the model
    if model_name == "mofa":
        Z, W, fit_obj = pm.fit_mofa(X, layer_vec.reindex(X.columns), N_FACTORS, seed)
    elif model_name == "cmf":
        Z, W, fit_obj = pm.fit_cmf(
            X, layer_vec.reindex(X.columns), N_FACTORS, seed,
            lambda_=0.1, niter=30, method="als", use_item_view_info=True
        )
    else:
        raise ValueError(f"Unsupported model: {model_name}")

    # Save it and score its latent space
    model_file = save_model_rds_and_pkl(model_name, style, seed, Z, W, fit_obj, imp=imp)
    cl_diag = score_clusters_compatible(
        Z, diag_labels.reindex(X.index), X_original=X, W_loadings=W, seed=seed
    )
    cl_miss = score_clusters_compatible(
        Z, miss_labels.reindex(X.index), X_original=X, W_loadings=W, seed=seed
    )
    return {
        "model": model_name, "model_file": model_file,
        "rand_diag": cl_diag["rand"], "ari_diag": cl_diag["ari"], "chi_diag": cl_diag["chi"],
        "db_diag": cl_diag["db"], "pve_diag": cl_diag["pve"],
        "R2": pm.score_regression(Z, mmse_vec.reindex(X.index)),
        "rand_miss": cl_miss["rand"], "ari_miss": cl_miss["ari"], "chi_miss": cl_miss["chi"],
        "db_miss": cl_miss["db"], "pve_miss": cl_miss["pve"],
        "imp": imp, "style": style, "run": "full",
    }


# MAIN


def main() -> None:
    print("============================================================")
    print("PYTHON FULL RUNS: MOFA + CMF")
    print("BENCH_DIR:", BENCH_DIR)
    print("DATA_DIR:", DATA_DIR)
    print("RESULTS_FULL_DIR:", RESULTS_FULL_DIR)
    print("FULL_MODEL_DIR:", FULL_MODEL_DIR)
    print("PYTHON_METHODS_FILE:", PYTHON_METHODS_FILE)
    print("============================================================")
    RESULTS_FULL_DIR.mkdir(parents=True, exist_ok=True)
    FULL_MODEL_DIR.mkdir(parents=True, exist_ok=True)
    backup_existing_results()

    # Settings and inputs
    settings = pm.Settings(
        data_dir=DATA_DIR, out_dir=RESULTS_FULL_DIR, model_dir=FULL_MODEL_DIR, save_models=True,
        n_factors=N_FACTORS, m_imputations=M_IMPUTATIONS, imputation_seed=IMPUTATION_SEED,
        bootstrap_B=0, bootstrap_seed=100, n_jobs_bootstrap=1,
        styles_to_run=("NA", "CC_SF", "MI_SF"), models=MODELS,
        sf_frequency_threshold=SF_FREQ_THRESHOLD, use_r_bootstrap_sampling=True,
        cmf_lambda=0.1, cmf_niter=30, cmf_method="als", cmf_use_item_view_info=True
    )
    inputs = pm.load_and_process_data(settings)
    X_stand: pd.DataFrame = inputs["X_stand"]
    layer_vec: pd.Series = inputs["layer_vec"]
    diag_labels: pd.Series = inputs["diag_labels"]
    miss_labels: pd.Series = inputs["miss_labels"]
    mmse_vec: pd.Series = inputs["mmse_vec"]
    common_ids: pd.Index = inputs["common_ids"]

    # NA: full incomplete matrix
    X_na = X_stand.copy()
    rows_na = [
        fit_one_model(
            m, X_na, layer_vec, diag_labels, miss_labels, mmse_vec, "NA", SEEDS[("NA", m, None)]
        )
        for m in MODELS
    ]
    full_na_new = pd.DataFrame(rows_na)[RESULT_COLS]

    # CC_SF: complete cases with SIS screening
    X_cc = X_stand.loc[common_ids].dropna(axis=1, how="any")
    feats_cc = pm.sel_SIS(X_cc, diag_labels.reindex(X_cc.index), seed=100)
    X_cc_sf = X_cc.loc[:, feats_cc]
    rows_ccsf = [
        fit_one_model(
            m, X_cc_sf, layer_vec.reindex(X_cc_sf.columns), diag_labels, miss_labels, mmse_vec,
            "CC_SF", SEEDS[("CC_SF", m, None)]
        )
        for m in MODELS
    ]
    full_ccsf_new = pd.DataFrame(rows_ccsf)[RESULT_COLS]

    # MI_SF: SIS on each imputation, features selected in >= SF_FREQ_THRESHOLD imputations
    print("Creating/recreating five missForest imputations for MI_SF...", flush=True)
    imp_list = pm.create_imputations(X_na, M_IMPUTATIONS, IMPUTATION_SEED)
    feat_lists = [pm.sel_SIS(Xi, diag_labels.reindex(Xi.index), seed=100) for Xi in imp_list]
    freq: Dict[str, int] = {}
    for feature_list in feat_lists:
        for feature in feature_list:
            freq[feature] = freq.get(feature, 0) + 1
    common_feats = [feature for feature, count in freq.items() if count >= SF_FREQ_THRESHOLD]
    if len(common_feats) == 0:
        raise RuntimeError("MI_SF selected zero common features.")

    rows_misf_raw = []
    for imp_idx, Xi in enumerate(imp_list, start=1):
        Xi_sf = Xi.loc[:, common_feats]
        for m in MODELS:
            rows_misf_raw.append(fit_one_model(
                m, Xi_sf, layer_vec.reindex(Xi_sf.columns), diag_labels, miss_labels, mmse_vec,
                "MI_SF", SEEDS[("MI_SF", m, imp_idx)], imp=imp_idx
            ))
    full_misf_raw_new = pd.DataFrame(rows_misf_raw)[RAW_RESULT_COLS]

    # Average the imputations per model
    metric_cols = [
        "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag",
        "rand_miss", "ari_miss", "chi_miss", "db_miss", "pve_miss", "R2",
    ]
    full_misf_new = full_misf_raw_new.groupby("model", as_index=False)[metric_cols].mean(
        numeric_only=True
    )
    file_summary = full_misf_raw_new.groupby("model", as_index=False)["model_file"].agg(
        lambda x: ";".join([str(v) for v in x.dropna()])
    )
    full_misf_new = full_misf_new.merge(file_summary, on="model", how="left")
    full_misf_new["style"] = "MI_SF"
    full_misf_new["run"] = "full"
    full_misf_new = full_misf_new[RESULT_COLS]

    # Update the full-run tables
    full_na = replace_models(
        read_csv_or_empty(RESULTS_FULL_DIR / "full_NA_all.csv", RESULT_COLS), full_na_new,
        RESULT_COLS
    )
    save_csv_and_rds(
        full_na, RESULTS_FULL_DIR / "full_NA_all.csv", RESULTS_FULL_DIR / "full_NA_all.rds"
    )
    full_ccsf = replace_models(
        read_csv_or_empty(RESULTS_FULL_DIR / "full_CC_SF_all.csv", RESULT_COLS), full_ccsf_new,
        RESULT_COLS
    )
    save_csv_and_rds(
        full_ccsf, RESULTS_FULL_DIR / "full_CC_SF_all.csv", RESULTS_FULL_DIR / "full_CC_SF_all.rds"
    )
    full_misf_raw = replace_models(
        read_csv_or_empty(RESULTS_FULL_DIR / "full_MI_SF_raw_all.csv", RAW_RESULT_COLS),
        full_misf_raw_new, RAW_RESULT_COLS
    )
    save_csv_and_rds(
        full_misf_raw, RESULTS_FULL_DIR / "full_MI_SF_raw_all.csv",
        RESULTS_FULL_DIR / "full_MI_SF_raw_all.rds"
    )
    full_misf = replace_models(
        read_csv_or_empty(RESULTS_FULL_DIR / "full_MI_SF_all.csv", RESULT_COLS), full_misf_new,
        RESULT_COLS
    )
    save_csv_and_rds(
        full_misf, RESULTS_FULL_DIR / "full_MI_SF_all.csv", RESULTS_FULL_DIR / "full_MI_SF_all.rds"
    )
    full_all = pd.concat([full_na, full_ccsf, full_misf], ignore_index=True)[RESULT_COLS]
    save_csv_and_rds(
        full_all, RESULTS_FULL_DIR / "full_ALL_styles.csv", RESULTS_FULL_DIR / "full_ALL_styles.rds"
    )

    print("============================================================")
    print("UPDATED FULL-RUN FILES")
    print("============================================================")
    for path in [
        RESULTS_FULL_DIR / "full_NA_all.csv",
        RESULTS_FULL_DIR / "full_CC_SF_all.csv",
        RESULTS_FULL_DIR / "full_MI_SF_raw_all.csv",
        RESULTS_FULL_DIR / "full_MI_SF_all.csv",
        RESULTS_FULL_DIR / "full_ALL_styles.csv",
    ]:
        print(path)
    print("============================================================")
    print("MOFA/CMF rows in full_ALL_styles.csv")
    print("============================================================")
    check = pd.read_csv(RESULTS_FULL_DIR / "full_ALL_styles.csv")
    print(check.loc[check["model"].isin(MODELS), ["model", "style", "model_file"]])


if __name__ == "__main__":
    main()
