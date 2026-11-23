"""MOFA and CMF wrappers, evaluation metrics, and the ADNI observed-data bootstrap.

Clustering, SIS screening, imputation, and bootstrap resampling call R through rpy2 so that
seeds and results match the R pipeline (code/methods/adni.R). mofa_cmf_tcga.py adapts this
module to TCGA, and the semi-synthetic driver reuses its model wrappers.
"""

from __future__ import annotations

import os
import pickle
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from project_paths import input_file

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.cluster import KMeans
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LinearRegression
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    rand_score,
)

# SETTINGS

CLEAN_ROOT = next(
    p for p in Path(__file__).resolve().parents if (p / "config/default.json").is_file()
)
SCRIPT_DIR = Path(os.environ["BENCH_WORK_DIR"])
LONG_DIR = CLEAN_ROOT


@dataclass
class Settings:
    data_dir: Path = Path(os.environ["ADNI_INPUT_DIR"])
    out_dir: Path = SCRIPT_DIR / "results_boot_python"
    model_dir: Path = SCRIPT_DIR / "Bootstraps"
    save_models: bool = True
    n_factors: int = 3
    m_imputations: int = 5
    imputation_seed: int = 10
    # Number of bootstrap replicates handled by this array element
    bootstrap_B: int = int(os.environ.get("B_BOOT", "5"))
    bootstrap_start: int = int(os.environ.get("BOOT_START", "1"))
    bootstrap_seed: int = 100
    n_jobs_bootstrap: int = int(os.environ.get("LSB_DJOB_NUMPROC", "1"))
    styles_to_run: Tuple[str, ...] = ("NA", "CC_SF", "MI_SF")
    models: Tuple[str, ...] = ("mofa", "cmf")
    sf_frequency_threshold: int = 4
    use_r_bootstrap_sampling: bool = True
    cmf_lambda: float = 0.1
    cmf_niter: int = 30
    cmf_method: str = "als"
    cmf_use_item_view_info: bool = True


def r_kmeans_labels(latent_matrix: np.ndarray, n_clusters: int, seed: int) -> np.ndarray:
    """Cluster with R's kmeans (25 starts) after set.seed(seed), as in the R pipeline."""
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.conversion import localconverter

    Z = np.asarray(latent_matrix, dtype=float)
    ro.r["set.seed"](int(seed))
    with localconverter(ro.default_converter + numpy2ri.converter):
        r_Z = ro.conversion.py2rpy(Z)
    fit = ro.r["kmeans"](r_Z, centers=int(n_clusters), nstart=25)
    return np.asarray(fit.rx2("cluster"), dtype=int)


# DATA LOADING


def standardize_with_nan(X: pd.DataFrame) -> pd.DataFrame:
    """Standardize columns using observed values only; missing cells stay missing."""
    means = X.mean(axis=0, skipna=True)
    sds = X.std(axis=0, skipna=True, ddof=1).replace(0.0, np.nan)
    return (X - means) / sds


def load_and_process_data(settings: Settings) -> Dict[str, object]:
    """Load ADNI inputs as in code/preprocessing/load_prepared_adni.R (top 500 genes)."""

    # Read the inputs, indexed by RID
    imeta_raw = pd.read_csv(settings.data_dir / "meta.csv")
    igene_raw = pd.read_csv(input_file(settings.data_dir, "rnaseq.csv"))
    clin = pd.read_csv(settings.data_dir / "clin.csv")
    for df in (imeta_raw, igene_raw, clin):
        df["RID"] = df["RID"].astype(str)
    imeta_raw = imeta_raw.set_index("RID")
    igene_raw = igene_raw.set_index("RID")
    clin = clin.set_index("RID")

    # Retain the 500 gene-expression features with the highest variance
    top_gene_cols = (
        igene_raw.var(axis=0, skipna=True, ddof=1).sort_values(ascending=False).head(500).index
    )
    igene_raw = igene_raw.loc[:, top_gene_cols]

    # Align patients across modalities and standardize
    ids_gene, ids_meta = igene_raw.index, imeta_raw.index
    all_ids = ids_gene.union(ids_meta, sort=False)
    igene = igene_raw.reindex(all_ids)
    imeta = imeta_raw.reindex(all_ids)
    X_all = pd.concat([igene, imeta], axis=1)
    X_stand = standardize_with_nan(X_all)
    layer_vec = pd.Series(
        ["gene"] * igene.shape[1] + ["meta"] * imeta.shape[1], index=X_stand.columns, name="layer"
    )

    # Label patients by observed modality availability
    missingness_label = pd.Series(index=all_ids, dtype="object", name="missingness")
    common = ids_gene.intersection(ids_meta)
    missingness_label.loc[common] = "CC"
    missingness_label.loc[ids_gene.difference(ids_meta)] = "missing_meta"
    missingness_label.loc[ids_meta.difference(ids_gene)] = "missing_gene"
    missingness_label = pd.Series(
        pd.Categorical(missingness_label, categories=["CC", "missing_meta", "missing_gene"]),
        index=all_ids, name="missingness"
    )

    # Restrict to patients with clinical data and drop empty features
    diag_labels = clin["diagnosis"].copy()
    mmse_vec = clin["MMSE"].copy()
    ids_score = all_ids.intersection(diag_labels.index, sort=False)
    X_stand = X_stand.reindex(ids_score)
    diag_labels = diag_labels.reindex(ids_score)
    mmse_vec = mmse_vec.reindex(ids_score)
    missingness_label = missingness_label.reindex(ids_score)
    common_ids = common.intersection(ids_score, sort=False)

    keep_cols = X_stand.notna().any(axis=0)
    X_stand = X_stand.loc[:, keep_cols]
    layer_vec = layer_vec.reindex(X_stand.columns)
    return {
        "X_stand": X_stand, "layer_vec": layer_vec, "diag_labels": diag_labels,
        "miss_labels": missingness_label, "mmse_vec": mmse_vec, "common_ids": common_ids
    }


# EVALUATION METRICS


def score_clusters(
    latent_matrix: pd.DataFrame, true_labels: pd.Series, X_original: Optional[pd.DataFrame] = None,
    W_loadings: Optional[pd.DataFrame] = None, seed: int = 100
) -> Dict[str, float]:
    """k-means agreement with the labels (Rand, ARI) plus internal criteria and PVE.

    Patients with missing labels are excluded. PVE is the variance of Z @ W.T relative to the
    observed variance of X.
    """

    # Align labels and drop missing ones
    labels = true_labels.reindex(latent_matrix.index)
    keep = labels.notna()
    Z = latent_matrix.loc[keep]
    labels = labels.loc[keep]
    k = int(labels.nunique())
    empty = {"rand": np.nan, "ari": np.nan, "chi": np.nan, "db": np.nan, "pve": np.nan}
    if k < 2 or len(labels) <= k:
        return empty

    # Cluster and compare with the labels
    pred = r_kmeans_labels(latent_matrix=Z.to_numpy(), n_clusters=k, seed=seed)
    out = {
        "rand": rand_score(labels, pred), "ari": adjusted_rand_score(labels, pred),
        "chi": calinski_harabasz_score(Z.to_numpy(), pred),
        "db": davies_bouldin_score(Z.to_numpy(), pred), "pve": np.nan
    }

    # Variance of the reconstruction relative to the observed variance
    if X_original is not None and W_loadings is not None:
        feats = X_original.columns.intersection(W_loadings.index, sort=False)
        if len(feats) > 0:
            X_use = X_original.loc[Z.index, feats].to_numpy(dtype=float)
            W_use = W_loadings.loc[feats].to_numpy(dtype=float)
            X_hat = Z.to_numpy(dtype=float) @ W_use.T
            var_x = np.nansum(np.nanvar(X_use, axis=0, ddof=1))
            var_hat = np.nansum(np.var(X_hat, axis=0, ddof=1))
            if np.isfinite(var_x) and var_x > 0:
                out["pve"] = float(var_hat / var_x)
    return out


def score_regression(latent_matrix: pd.DataFrame, mmse_vec: Optional[pd.Series]) -> float:
    """In-sample R^2 of the clinical score regressed on all latent factors."""
    if mmse_vec is None:
        return np.nan
    y = mmse_vec.reindex(latent_matrix.index)
    keep = y.notna()
    if keep.sum() <= latent_matrix.shape[1] + 1:
        return np.nan
    model = LinearRegression().fit(latent_matrix.loc[keep].to_numpy(), y.loc[keep].to_numpy())
    return float(model.score(latent_matrix.loc[keep].to_numpy(), y.loc[keep].to_numpy()))


# FEATURE SCREENING (*_SF STRATEGIES)


def sel_SIS(X_complete: pd.DataFrame, label_vec: pd.Series, seed: int = 100) -> List[str]:
    """Binary sure independence screening (binomial, BIC tuning) through the R SIS package.

    ADNI codes AD as 1 and every other diagnosis as 0. For TCGA (no AD label, exactly two
    disease groups) the second sorted level is coded 1. This keeps the original ADNI behavior
    while letting the same bootstrap code run on the binary TCGA disease_code labels.
    """
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.vectors import IntVector
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    # Complete rows with a label
    y_raw = label_vec.reindex(X_complete.index)
    keep = y_raw.notna()
    X_use = X_complete.loc[keep].copy()
    labels = y_raw.loc[keep].astype(str)
    if X_use.isna().any().any():
        raise ValueError("SIS requires a completed matrix.")
    if X_use.shape[1] == 0:
        raise ValueError("SIS received zero features.")

    levels = sorted(pd.unique(labels).tolist())
    if "AD" in levels:
        # ADNI: AD versus all other diagnoses
        ybin = (labels == "AD").astype(int).to_numpy()
        if np.unique(ybin).size < 2:
            raise ValueError("SIS requires both AD and non-AD labels.")
    elif len(levels) == 2:
        # TCGA: two disease groups
        ybin = (labels == levels[1]).astype(int).to_numpy()
        if np.unique(ybin).size < 2:
            raise ValueError("SIS requires both TCGA disease groups.")
    else:
        raise ValueError(
            "SIS requires either ADNI labels containing AD "
            "or exactly two outcome classes. "
            f"Observed labels: {levels}"
        )

    # Run SIS in R
    SIS_pkg = importr("SIS")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_X_df = ro.conversion.py2rpy(X_use)
    r_y = IntVector(ybin.astype(int).tolist())
    ro.r["set.seed"](int(seed))
    sis_fit = SIS_pkg.SIS(
        x=ro.r["as.matrix"](r_X_df), y=r_y, family="binomial", tune="bic", nsis=int(X_use.shape[1])
    )
    selected_r = np.array(sis_fit.rx2("sis.ix0"), dtype=int)
    if selected_r.size == 0:
        raise ValueError("SIS selected zero features.")
    return X_use.columns[selected_r - 1].tolist()


# IMPUTATION


def impute_once(X_na: pd.DataFrame, seed: int) -> pd.DataFrame:
    """One missForest imputation in R after set.seed(seed)."""
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    # Seed, convert to an R data frame, impute, and convert back
    missForest = importr("missForest")
    ro.r["set.seed"](int(seed))
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_X = ro.conversion.py2rpy(X_na)
    fit = missForest.missForest(r_X, verbose=False)
    with localconverter(ro.default_converter + pandas2ri.converter):
        out = ro.conversion.rpy2py(fit.rx2("ximp"))
    return pd.DataFrame(out, index=X_na.index, columns=X_na.columns).astype(float)


def create_imputations(X_na: pd.DataFrame, M: int, seed: int) -> List[pd.DataFrame]:
    """M imputations seeded seed + 1, ..., seed + M, matching prepare_adni.R."""
    return [impute_once(X_na, seed + m + 1) for m in range(M)]


def load_or_create_imputations(X_na: pd.DataFrame, settings: Settings) -> List[pd.DataFrame]:
    """Reuse the imputations cached by the R preparation step, or create and cache them."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.conversion import localconverter

    # Cache file, preferring a copy in CACHE_INPUT_DIR
    var_dir = SCRIPT_DIR / "variables"
    var_dir.mkdir(parents=True, exist_ok=True)
    imp_file = var_dir / f"imp_list_M{settings.m_imputations}_seed{settings.imputation_seed}.rds"
    external = Path(os.environ.get("CACHE_INPUT_DIR") or var_dir) / imp_file.name
    if external.exists():
        imp_file = external

    # Read the cached list of imputed matrices
    if imp_file.exists():
        print(f"Loading saved imp_list from: {imp_file}", flush=True)
        r_list = ro.r["readRDS"](str(imp_file))
        imp_list: List[pd.DataFrame] = []
        for i in range(len(r_list)):
            with localconverter(ro.default_converter + numpy2ri.converter):
                arr = ro.conversion.rpy2py(r_list[i])
            imp_list.append(
                pd.DataFrame(np.asarray(arr, dtype=float), index=X_na.index, columns=X_na.columns)
            )
        return imp_list

    # Otherwise impute now and save the list as RDS
    print(f"Computing imp_list and saving to: {imp_file}", flush=True)
    imp_list = create_imputations(X_na, settings.m_imputations, settings.imputation_seed)
    r_objects = {}
    for i, Xi in enumerate(imp_list, start=1):
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_df = ro.conversion.py2rpy(Xi)
        r_objects[str(i)] = ro.r["as.matrix"](r_df)
    r_list = ro.ListVector(r_objects)
    ro.r["saveRDS"](r_list, str(imp_file))
    return imp_list


# MODELS


def get_view_columns(
    X: pd.DataFrame, layer_vec: pd.Series, require_observed: bool = False
) -> Dict[str, List[str]]:
    """Split the columns of X into omics views using layer_vec."""

    # Keep views with at least one column (and, if required, one observed value)
    layer_use = layer_vec.reindex(X.columns)
    view_columns: Dict[str, List[str]] = {}
    for view in pd.unique(layer_use.dropna()):
        cols = X.columns[layer_use.to_numpy() == view].tolist()
        if cols and (
            not require_observed or np.isfinite(X.loc[:, cols].to_numpy(dtype=float)).sum() > 0
        ):
            view_columns[str(view)] = cols
    if not view_columns:
        raise ValueError("No usable omics views found.")
    return view_columns


# Orient a MOFA output array as samples x factors
def _as_samples_by_factors(arr: np.ndarray, n_samples: int, n_factors: int) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 3:
        arr = arr[0]
    if arr.shape[0] == n_samples:
        return arr
    if arr.shape[1] == n_samples:
        return arr.T
    if arr.shape[0] == n_factors:
        return arr.T
    return arr


# Orient a MOFA output array as features x factors
def _as_features_by_factors(arr: np.ndarray, n_features: int, n_factors: int) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.shape[0] == n_features:
        return arr
    if arr.shape[1] == n_features:
        return arr.T
    if arr.shape[0] == n_factors:
        return arr.T
    return arr


def fit_mofa(
    X: pd.DataFrame, layer_vec: pd.Series, n_factors: int, seed: int
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Any]:
    """MOFA (mofapy2) with Gaussian views and view/group scaling disabled; NaNs stay missing."""
    from mofapy2.run.entry_point import entry_point

    # One matrix per view with missing cells left as NaN
    view_columns = get_view_columns(X, layer_vec, require_observed=True)
    views = list(view_columns.keys())
    data = [[X.loc[:, view_columns[view]].to_numpy(dtype=float)] for view in views]

    # Build and train the model
    ent = entry_point()
    ent.set_data_options(scale_views=False, scale_groups=False, center_groups=False)
    ent.set_data_matrix(
        data=data, likelihoods=["gaussian"] * len(views), views_names=views,
        groups_names=["group1"], samples_names=[X.index.astype(str).tolist()],
        features_names=[view_columns[view] for view in views]
    )
    ent.set_model_options(factors=int(n_factors))
    ent.set_train_options(seed=int(seed), verbose=False, quiet=True, gpu_mode=False)
    ent.build()
    ent.run()

    # Posterior mean factors and weights
    Z_raw = ent.model.nodes["Z"].getExpectation()
    Z_arr = _as_samples_by_factors(np.asarray(Z_raw), n_samples=X.shape[0], n_factors=n_factors)
    Z = pd.DataFrame(
        Z_arr[:, : min(n_factors, Z_arr.shape[1])], index=X.index,
        columns=[f"Factor{i + 1}" for i in range(min(n_factors, Z_arr.shape[1]))]
    )

    # Stack the view weights in feature order; W stays None if they cannot be read
    W: Optional[pd.DataFrame] = None
    try:
        W_raw = ent.model.nodes["W"].getExpectation()
        W_parts = []
        for i, view in enumerate(views):
            W_i = _as_features_by_factors(
                np.asarray(W_raw[i]), n_features=len(view_columns[view]), n_factors=Z.shape[1]
            )
            W_parts.append(
                pd.DataFrame(W_i[:, : Z.shape[1]], index=view_columns[view], columns=Z.columns)
            )
        W = pd.concat(W_parts, axis=0).reindex(X.columns)
    except Exception:
        W = None
    return Z, W, ent


def _make_cmf_item_info(X: pd.DataFrame, layer_vec: pd.Series) -> Optional[np.ndarray]:
    """One-hot view membership of each feature, used as CMF item side information."""

    layers = layer_vec.reindex(X.columns).astype(str)
    unique_layers = list(pd.unique(layers))
    if len(unique_layers) <= 1:
        return None
    layer_to_idx = {layer: i for i, layer in enumerate(unique_layers)}

    I = np.zeros((X.shape[1], len(unique_layers)), dtype=float)
    for j, feature in enumerate(X.columns):
        I[j, layer_to_idx[layers.loc[feature]]] = 1.0
    return I


def fit_cmf(
    X: pd.DataFrame, layer_vec: pd.Series, n_factors: int, seed: int, lambda_: float = 0.1,
    niter: int = 30, method: str = "als", use_item_view_info: bool = True
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Any]:
    """Collective matrix factorization (cmfrec) without biases or centering; NaNs are unobserved."""
    from cmfrec import CMF

    # Every row and feature needs at least one observed value
    X_arr = X.to_numpy(dtype=float)
    if not np.isfinite(X_arr).any():
        raise ValueError("CMF received a matrix with no observed values.")
    observed_by_row = np.isfinite(X_arr).sum(axis=1)
    observed_by_col = np.isfinite(X_arr).sum(axis=0)
    if np.any(observed_by_row == 0):
        raise ValueError("CMF: at least one bootstrap row has no observed entries.")
    if np.any(observed_by_col == 0):
        raise ValueError("CMF: at least one feature has no observed entries.")

    # Fit CMF with view membership as item side information
    I = _make_cmf_item_info(X, layer_vec) if use_item_view_info else None
    kwargs = dict(
        k=int(n_factors), lambda_=float(lambda_), method=str(method), niter=int(niter),
        user_bias=False, item_bias=False, center=False, center_U=False, center_I=False,
        random_state=int(seed), verbose=False, nthreads=1, produce_dicts=False
    )
    try:
        fit = CMF(**kwargs)
    except TypeError:
        # Older cmfrec versions use different argument names for centering and biases
        kwargs.pop("center", None)
        kwargs["center_ratings"] = False
        kwargs["add_user_bias"] = False
        kwargs["add_item_bias"] = False
        fit = CMF(**kwargs)
    fit.fit(X_arr, I=I)

    # Sample factors become Z, feature factors become W
    A = np.asarray(fit.A_, dtype=float)
    B = np.asarray(fit.B_, dtype=float)
    if A.shape[0] != X.shape[0]:
        raise ValueError(f"CMF: expected {X.shape[0]} sample factors but got {A.shape[0]}.")
    if B.shape[0] != X.shape[1]:
        raise ValueError(f"CMF: expected {X.shape[1]} feature factors but got {B.shape[0]}.")
    k = min(int(n_factors), A.shape[1], B.shape[1])
    if k < 1:
        raise ValueError("CMF returned zero usable factors.")
    factor_names = [f"Factor{i + 1}" for i in range(k)]
    Z = pd.DataFrame(A[:, :k], index=X.index, columns=factor_names)
    W = pd.DataFrame(B[:, :k], index=X.columns, columns=factor_names)
    return Z, W, fit


# RUNNING MODELS


def save_boot_model(
    obj: Dict[str, Any], model_name: str, style: str, boot: int, settings: Settings,
    imp: Optional[int] = None
) -> str:
    model_name = model_name.lower()
    style = style.upper()
    model_dir = settings.model_dir / model_name
    model_dir.mkdir(parents=True, exist_ok=True)
    imp_part = f"_imp{imp}" if imp is not None else ""
    path = model_dir / f"{model_name}_{style}_{boot}{imp_part}.pkl"
    try:
        with open(path, "wb") as f:
            pickle.dump(obj, f)
    except Exception:
        # Some package model objects do not pickle; keep Z, W, and metadata so the fit can
        # still be inspected
        obj = dict(obj)
        obj["fitted_model"] = None
        obj["pickle_note"] = "Raw fitted object could not be pickled; saved Z/W/metadata only."
        with open(path, "wb") as f:
            pickle.dump(obj, f)
    return str(path)


def run_model(
    model_name: str, X: pd.DataFrame, layer_vec: pd.Series, y_diag: pd.Series,
    y_miss: Optional[pd.Series], mmse_vec: Optional[pd.Series], n_factors: int, model_seed: int,
    settings: Settings, style: str, boot: int, imp: Optional[int] = None
) -> Dict[str, object]:
    """Fit MOFA or CMF and score it; k-means seeds are 100 + boot (diagnosis), 200 + boot."""
    model_name = model_name.lower()
    if model_name == "mofa":
        Z, W, fit_obj = fit_mofa(X, layer_vec, n_factors, model_seed)
    elif model_name == "cmf":
        Z, W, fit_obj = fit_cmf(
            X, layer_vec, n_factors, model_seed, settings.cmf_lambda, settings.cmf_niter,
            settings.cmf_method, settings.cmf_use_item_view_info
        )
    else:
        raise ValueError(f"Unsupported model: {model_name}. This pipeline runs only mofa and cmf.")

    # Optionally save the fit
    model_file: Optional[str] = None
    if settings.save_models:
        model_file = save_boot_model(
            obj={
                "model": model_name, "style": style.upper(), "boot": boot, "imp": imp,
                "seed": model_seed, "n_factors": n_factors, "Z": Z, "W": W,
                "fitted_model": fit_obj, "rows": list(X.index), "cols": list(X.columns),
            },
            model_name=model_name, style=style, boot=boot, imp=imp, settings=settings
        )

    # Score the latent space
    cl_diag = score_clusters(Z, y_diag, X_original=X, W_loadings=W, seed=100 + boot)

    result: Dict[str, object] = {
        "model": model_name, "model_file": model_file, "rand_diag": cl_diag["rand"],
        "ari_diag": cl_diag["ari"], "chi_diag": cl_diag["chi"], "db_diag": cl_diag["db"],
        "pve_diag": cl_diag["pve"], "R2": score_regression(Z, mmse_vec)
    }

    # Missingness-category clustering, when labels are given
    if y_miss is not None:
        cl_miss = score_clusters(Z, y_miss, X_original=X, W_loadings=W, seed=200 + boot)
        result.update({
            "rand_miss": cl_miss["rand"], "ari_miss": cl_miss["ari"], "chi_miss": cl_miss["chi"],
            "db_miss": cl_miss["db"], "pve_miss": cl_miss["pve"],
        })
    return result


def run_complete(
    X: pd.DataFrame, layer_vec: pd.Series, diag_labels: pd.Series, miss_labels: pd.Series,
    mmse_vec: pd.Series, models: Sequence[str], n_factors: int, model_seed: int, settings: Settings,
    style: str, boot: int, imp: Optional[int] = None
) -> pd.DataFrame:
    # One result row per model
    rows = [
        run_model(
            model, X, layer_vec.reindex(X.columns), diag_labels.reindex(X.index),
            miss_labels.reindex(X.index), mmse_vec.reindex(X.index), n_factors, model_seed,
            settings, style=style, boot=boot, imp=imp
        )
        for model in models
    ]
    return pd.DataFrame(rows)


def run_NA(
    X_na: pd.DataFrame, layer_vec: pd.Series, diag_labels: pd.Series, miss_labels: pd.Series,
    mmse_vec: pd.Series, models: Sequence[str], n_factors: int, model_seed: int, settings: Settings,
    style: str, boot: int
) -> pd.DataFrame:
    """Fit methods that accept missing values on the incomplete matrix (NA strategy)."""
    return run_complete(
        X_na, layer_vec, diag_labels, miss_labels, mmse_vec, models, n_factors, model_seed,
        settings, style=style, boot=boot
    )


METRIC_COLUMNS = [
    "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag",
    "rand_miss", "ari_miss", "chi_miss", "db_miss", "pve_miss", "R2",
]


def run_multiple_complete(
    imp_list: Sequence[pd.DataFrame], layer_vec: pd.Series, fs_diag_labels: pd.Series,
    cluster_diag_labels: pd.Series, cluster_miss_labels: pd.Series, mmse_vec: pd.Series,
    models: Sequence[str], fs: bool, n_factors: int, model_seed: int, settings: Settings,
    style: str, boot: int
) -> pd.DataFrame:
    """Fit methods on each imputation.

    With fs, SIS runs on every imputation and only features selected in at least
    settings.sf_frequency_threshold imputations are kept.
    """

    common_feats: Optional[List[str]] = None

    # Count how often each feature is selected across imputations
    if fs:
        feat_lists = [
            sel_SIS(Xi, fs_diag_labels.reindex(Xi.index), seed=model_seed) for Xi in imp_list
        ]
        freq: Dict[str, int] = {}
        for feature_list in feat_lists:
            for feature in feature_list:
                freq[feature] = freq.get(feature, 0) + 1
        threshold = settings.sf_frequency_threshold
        common_feats = [feature for feature, count in freq.items() if count >= threshold]
        if not common_feats:
            raise ValueError(f"MI_SF: no features met frequency threshold (>={threshold}).")

    # Fit the methods on each imputation
    rows: List[pd.DataFrame] = []

    for i, Xi in enumerate(imp_list, start=1):
        Xi_use = Xi.loc[:, common_feats] if fs else Xi
        res = run_complete(
            Xi_use, layer_vec.reindex(Xi_use.columns), cluster_diag_labels, cluster_miss_labels,
            mmse_vec, models, n_factors, model_seed, settings, style=style, boot=boot, imp=i
        )
        res["imp"] = i
        rows.append(res)
    return pd.concat(rows, ignore_index=True)


# BOOTSTRAP


def _bootstrap_series(
    values: pd.Series, sampled_ids: Sequence[str], new_index: pd.Index
) -> pd.Series:
    return pd.Series(values.reindex(sampled_ids).to_numpy(), index=new_index, name=values.name)


def sample_ids_for_bootstrap(ids_all: Sequence[str], seed: int, use_r: bool = True) -> np.ndarray:
    """Resample patient IDs with R's set.seed + sample() so draws match the R bootstrap."""
    if not use_r:
        rng = np.random.default_rng(seed)
        return rng.choice(np.asarray(ids_all), size=len(ids_all), replace=True)

    import rpy2.robjects as ro
    from rpy2.robjects.vectors import StrVector

    ro.r["set.seed"](int(seed))
    sampled = ro.r["sample"](StrVector([str(x) for x in ids_all]), int(len(ids_all)), replace=True)
    return np.asarray(list(sampled), dtype=object)


def _run_one_boot(
    b: int, X_stand: pd.DataFrame, layer_vec: pd.Series, diag_labels: pd.Series,
    miss_labels: pd.Series, mmse_vec: pd.Series, common_ids: pd.Index, style: str,
    models: Sequence[str], n_factors: int, seed: int, settings: Settings,
    imputed_list: Optional[Sequence[pd.DataFrame]], single_imputed: Optional[pd.DataFrame]
) -> pd.DataFrame:
    """Replicate b: resample with seed + b, then apply one strategy (SIS is rerun per replicate)."""

    # Resample patients and carry labels along
    ids_all = X_stand.index.astype(str).to_numpy()
    sampled_ids = sample_ids_for_bootstrap(
        ids_all, seed + b, use_r=settings.use_r_bootstrap_sampling
    )
    new_index = pd.Index([f"{rid}_b{b}_{i + 1}" for i, rid in enumerate(sampled_ids)])
    X_b = X_stand.reindex(sampled_ids).copy()
    X_b.index = new_index
    diag_b = _bootstrap_series(diag_labels, sampled_ids, new_index)
    miss_b = _bootstrap_series(miss_labels, sampled_ids, new_index)
    mmse_b = _bootstrap_series(mmse_vec, sampled_ids, new_index)

    fs_diag_b = diag_b.copy()
    is_cc = np.isin(sampled_ids, common_ids.to_numpy())
    style = style.upper()

    # Fit the requested strategy on the replicate
    if style == "CC":
        # Complete cases with fully observed features
        X_use = X_b.loc[is_cc].dropna(axis=1, how="any")
        res = run_complete(
            X_use, layer_vec, diag_b.loc[X_use.index], miss_b.loc[X_use.index],
            mmse_b.loc[X_use.index], models, n_factors, seed + b, settings, style=style, boot=b
        )
    elif style == "CC_SF":
        # Complete cases followed by SIS screening
        X_cc = X_b.loc[is_cc].dropna(axis=1, how="any")
        selected = sel_SIS(X_cc, fs_diag_b.reindex(X_cc.index), seed=seed + b)
        X_use = X_cc.loc[:, selected]
        res = run_complete(
            X_use, layer_vec.reindex(selected), diag_b.loc[X_use.index], miss_b.loc[X_use.index],
            mmse_b.loc[X_use.index], models, n_factors, seed + b, settings, style=style, boot=b
        )
    elif style == "SI":
        # Resample the precomputed single imputation
        if single_imputed is None:
            raise ValueError("SI requires a precomputed single imputation.")
        X_use = single_imputed.reindex(sampled_ids).copy()
        X_use.index = new_index
        res = run_complete(
            X_use, layer_vec, diag_b, miss_b, mmse_b, models, n_factors, seed + b, settings,
            style=style, boot=b
        )
    elif style == "SI_SF":
        # Resample the precomputed single imputation, then screen
        if single_imputed is None:
            raise ValueError("SI_SF requires a precomputed single imputation.")
        X_si = single_imputed.reindex(sampled_ids).copy()
        X_si.index = new_index
        selected = sel_SIS(X_si, fs_diag_b, seed=seed + b)
        X_use = X_si.loc[:, selected]
        res = run_complete(
            X_use, layer_vec.reindex(selected), diag_b, miss_b, mmse_b, models, n_factors, seed + b,
            settings, style=style, boot=b
        )
    elif style in {"MI", "MI_SF"}:
        # Resample every imputation with the same patient draw and average over imputations
        if imputed_list is None:
            raise ValueError(f"{style} requires precomputed multiple imputations.")
        imp_list_b: List[pd.DataFrame] = []
        for Xi in imputed_list:
            Xi_b = Xi.reindex(sampled_ids).copy()
            Xi_b.index = new_index
            imp_list_b.append(Xi_b)
        raw = run_multiple_complete(
            imp_list_b, layer_vec, fs_diag_b, diag_b, miss_b, mmse_b, models,
            fs=(style == "MI_SF"), n_factors=n_factors, model_seed=seed + b, settings=settings,
            style=style, boot=b
        )

        res = raw.groupby("model", as_index=False)[METRIC_COLUMNS].mean(numeric_only=True)
        if "model_file" in raw.columns:
            file_summary = raw.groupby("model", as_index=False)["model_file"].agg(
                lambda x: ";".join([str(v) for v in x.dropna()])
            )
            res = res.merge(file_summary, on="model", how="left")
    elif style == "NA":
        # Incomplete matrix for methods that accept missing values
        res = run_NA(
            X_b, layer_vec, diag_b, miss_b, mmse_b, models, n_factors, seed + b, settings,
            style=style, boot=b
        )
    else:
        raise ValueError(f"Unknown style: {style}")

    res["boot"] = b
    return res


def bootstrap_pipe(
    X_stand: pd.DataFrame, layer_vec: pd.Series, diag_labels: pd.Series, miss_labels: pd.Series,
    mmse_vec: pd.Series, common_ids: pd.Index, style: str, B: int, models: Sequence[str],
    n_factors: int, seed: int, settings: Settings,
    imputed_list: Optional[Sequence[pd.DataFrame]] = None, n_jobs: int = 1, verbose: bool = True,
    boot_start: int = 1
) -> pd.DataFrame:
    """Run replicates boot_start, ..., boot_start + B - 1 for one strategy."""

    style = style.upper()
    single_imputed = None
    if style in {"SI", "SI_SF"}:
        single_imputed = (
            imputed_list[0]
            if imputed_list is not None
            else impute_once(X_stand, settings.imputation_seed)
        )

    boot_ids = range(boot_start, boot_start + B)

    if verbose:
        mode = "sequential" if n_jobs == 1 else f"parallel ({n_jobs} workers)"
        print(
            f"Starting {style} bootstrap for {list(models)}: {mode}; "
            f"boots {boot_start}-{boot_start + B - 1}",
            flush=True,
        )

    # Run the replicates sequentially or in parallel
    args = dict(
        X_stand=X_stand, layer_vec=layer_vec, diag_labels=diag_labels, miss_labels=miss_labels,
        mmse_vec=mmse_vec, common_ids=common_ids, style=style, models=models, n_factors=n_factors,
        seed=seed, settings=settings, imputed_list=imputed_list, single_imputed=single_imputed
    )
    if n_jobs == 1:
        results = [_run_one_boot(b=b, **args) for b in boot_ids]
    else:
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_run_one_boot)(b=b, **args) for b in boot_ids
        )
    return pd.concat(results, ignore_index=True)


# BENCHMARK DRIVERS


def run_benchmark(settings: Settings) -> pd.DataFrame:
    """Run one chunk of bootstrap replicates for every strategy and save the results."""

    inputs = load_and_process_data(settings)
    boot_start = settings.bootstrap_start
    boot_end = boot_start + settings.bootstrap_B - 1
    boot_tag = f"{boot_start:04d}-{boot_end:04d}"
    print(
        f"Running Python bootstrap chunk {boot_tag} with {settings.n_jobs_bootstrap} workers.",
        flush=True,
    )

    X_stand: pd.DataFrame = inputs["X_stand"]
    layer_vec: pd.Series = inputs["layer_vec"]
    diag_labels: pd.Series = inputs["diag_labels"]
    miss_labels: pd.Series = inputs["miss_labels"]
    mmse_vec: pd.Series = inputs["mmse_vec"]
    common_ids: pd.Index = inputs["common_ids"]
    settings.out_dir.mkdir(parents=True, exist_ok=True)
    settings.model_dir.mkdir(parents=True, exist_ok=True)

    needs_imputation = any(
        s.upper() in {"MI", "MI_SF", "SI", "SI_SF"} for s in settings.styles_to_run
    )
    imp_list = load_or_create_imputations(X_stand, settings) if needs_imputation else None

    # Run each strategy and save its results
    style_outputs: List[pd.DataFrame] = []
    for style in settings.styles_to_run:
        result = bootstrap_pipe(
            X_stand=X_stand, layer_vec=layer_vec, diag_labels=diag_labels, miss_labels=miss_labels,
            mmse_vec=mmse_vec, common_ids=common_ids, style=style, B=settings.bootstrap_B,
            models=settings.models, n_factors=settings.n_factors, seed=settings.bootstrap_seed,
            settings=settings, imputed_list=imp_list, n_jobs=settings.n_jobs_bootstrap,
            verbose=True, boot_start=settings.bootstrap_start
        )
        result["style"] = style.upper()
        result.to_csv(
            settings.out_dir / f"boot_{style.upper()}_mofa_cmf_{boot_tag}.csv", index=False
        )
        result.to_pickle(settings.out_dir / f"boot_{style.upper()}_mofa_cmf_{boot_tag}.pkl")
        style_outputs.append(result)

    boot_all = pd.concat(style_outputs, ignore_index=True)
    boot_all.to_csv(settings.out_dir / f"boot_ALL_styles_mofa_cmf_{boot_tag}.csv", index=False)
    boot_all.to_pickle(settings.out_dir / f"boot_ALL_styles_mofa_cmf_{boot_tag}.pkl")
    print(boot_all)
    return boot_all


def run_full_benchmark(settings: Settings) -> pd.DataFrame:
    """Fit each strategy once on the full cohort (boot = 0)."""
    inputs = load_and_process_data(settings)

    X_stand: pd.DataFrame = inputs["X_stand"]
    layer_vec: pd.Series = inputs["layer_vec"]
    diag_labels: pd.Series = inputs["diag_labels"]
    miss_labels: pd.Series = inputs["miss_labels"]
    mmse_vec: pd.Series = inputs["mmse_vec"]
    common_ids: pd.Index = inputs["common_ids"]
    settings.out_dir.mkdir(parents=True, exist_ok=True)
    settings.model_dir.mkdir(parents=True, exist_ok=True)

    model_seed = settings.bootstrap_seed
    full_outputs: List[pd.DataFrame] = []

    # NA: the full incomplete matrix
    if "NA" in [s.upper() for s in settings.styles_to_run]:
        full_na = run_NA(
            X_na=X_stand, layer_vec=layer_vec, diag_labels=diag_labels, miss_labels=miss_labels,
            mmse_vec=mmse_vec, models=settings.models, n_factors=settings.n_factors,
            model_seed=model_seed, settings=settings, style="NA", boot=0
        )
        full_na["style"] = "NA"
        full_na["run"] = "full"
        full_na["boot"] = 0
        full_na.to_csv(settings.out_dir / "full_NA_mofa_cmf.csv", index=False)
        full_na.to_pickle(settings.out_dir / "full_NA_mofa_cmf.pkl")
        full_outputs.append(full_na)

    # CC_SF: complete cases, with feature screening run once
    if "CC_SF" in [s.upper() for s in settings.styles_to_run]:
        X_cc = X_stand.loc[common_ids].dropna(axis=1, how="any")
        selected = sel_SIS(X_cc, diag_labels.reindex(X_cc.index), seed=model_seed)
        X_cc_sf = X_cc.loc[:, selected]
        full_cc_sf = run_complete(
            X=X_cc_sf, layer_vec=layer_vec.reindex(selected), diag_labels=diag_labels,
            miss_labels=miss_labels, mmse_vec=mmse_vec, models=settings.models,
            n_factors=settings.n_factors, model_seed=model_seed, settings=settings, style="CC_SF",
            boot=0
        )
        full_cc_sf["style"] = "CC_SF"
        full_cc_sf["run"] = "full"
        full_cc_sf["boot"] = 0
        full_cc_sf.to_csv(settings.out_dir / "full_CC_SF_mofa_cmf.csv", index=False)
        full_cc_sf.to_pickle(settings.out_dir / "full_CC_SF_mofa_cmf.pkl")
        full_outputs.append(full_cc_sf)

    # MI_SF: five imputations of the full matrix, with feature screening across imputations
    if "MI_SF" in [s.upper() for s in settings.styles_to_run]:
        imp_list = load_or_create_imputations(X_stand, settings)
        full_mi_sf_raw = run_multiple_complete(
            imp_list=imp_list, layer_vec=layer_vec, fs_diag_labels=diag_labels,
            cluster_diag_labels=diag_labels, cluster_miss_labels=miss_labels, mmse_vec=mmse_vec,
            models=settings.models, fs=True, n_factors=settings.n_factors, model_seed=model_seed,
            settings=settings, style="MI_SF", boot=0
        )
        full_mi_sf_raw["style"] = "MI_SF"
        full_mi_sf_raw["run"] = "full"
        full_mi_sf_raw["boot"] = 0

        full_mi_sf = full_mi_sf_raw.groupby("model", as_index=False)[METRIC_COLUMNS].mean(
            numeric_only=True
        )

        file_summary = full_mi_sf_raw.groupby("model", as_index=False)["model_file"].agg(
            lambda x: ";".join([str(v) for v in x.dropna()])
        )
        full_mi_sf = full_mi_sf.merge(file_summary, on="model", how="left")
        full_mi_sf["style"] = "MI_SF"
        full_mi_sf["run"] = "full"
        full_mi_sf["boot"] = 0
        full_mi_sf_raw.to_csv(settings.out_dir / "full_MI_SF_raw_mofa_cmf.csv", index=False)
        full_mi_sf_raw.to_pickle(settings.out_dir / "full_MI_SF_raw_mofa_cmf.pkl")
        full_mi_sf.to_csv(settings.out_dir / "full_MI_SF_mofa_cmf.csv", index=False)
        full_mi_sf.to_pickle(settings.out_dir / "full_MI_SF_mofa_cmf.pkl")
        full_outputs.append(full_mi_sf)

    # Combine and save all strategies
    full_all = pd.concat(full_outputs, ignore_index=True)
    full_all.to_csv(settings.out_dir / "full_ALL_styles_mofa_cmf.csv", index=False)
    full_all.to_pickle(settings.out_dir / "full_ALL_styles_mofa_cmf.pkl")
    print(full_all)
    return full_all
