"""Prepare TCGA inputs and evaluate the Python integration methods."""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

import python_methods as base


# PATHS

BENCH_DIR = Path(os.environ["BENCH_WORK_DIR"])

TCGA_INPUT_DIR = Path(os.environ.get("TCGA_INPUT_DIR", os.environ["TCGA_INPUT_DIR"]))

TCGA_MRNA_FILE = TCGA_INPUT_DIR / "irsem_fs.csv"
TCGA_METHYLATION_FILE = TCGA_INPUT_DIR / "imeth_fs.csv"
TCGA_RPPA_FILE = TCGA_INPUT_DIR / "irppa_fs.csv"
TCGA_MIRNA_FILE = TCGA_INPUT_DIR / "imirna_fs.csv"
TCGA_CLINICAL_FILE = TCGA_INPUT_DIR / "clin_common.csv"

VARIABLE_DIR = BENCH_DIR / "variables"
VARIABLE_DIR.mkdir(parents=True, exist_ok=True)


# TCGA INPUT HELPERS

NA_VALUES = ["", "NA", "N/A", "NaN", "NULL"]


def read_prepared_omics(path: Path, modality: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required {modality} file was not found: {path}")

    dat = pd.read_csv(path, na_values=NA_VALUES, keep_default_na=True)

    if dat.shape[1] < 2:
        raise ValueError(f"{modality} must contain an ID column plus features.")

    raw_ids = dat.iloc[:, 0]

    if raw_ids.isna().any():
        raise ValueError(f"{modality} contains missing patient IDs.")

    ids = raw_ids.astype(str)

    if (ids.str.len() == 0).any():
        raise ValueError(f"{modality} contains blank patient IDs.")

    if ids.duplicated().any():
        raise ValueError(f"{modality} contains duplicated patient IDs.")

    X = dat.iloc[:, 1:].copy()

    for col in X.columns:
        X[col] = pd.to_numeric(X[col], errors="raise")

    X.index = ids.to_numpy()
    X.index.name = "patient_id"

    return X.astype(float)


def read_prepared_clinical(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required TCGA clinical file was not found: {path}")

    clin = pd.read_csv(path, na_values=NA_VALUES, keep_default_na=True)

    if clin.shape[1] < 2:
        raise ValueError("clin_common.csv must contain an ID column and clinical variables.")

    raw_ids = clin.iloc[:, 0]

    if raw_ids.isna().any():
        raise ValueError("clin_common.csv contains missing patient IDs.")

    ids = raw_ids.astype(str)

    if (ids.str.len() == 0).any():
        raise ValueError("clin_common.csv contains blank patient IDs.")

    if ids.duplicated().any():
        raise ValueError("clin_common.csv contains duplicated patient IDs.")

    clin.index = ids.to_numpy()
    clin.index.name = "patient_id"

    required = [
        "disease_code",
        "karnofsky_performance_score",
        "missingness",
    ]

    missing = [x for x in required if x not in clin.columns]

    if missing:
        raise ValueError("clin_common.csv is missing required columns: " + ", ".join(missing))

    return clin


def ordered_union(indexes) -> pd.Index:
    out = []
    seen = set()

    for idx in indexes:
        for value in idx:
            value = str(value)
            if value not in seen:
                seen.add(value)
                out.append(value)

    return pd.Index(out)


def r_make_unique(names: List[str]) -> List[str]:
    """
    Match R's make.unique(names, sep=".") so duplicate feature names
    are handled exactly like the TCGA R pipeline.
    """
    import rpy2.robjects as ro
    from rpy2.robjects.vectors import StrVector

    result = ro.r["make.unique"](StrVector(names), ".")
    return [str(x) for x in result]


# TCGA DATA LOADER

def load_and_process_data_tcga(settings: base.Settings) -> Dict[str, object]:

    print("=" * 70)
    print("LOADING PREPARED TCGA DATA")
    print("=" * 70)
    print("Input directory:", TCGA_INPUT_DIR)

    omics = {
        "mrna": read_prepared_omics(TCGA_MRNA_FILE, "mRNA"),
        "methylation": read_prepared_omics(TCGA_METHYLATION_FILE, "methylation"),
        "rppa": read_prepared_omics(TCGA_RPPA_FILE, "RPPA"),
    }

    # Match final_tcga.R: include prepared miRNA only if it exists.
    if TCGA_MIRNA_FILE.exists():
        print("Including optional prepared miRNA:", TCGA_MIRNA_FILE)
        omics["mirna"] = read_prepared_omics(TCGA_MIRNA_FILE, "miRNA")
    else:
        print("No optional miRNA file; using mRNA + methylation + RPPA.")

    clin = read_prepared_clinical(TCGA_CLINICAL_FILE)

    ids_by_modality = {
        name: df.index
        for name, df in omics.items()
    }

    all_omics_ids = ordered_union([df.index for df in omics.values()])

    # Keep patients with matching clinical data.
    ids_score = pd.Index([x for x in all_omics_ids if x in clin.index])

    diagnosis_raw = clin.loc[
        ids_score,
        "disease_code",
    ]

    keep_diag = (diagnosis_raw.notna() & diagnosis_raw.astype(str).str.len().gt(0))

    ids_score = ids_score[keep_diag.to_numpy()]
    diagnosis_raw = diagnosis_raw.loc[ids_score].astype(str)

    if len(ids_score) < 2:
        raise ValueError("Fewer than two TCGA patients remain after clinical matching.")

    # Align every modality to the same union sample set.
    aligned = {
        name: df.reindex(ids_score)
        for name, df in omics.items()
    }

    # Match R make.unique() for duplicate biological names across modalities.
    original_names = []
    feature_modalities = []

    for name, df in aligned.items():
        original_names.extend([str(x) for x in df.columns])
        feature_modalities.extend([name] * df.shape[1])

    internal_names = r_make_unique(original_names)

    cursor = 0

    for name in aligned:
        width = aligned[name].shape[1]
        aligned[name].columns = internal_names[
            cursor:cursor + width
        ]
        cursor += width

    X_all = pd.concat(list(aligned.values()), axis=1).astype(float)

    # Match final_tcga.R:
    # standardize once using observed means and sample SDs.
    feature_means = X_all.mean(axis=0, skipna=True)

    feature_sds = X_all.std(axis=0, skipna=True, ddof=1)

    bad_sd = (~np.isfinite(feature_sds.to_numpy()) | (feature_sds.to_numpy() <= 0))

    if np.any(bad_sd):
        bad_names = list(X_all.columns[np.where(bad_sd)[0]])

        raise ValueError(
            "TCGA contains zero-variance or invalid features "
            "before standardization. First examples: "
            + ", ".join(bad_names[:10])
        )

    X_stand = ((X_all - feature_means) / feature_sds)

    layer_vec = pd.Series(feature_modalities, index=internal_names, name="layer")

    diag_labels = pd.Series(pd.Categorical(diagnosis_raw), index=ids_score, name="disease_code")

    # Require exactly two disease groups, matching final_tcga.R.
    disease_levels = sorted(diagnosis_raw.unique().tolist())

    if len(disease_levels) != 2:
        raise ValueError(
            "TCGA SIS requires exactly two disease_code groups. Found: "
            + ", ".join(disease_levels)
        )

    missingness_raw = clin.loc[
        ids_score,
        "missingness",
    ]

    if (missingness_raw.isna().any() or missingness_raw.astype(str).str.len().eq(0).any()):
        raise ValueError("TCGA clinical missingness contains missing/blank values.")

    missingness_raw = missingness_raw.astype(str)

    missingness_levels = (["CC"] + sorted(x for x in missingness_raw.unique() if x != "CC"))

    miss_labels = pd.Series(
        pd.Categorical(missingness_raw, categories=missingness_levels),
        index=ids_score,
        name="missingness",
    )

    karnofsky = pd.to_numeric(
        clin.loc[
            ids_score,
            "karnofsky_performance_score",
        ],
        errors="coerce",
    )

    mmse_vec = pd.Series(karnofsky.to_numpy(), index=ids_score, name="karnofsky_performance_score")

    # Patients present in every included omics modality.
    common_set = set(ids_score)

    for idx in ids_by_modality.values():
        common_set &= set(idx.astype(str))

    common_ids = pd.Index([x for x in ids_score if x in common_set])

    if len(common_ids) < 2:
        raise ValueError("Fewer than two TCGA patients are present in every modality.")

    if X_stand.index.duplicated().any():
        raise ValueError("Duplicate TCGA sample IDs remain.")

    if X_stand.columns.duplicated().any():
        raise ValueError("Duplicate TCGA internal feature names remain.")

    print()
    print("=" * 70)
    print("TCGA PYTHON INPUT SUMMARY")
    print("=" * 70)
    print("Modalities:", ", ".join(omics.keys()))
    print("Patients used:", X_stand.shape[0])
    print("Features used:", X_stand.shape[1])
    print("Complete-modality patients:", len(common_ids))
    print("Disease counts:")
    print(diag_labels.value_counts(dropna=False))
    print("Missingness counts:")
    print(miss_labels.value_counts(dropna=False))
    print("=" * 70)

    return {
        "X_stand": X_stand,
        "layer_vec": layer_vec,
        "diag_labels": diag_labels,
        "miss_labels": miss_labels,
        "mmse_vec": mmse_vec,
        "common_ids": common_ids,
    }


# TCGA SIS

def sel_SIS_tcga(X_complete: pd.DataFrame, label_vec: pd.Series, seed: int = 100) -> List[str]:

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.vectors import IntVector
    from rpy2.robjects.packages import importr

    y_raw = label_vec.reindex(X_complete.index)

    keep = (y_raw.notna() & ~X_complete.isna().any(axis=1))

    X_use = X_complete.loc[keep].copy()
    labels = y_raw.loc[keep].astype(str)

    levels = sorted(labels.unique().tolist())

    if len(levels) != 2:
        raise ValueError(
            "TCGA SIS requires exactly two disease groups. Found: "
            + ", ".join(levels)
        )

    if min(labels.value_counts()) < 2:
        raise ValueError("TCGA SIS requires at least two observations " "in each disease group.")

    # Match final_tcga.R:
    # second factor level = 1, first = 0.
    ybin = (labels == levels[1]).astype(int).to_numpy()

    SIS_pkg = importr("SIS")

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_X_df = ro.conversion.py2rpy(X_use)

    r_y = IntVector(ybin.astype(int).tolist())

    ro.r["set.seed"](int(seed))

    fit = SIS_pkg.SIS(
        x=ro.r["as.matrix"](r_X_df),
        y=r_y,
        family="binomial",
        tune="bic",
        nsis=int(X_use.shape[1]),
    )

    selected_r = np.asarray(fit.rx2("sis.ix0"), dtype=int)

    if selected_r.size == 0:
        raise ValueError("TCGA SIS selected zero features.")

    return X_use.columns[
        selected_r - 1
    ].tolist()


# TCGA-SPECIFIC MULTIPLE IMPUTATION CACHE

def load_or_create_tcga_imputations(X_na: pd.DataFrame, settings: base.Settings):

    import rpy2.robjects as ro
    from rpy2.robjects import (pandas2ri, numpy2ri)
    from rpy2.robjects.conversion import localconverter

    M = int(settings.m_imputations)
    seed = int(settings.imputation_seed)

    imp_file = VARIABLE_DIR / (f"tcga_imp_list_M{M}_seed{seed}.rds")

    local_imp_file = imp_file
    external = Path(os.environ.get("CACHE_INPUT_DIR") or VARIABLE_DIR) / imp_file.name
    if external.exists():
        imp_file = external

    if imp_file.exists():
        print("Loading existing TCGA imputation cache:", imp_file, flush=True)

        r_list = ro.r["readRDS"](str(imp_file))

        if len(r_list) == M:
            converted = []
            valid = True

            for i in range(M):
                with localconverter(ro.default_converter + numpy2ri.converter):
                    arr = ro.conversion.rpy2py(r_list[i])

                arr = np.asarray(arr, dtype=float)

                if (arr.shape != X_na.shape or not np.isfinite(arr).all()):
                    valid = False
                    break

                converted.append(pd.DataFrame(arr, index=X_na.index, columns=X_na.columns))

            if valid:
                print("Existing TCGA imputations passed validation.", flush=True)
                return converted

        print("Existing TCGA imputation cache is incompatible; " "recomputing it.", flush=True)

    print("Computing 5 TCGA missForest imputations...", flush=True)

    imp_file = local_imp_file
    imputations = base.create_imputations(X_na, M=M, seed=seed)

    r_objects = {}

    for i, Xi in enumerate(imputations, start=1):
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_df = ro.conversion.py2rpy(Xi)

        r_objects[str(i)] = ro.r[
            "as.matrix"
        ](r_df)

    r_list = ro.ListVector(r_objects)

    ro.r["saveRDS"](r_list, str(imp_file))

    print("Saved TCGA imputations to:", imp_file, flush=True)

    return imputations


# PATCH ONLY THE DATASET-SPECIFIC PIECES

base.load_and_process_data = load_and_process_data_tcga
base.sel_SIS = sel_SIS_tcga
base.load_or_create_imputations = load_or_create_tcga_imputations


# MAIN

