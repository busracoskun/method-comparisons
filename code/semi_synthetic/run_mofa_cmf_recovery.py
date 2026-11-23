#!/usr/bin/env python3
"""Run the semi-synthetic modality-loss experiment for MOFA and CMF.

The data preparation, missingness masks, missForest imputations, and recovery metrics come
from code/semi_synthetic/modality_loss_and_recovery.R (called through rpy2), with the same
scenario seeds, so MOFA and CMF are evaluated on exactly the scenarios used for the R methods.
"""

from __future__ import annotations

import importlib.util
import os
import pickle
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector, StrVector

# SETTINGS

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", os.environ["CLEAN_ROOT"]))
BENCH_DIR = Path(os.environ["BENCH_WORK_DIR"])
SEMI_R_FILE = Path(__file__).resolve().parent / "modality_loss_and_recovery.R"
PYTHON_METHODS_FILE = Path(__file__).resolve().parents[2] / "code/methods" / "mofa_cmf.py"
OUT_DIR = Path(
    os.environ.get(
        "PY_SEMISYNTH_OUT_DIR", str(BENCH_DIR / "results_semi_synthetic_python_mofa_cmf_v1")
    )
)

# Run mode and the fixed design settings
MODE = os.environ.get("PY_SEMISYNTH_MODE", "run").lower()
TASK_ID = int(os.environ.get("PY_SEMISYNTH_TASK_ID", "0"))
N_FACTORS = 3
BETA = 2.0
MISSFOREST_NTREE = 100
MISSFOREST_MAXITER = 10
MODELS = tuple(os.environ.get("PY_METHODS", "mofa,cmf").split(","))

# The R experiment runs 12 methods; MOFA and CMF are methods 13 and 14, so their fitting seeds
# (scenario seed + index) are distinct from the R methods while the scenario seed is shared
MODEL_INDEX = {"mofa": 13, "cmf": 14}

# Output locations
GROUND_TRUTH_DIR = OUT_DIR / "ground_truth"
GROUND_TRUTH_FILE = GROUND_TRUTH_DIR / "ground_truth_mofa_cmf.pkl"
GROUND_TRUTH_STATUS_FILE = GROUND_TRUTH_DIR / "ground_truth_status_mofa_cmf.csv"
TASK_DIR = OUT_DIR / "tasks"
COMBINED_DIR = OUT_DIR / "combined"

for directory in (OUT_DIR, GROUND_TRUTH_DIR, TASK_DIR, COMBINED_DIR):
    directory.mkdir(parents=True, exist_ok=True)

# Load the shared MOFA/CMF module; @dataclass needs the module registered in sys.modules
if not PYTHON_METHODS_FILE.exists():
    raise FileNotFoundError(f"Missing mofa_cmf.py: {PYTHON_METHODS_FILE}")
spec = importlib.util.spec_from_file_location("python_methods_existing", PYTHON_METHODS_FILE)
pm = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = pm
spec.loader.exec_module(pm)


# R HELPERS


def load_r_helpers():
    """Source modality_loss_and_recovery.R up to its final marker line into an R environment.

    This provides prepare_dataset(), make_mar_score(), make_mnar_score(), make_missing(), and
    latent_metrics() exactly as used by the R experiment.
    """
    if not SEMI_R_FILE.exists():
        raise FileNotFoundError(f"Missing modality_loss_and_recovery.R: {SEMI_R_FILE}")

    # The part of the R script above its marker line
    text = SEMI_R_FILE.read_text()
    marker = "# ENTRY POINT"
    marker_position = text.find(marker)
    if marker_position < 0:
        raise RuntimeError("Could not find '# ENTRY POINT' in modality_loss_and_recovery.R")
    prefix = text[:marker_position]

    os.environ["PROJECT_ROOT"] = str(PROJECT_ROOT)
    os.environ["SEMISYNTH_OUT_DIR"] = str(OUT_DIR)
    # Read by the R script when it is sourced; they mirror the fixed design settings
    os.environ["SEMISYNTH_METHODS"] = "mofa,cmf"
    os.environ["SEMISYNTH_REPETITIONS"] = "100"
    os.environ["SEMISYNTH_N_FACTORS"] = "3"
    os.environ["SEMISYNTH_BETA"] = "2"
    os.environ["MISSFOREST_NTREE"] = "100"
    os.environ["MISSFOREST_MAXITER"] = "10"
    os.environ["SAVE_SEMISYNTH_DATA"] = "0"
    os.environ["RETRY_FAILED"] = "0"

    # Source it into a fresh R environment from the project root
    ro.r["source"](str(Path(__file__).resolve().parents[2] / "code" / "load_dataset.R"))
    helper_env = ro.r["new.env"](parent=ro.globalenv)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as handle:
        handle.write(prefix)
        temp_r = handle.name
    old_wd = os.getcwd()
    try:
        os.chdir(PROJECT_ROOT)
        ro.r["sys.source"](temp_r, envir=helper_env)
    finally:
        os.chdir(old_wd)
        try:
            os.unlink(temp_r)
        except OSError:
            pass
    # The same missForest call as run_task() in the R script
    extra_r = r"""
    python_semisynth_missforest <- function(X, seed, ntree = 100L, maxiter = 10L) {
      started <- proc.time()[["elapsed"]]
      tryCatch(
        {
          set.seed(as.integer(seed))
          fitted <- missForest::missForest(
            xmis = X, maxiter = as.integer(maxiter), ntree = as.integer(ntree), verbose = FALSE
          )
          list(
            status = "success", error_message = NA_character_,
            elapsed_seconds = proc.time()[["elapsed"]] - started,
            X = as.matrix(fitted$ximp), OOBerror = fitted$OOBerror
          )
        },
        error = function(e) {
          list(
            status = "failed", error_message = conditionMessage(e),
            elapsed_seconds = proc.time()[["elapsed"]] - started, X = NULL, OOBerror = NA
          )
        }
      )
    }
    """

    parsed = ro.r["parse"](text=extra_r)
    ro.r["eval"](parsed, envir=helper_env)
    return helper_env


R_HELPERS = load_r_helpers()


# CONVERSION HELPERS


def r_matrix_to_df(x) -> pd.DataFrame:
    arr = np.asarray(x, dtype=float)
    rownames = [str(v) for v in ro.r["rownames"](x)]
    colnames = [str(v) for v in ro.r["colnames"](x)]
    return pd.DataFrame(arr, index=pd.Index(rownames), columns=pd.Index(colnames))


def r_named_character_to_series(x) -> pd.Series:
    names = [str(v) for v in x.names]
    values = [str(v) for v in x]
    return pd.Series(values, index=pd.Index(names), dtype="object")


def df_to_r_matrix(df: pd.DataFrame):
    arr = df.to_numpy(dtype=float)
    r_matrix = ro.r["matrix"](
        FloatVector(arr.ravel(order="C")), nrow=arr.shape[0], ncol=arr.shape[1], byrow=True
    )
    r_matrix = ro.r["rownames<-"](r_matrix, StrVector(df.index.astype(str).tolist()))
    r_matrix = ro.r["colnames<-"](r_matrix, StrVector(df.columns.astype(str).tolist()))
    return r_matrix


# DATA AND TASK SETUP


def prepare_dataset(dataset: str) -> Dict[str, Any]:
    """Complete-case truth matrix from the R prepare_dataset(), in R and pandas form."""

    setup_r = R_HELPERS["prepare_dataset"](dataset)
    X_truth_r = setup_r.rx2("X_truth")
    layer_r = setup_r.rx2("layer_vec")
    X_truth = r_matrix_to_df(X_truth_r)
    layer_vec = r_named_character_to_series(layer_r)
    modalities = [str(v) for v in setup_r.rx2("modalities")]
    dataset_label = str(setup_r.rx2("dataset_label")[0])
    return {
        "dataset": dataset, "dataset_label": dataset_label, "setup_r": setup_r,
        "X_truth_r": X_truth_r, "layer_r": layer_r, "X_truth": X_truth, "layer_vec": layer_vec,
        "modalities": modalities,
    }


# Task mapping as in the R script: tasks 1-100 MCAR 20%, 101-200 MAR 20%, 201-300 MNAR 20%,
# 301-400 MCAR 70%, 401-500 MAR 70%, 501-600 MNAR 70%
CONDITIONS = [
    ("MCAR", 0.20), ("MAR", 0.20), ("MNAR", 0.20), ("MCAR", 0.70), ("MAR", 0.70), ("MNAR", 0.70)
]


def map_task(task_id: int) -> Dict[str, Any]:
    if task_id < 1 or task_id > 600:
        raise ValueError("Task ID must be between 1 and 600.")
    condition_index = (task_id - 1) // 100
    repetition = ((task_id - 1) % 100) + 1
    mechanism, rate = CONDITIONS[condition_index]
    return {
        "task_id": task_id, "condition_id": condition_index + 1, "repetition": repetition,
        "mechanism": mechanism, "target_rate": rate,
    }


# FITTING AND METRICS


def fit_model(model: str, X: pd.DataFrame, layer_vec: pd.Series, seed: int) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        if model == "mofa":
            Z, W, _ = pm.fit_mofa(
                X=X, layer_vec=layer_vec.reindex(X.columns), n_factors=N_FACTORS, seed=seed
            )
        elif model == "cmf":
            Z, W, _ = pm.fit_cmf(
                X=X, layer_vec=layer_vec.reindex(X.columns), n_factors=N_FACTORS, seed=seed,
                lambda_=0.1, niter=30, method="als", use_item_view_info=True
            )
        else:
            raise ValueError(f"Unsupported model: {model}")

        # Check the fitted scores and loadings
        Z = pd.DataFrame(Z).reindex(X.index)
        if Z.shape[1] < 1:
            raise RuntimeError("Latent score matrix has zero factors.")
        if not np.isfinite(Z.to_numpy(dtype=float)).all():
            raise RuntimeError("Latent score matrix contains NA/NaN/Inf.")
        if W is not None:
            W = pd.DataFrame(W).reindex(X.columns)
            if not np.isfinite(W.to_numpy(dtype=float)).all():
                raise RuntimeError("Loading matrix contains NA/NaN/Inf.")
        return {
            "status": "success", "error_message": "",
            "elapsed_seconds": time.perf_counter() - started, "Z": Z, "W": W,
        }
    except Exception as exc:
        return {
            "status": "failed", "error_message": str(exc),
            "elapsed_seconds": time.perf_counter() - started, "Z": None, "W": None,
        }


def compare_fit(baseline: Dict[str, Any], fitted: Dict[str, Any], layer_r) -> Dict[str, Any]:
    """Recovery metrics from the R latent_metrics(), so they match the R experiment exactly."""

    # A comparison needs both the reference fit and the new fit
    if baseline["status"] != "success":
        return {
            "status": "baseline_failed",
            "error_message": baseline.get("error_message", "Ground-truth fit failed."),
        }
    if fitted["status"] != "success":
        return {
            "status": "fit_failed",
            "error_message": fitted.get("error_message", "Corrupted fit failed."),
        }
    try:
        # Convert both fits to R matrices and compute the metrics
        Z0 = df_to_r_matrix(baseline["Z"])
        Z1 = df_to_r_matrix(fitted["Z"])
        W0 = df_to_r_matrix(baseline["W"]) if baseline["W"] is not None else ro.NULL
        W1 = df_to_r_matrix(fitted["W"]) if fitted["W"] is not None else ro.NULL
        metrics_r = R_HELPERS["latent_metrics"](
            Z0=Z0, Z1=Z1, W0=W0, W1=W1, block=layer_r, gw=True, top_k=50
        )

        # Named R vector to a plain dict
        metric_names = [str(v) for v in metrics_r.names]
        metric_values = np.asarray(metrics_r, dtype=float)
        out = {"status": "success", "error_message": ""}
        for name, value in zip(metric_names, metric_values):
            out[name] = float(value)
        return out
    except Exception as exc:
        return {"status": "metric_failed", "error_message": str(exc)}


# GROUND TRUTH


def build_ground_truth():
    """Fit MOFA and CMF on the complete-case X_truth of this experiment (seed 10000 + index).

    These reference fits are specific to the semi-synthetic design; they are not the
    observed-data full-cohort fits.
    """
    print("============================================================")
    print("BUILDING PYTHON MOFA + CMF SEMI-SYNTHETIC GROUND TRUTH")
    print("============================================================")
    ground_truth = {
        "methods": list(MODELS), "n_factors": N_FACTORS, "model_index": MODEL_INDEX, "datasets": {}
    }

    status_rows = []
    for dataset in ("adni", "tcga"):
        if dataset != os.environ["BENCH_DATASET"]:
            continue

        setup = prepare_dataset(dataset)
        print()
        print(setup["dataset_label"])
        print("Samples:", setup["X_truth"].shape[0])
        print("Features:", setup["X_truth"].shape[1])
        print("Modalities:", ", ".join(setup["modalities"]))

        # Fit each model on the complete data with its reference seed
        fits = {}
        for model in MODELS:
            seed = 10000 + MODEL_INDEX[model]
            print(
                f"Ground truth: {setup['dataset_label']} {model} seed={seed}",
                flush=True,
            )
            fit = fit_model(
                model=model, X=setup["X_truth"], layer_vec=setup["layer_vec"], seed=seed
            )
            fits[model] = fit
            status_rows.append({
                "dataset": setup["dataset_label"], "model": model, "seed": seed,
                "status": fit["status"], "error_message": fit["error_message"],
                "elapsed_seconds": fit["elapsed_seconds"], "n_samples": setup["X_truth"].shape[0],
                "n_features": setup["X_truth"].shape[1],
                "modalities": ";".join(setup["modalities"]),
            })

        # Keep the sample and feature order so later tasks can check them
        ground_truth["datasets"][dataset] = {
            "dataset_label": setup["dataset_label"],
            "sample_ids": list(setup["X_truth"].index.astype(str)),
            "feature_names": list(setup["X_truth"].columns.astype(str)),
            "layer_vec": setup["layer_vec"], "modalities": setup["modalities"], "fits": fits
        }

    # Save the reference fits and their status
    with open(GROUND_TRUTH_FILE, "wb") as handle:
        pickle.dump(ground_truth, handle)
    status_df = pd.DataFrame(status_rows)
    status_df.to_csv(GROUND_TRUTH_STATUS_FILE, index=False)
    print()
    print(status_df.to_string(index=False))
    print()
    print("Ground truth:", GROUND_TRUTH_FILE)
    print("Status:", GROUND_TRUTH_STATUS_FILE)
    if (status_df["status"] != "success").any():
        raise RuntimeError("One or more MOFA/CMF ground-truth fits failed.")


# SCENARIOS AND TASKS


def generate_scenario(
    setup: Dict[str, Any], task: Dict[str, Any], dataset_index: int, modality_index: int,
    modality: str
) -> Dict[str, Any]:
    """Mask one modality and impute with missForest using the R functions and seeds."""
    # Same scenario-seed formula as the R script; imputation uses scenario_seed + 500000
    scenario_seed = int(1_000_000 + task["task_id"] * 100 + dataset_index * 10 + modality_index)

    # Mechanism scores, then the mask
    block_cols = setup["layer_vec"].index[setup["layer_vec"] == modality].astype(str).tolist()
    if task["mechanism"] == "MAR":
        mar_score = R_HELPERS["make_mar_score"](setup["X_truth_r"], setup["layer_r"], modality)
    else:
        mar_score = ro.NULL
    if task["mechanism"] == "MNAR":
        mnar_score = R_HELPERS["make_mnar_score"](setup["X_truth_r"], setup["layer_r"], modality)
    else:
        mnar_score = ro.NULL
    simulated = R_HELPERS["make_missing"](
        X=setup["X_truth_r"], block_cols=StrVector(block_cols), mechanism=task["mechanism"],
        target_rate=float(task["target_rate"]), observed_score=mar_score, latent_score=mnar_score,
        beta=BETA, seed=scenario_seed
    )

    X_miss_r = simulated.rx2("X_miss")
    missing = np.asarray(simulated.rx2("missing"), dtype=bool)
    realized_rate = float(simulated.rx2("realized_rate")[0])

    # Impute with missForest in R
    imputation_seed = scenario_seed + 500_000
    imputation = R_HELPERS["python_semisynth_missforest"](
        X_miss_r, int(imputation_seed), int(MISSFOREST_NTREE), int(MISSFOREST_MAXITER)
    )
    imp_status = str(imputation.rx2("status")[0])
    imp_elapsed = float(imputation.rx2("elapsed_seconds")[0])
    if imp_status == "success":
        X_imputed = r_matrix_to_df(imputation.rx2("X"))
        imp_error = ""
    else:
        X_imputed = None
        error_r = imputation.rx2("error_message")
        imp_error = "" if len(error_r) == 0 else str(error_r[0])
    return {
        "scenario_seed": scenario_seed, "imputation_seed": imputation_seed,
        "realized_rate": realized_rate, "n_missing_samples": int(missing.sum()),
        "imputation_status": imp_status, "imputation_error": imp_error,
        "imputation_elapsed_seconds": imp_elapsed, "X_imputed": X_imputed,
    }


def run_task(task_id: int):
    """Run one of the 600 tasks, resuming from a checkpoint after each completed block."""

    # Output files for this task
    task = map_task(task_id)
    tag = f"{task_id:04d}"
    result_csv = TASK_DIR / f"task_{tag}_results.csv"
    result_pkl = TASK_DIR / f"task_{tag}_results.pkl"
    checkpoint_csv = TASK_DIR / f"task_{tag}_checkpoint.csv"
    if result_csv.exists():
        print(f"Task {task_id} already has a final result; skipping.")
        return

    # Load the reference fits
    if not GROUND_TRUTH_FILE.exists():
        raise FileNotFoundError(
            "Python MOFA/CMF semi-synthetic ground truth does not exist: " + str(GROUND_TRUTH_FILE)
        )
    with open(GROUND_TRUTH_FILE, "rb") as handle:
        ground_truth = pickle.load(handle)
    if ground_truth["n_factors"] != N_FACTORS:
        raise RuntimeError("Ground-truth factor count does not match.")

    print("============================================================")
    print(f"PYTHON MOFA/CMF SEMI-SYNTHETIC TASK {task_id} OF 600")
    print("Condition:", task["mechanism"], task["target_rate"])
    print("Repetition:", task["repetition"])
    print("Models:", ", ".join(MODELS))
    print("============================================================")

    # Resume from a checkpoint if one exists
    rows = []
    if checkpoint_csv.exists():
        checkpoint = pd.read_csv(checkpoint_csv)
        rows = checkpoint.to_dict(orient="records")
        print(f"Resuming from checkpoint with {len(rows)} rows.")

    # Loop over datasets and dropped modalities
    for dataset_index, dataset in enumerate(("adni", "tcga"), start=1):
        if dataset != os.environ["BENCH_DATASET"]:
            continue
        setup = prepare_dataset(dataset)
        baseline_dataset = ground_truth["datasets"][dataset]
        current_samples = list(setup["X_truth"].index.astype(str))
        current_features = list(setup["X_truth"].columns.astype(str))
        if current_samples != baseline_dataset["sample_ids"]:
            raise RuntimeError(
                f"{setup['dataset_label']} sample IDs changed after ground truth was created."
            )
        if current_features != baseline_dataset["feature_names"]:
            raise RuntimeError(
                f"{setup['dataset_label']} features changed after ground truth was created."
            )

        for modality_index, modality in enumerate(setup["modalities"], start=1):
            completed = {
                str(row["model"]).lower()
                for row in rows
                if (
                    row.get("dataset") == setup["dataset_label"]
                    and row.get("dropped_modality") == modality
                )
            }
            if all(model in completed for model in MODELS):
                print(
                    "Checkpoint: skipping completed block:", setup["dataset_label"], "/", modality
                )
                continue
            print()
            print(
                "Scenario:", setup["dataset_label"], "/", modality, "/", task["mechanism"], "/",
                task["target_rate"], flush=True
            )

            scenario = generate_scenario(
                setup=setup, task=task, dataset_index=dataset_index, modality_index=modality_index,
                modality=modality
            )

            # Fit each model on the imputed data and compare with its reference fit
            for model in MODELS:
                # Skip models already recorded for a partially checkpointed block
                if model in completed:
                    continue
                baseline = baseline_dataset["fits"][model]
                model_seed = scenario["scenario_seed"] + MODEL_INDEX[model]
                if scenario["imputation_status"] != "success":
                    fitted = {
                        "status": "imputation_failed",
                        "error_message": scenario["imputation_error"],
                        "elapsed_seconds": 0.0, "Z": None, "W": None,
                    }
                elif baseline["status"] != "success":
                    fitted = {
                        "status": "not_run",
                        "error_message": "Skipped because complete-data baseline failed.",
                        "elapsed_seconds": 0.0, "Z": None, "W": None,
                    }
                else:
                    print(f"Running {model} | seed={model_seed}", flush=True)
                    fitted = fit_model(
                        model=model, X=scenario["X_imputed"], layer_vec=setup["layer_vec"],
                        seed=model_seed
                    )

                comparison = compare_fit(baseline=baseline, fitted=fitted, layer_r=setup["layer_r"])
                row = {
                    "task_id": task_id, "condition_id": task["condition_id"],
                    "repetition": task["repetition"], "dataset": setup["dataset_label"],
                    "dropped_modality": modality, "mechanism": task["mechanism"],
                    "target_rate": task["target_rate"], "realized_rate": scenario["realized_rate"],
                    "n_samples": setup["X_truth"].shape[0],
                    "n_missing_samples": scenario["n_missing_samples"],
                    "n_features": setup["X_truth"].shape[1], "model": model,
                    "n_factors": N_FACTORS, "beta": BETA, "imputation": "missForest",
                    "missforest_ntree": MISSFOREST_NTREE, "missforest_maxiter": MISSFOREST_MAXITER,
                    "imputation_status": scenario["imputation_status"],
                    "imputation_error": scenario["imputation_error"],
                    "imputation_elapsed_seconds": scenario["imputation_elapsed_seconds"],
                    "fit_status": fitted["status"], "fit_error": fitted["error_message"],
                    "fit_elapsed_seconds": fitted["elapsed_seconds"],
                    "status": comparison["status"],
                    "error_message": comparison.get("error_message", ""),
                    "n_attempted": 1, "n_successful": int(comparison["status"] == "success"),
                    "n_failed": int(comparison["status"] != "success"),
                    # As in the R results, "seed" is the missingness scenario seed
                    "seed": scenario["scenario_seed"],
                }

                for key, value in comparison.items():
                    if key not in ("status", "error_message"):
                        row[key] = value
                rows.append(row)

            pd.DataFrame(rows).to_csv(checkpoint_csv, index=False)

    # Save the task results and failures
    results = pd.DataFrame(rows)
    results = results.sort_values(["dataset", "dropped_modality", "model"]).reset_index(drop=True)
    results.to_csv(result_csv, index=False)
    results.to_pickle(result_pkl)
    if checkpoint_csv.exists():
        checkpoint_csv.unlink()

    failures = results[results["status"] != "success"]
    failures.to_csv(TASK_DIR / f"task_{tag}_failures.csv", index=False)
    print()
    print("============================================================")
    print(f"TASK {task_id} COMPLETE")
    print("Rows:", len(results))
    print("Successes:", int((results["status"] == "success").sum()))
    print("Failures:", int((results["status"] != "success").sum()))
    print("Saved:", result_csv)
    print("============================================================")


def combine_tasks():
    # Read and combine all completed tasks
    files = sorted(TASK_DIR.glob("task_*_results.csv"))
    if not files:
        raise RuntimeError("No completed Python semi-synthetic task files were found.")
    frames = [pd.read_csv(f) for f in files]
    results = pd.concat(frames, ignore_index=True, sort=False)
    results = results.sort_values(["task_id", "dataset", "dropped_modality", "model"]).reset_index(
        drop=True
    )
    combined_csv = COMBINED_DIR / "a_semi_synthetic_mofa_cmf_ALL.csv"
    combined_pkl = COMBINED_DIR / "a_semi_synthetic_mofa_cmf_ALL.pkl"
    results.to_csv(combined_csv, index=False)
    results.to_pickle(combined_pkl)

    # Fit counts per condition, failures, and task coverage
    counts = (
        results.groupby(
            ["dataset", "dropped_modality", "mechanism", "target_rate", "model"], dropna=False
        )
        .agg(
            n_attempted=("status", "size"),
            n_successful=("status", lambda x: int((x == "success").sum())),
        )
        .reset_index()
    )
    counts["n_failed"] = counts["n_attempted"] - counts["n_successful"]
    counts["success_rate"] = counts["n_successful"] / counts["n_attempted"]
    counts.to_csv(COMBINED_DIR / "a_semi_synthetic_mofa_cmf_run_counts.csv", index=False)
    failures = results[results["status"] != "success"]
    failures.to_csv(COMBINED_DIR / "a_semi_synthetic_mofa_cmf_failures.csv", index=False)
    completed_tasks = set(
        pd.to_numeric(results["task_id"], errors="coerce").dropna().astype(int).unique()
    )
    coverage = pd.DataFrame({"task_id": np.arange(1, 601, dtype=int)})
    coverage["completed"] = coverage["task_id"].isin(completed_tasks)
    coverage.to_csv(COMBINED_DIR / "a_semi_synthetic_mofa_cmf_task_coverage.csv", index=False)

    print("============================================================")
    print("COMBINED PYTHON MOFA/CMF SEMI-SYNTHETIC RESULTS")
    print("Completed tasks:", len(completed_tasks), "of 600")
    print("Rows:", len(results))
    print("Combined CSV:", combined_csv)
    print("============================================================")


if __name__ == "__main__":
    if MODE == "ground_truth":
        build_ground_truth()
    elif MODE == "run":
        if TASK_ID < 1:
            raise RuntimeError("Set PY_SEMISYNTH_TASK_ID=1..600")
        run_task(TASK_ID)
    elif MODE == "combine":
        combine_tasks()
    else:
        raise RuntimeError("PY_SEMISYNTH_MODE must be ground_truth, run, or combine.")
