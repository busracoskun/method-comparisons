#!/usr/bin/env python3
"""Recalculate the bootstrap and recovery summary tables from the saved run-level records.

With --verify, every recalculated table is checked against its archived copy in results/.
Outputs are written to build/recalculated by default.
"""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "config/default.json").is_file())

# Metrics summarized for the observed-data and semi-synthetic analyses
OBS = [
    "rand_diag", "ari_diag", "chi_diag", "db_diag", "pve_diag", "R2",
    "rand_miss", "ari_miss", "chi_miss", "db_miss", "pve_miss",
]
SEMI = [
    "rv", "gw", "recon_relerr_full", "recon_cor", "rv_W", "wgene_rankcor", "wgene_jaccard",
    "recon_relerr_gene", "recon_relerr_meta", "recon_relerr_mrna", "recon_relerr_methylation",
    "recon_relerr_rppa",
]


def read(path):
    return pd.read_csv(
        path, keep_default_na=False, na_values=["", "NaN", "nan", "NULL"], low_memory=False,
        float_precision="round_trip"
    )


def summarize(d, keys, metrics, observed=False):
    """Mean, median, and 2.5%/97.5% percentiles of the finite values of each metric.

    Observed-data records count when at least one fit succeeded (status success or partial);
    semi-synthetic records only when the status is success. Missingness-clustering metrics
    are undefined for CC_SF (every complete case has class "CC") and are left empty.
    """

    rows = []

    for key, g in d.groupby(keys, dropna=False, sort=True):
        # Records that count as successful
        base = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        good = (
            g.loc[g.status.isin(["success", "partial"]) & g.n_successful.gt(0)]
            if observed
            else g.loc[g.status.eq("success")]
        )

        # Finite values of each metric, then their summary statistics
        for metric in metrics:
            if metric not in g:
                continue
            x = pd.to_numeric(good[metric], errors="coerce").to_numpy(float)
            x = x[np.isfinite(x)]
            if base.get("style") == "CC_SF" and metric.endswith("_miss"):
                x = np.array([])
            rows.append(
                base
                | dict(
                    metric=metric, recorded_records=len(g), successful_records=len(good), n=len(x),
                    mean=np.mean(x) if len(x) else np.nan,
                    median=np.median(x) if len(x) else np.nan,
                    low=np.percentile(x, 2.5) if len(x) else np.nan,
                    high=np.percentile(x, 97.5) if len(x) else np.nan,
                )
            )
    return pd.DataFrame(rows)


def compare(a, b, keys, columns):
    """Assert that two summary tables have identical keys and numerically equal columns."""

    # Same keys in the same order, then numerically equal columns
    a = a.sort_values(keys).reset_index(drop=True)
    b = b.sort_values(keys).reset_index(drop=True)
    assert len(a) == len(b) and a[keys].equals(b[keys]), "Mismatched summary keys"
    for col in columns:
        np.testing.assert_allclose(
            a[col], b[col], rtol=1e-10, atol=1e-12, equal_nan=True, err_msg=col
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=ROOT / "build/recalculated")
    p.add_argument("--verify", action="store_true")
    a = p.parse_args()

    # Output folder (never inside results/)
    out = a.output.resolve()
    if out == ROOT / "results" or (ROOT / "results") in out.parents:
        p.error("Choose an output outside frozen results/")
    out.mkdir(parents=True, exist_ok=True)

    # Metric summaries
    outputs = {}
    summaries = [
        (
            "observed_data", "bootstrap_metrics_by_replicate", ["dataset", "model", "style"],
            OBS, True,
        ),
        (
            "semi_synthetic", "recovery_metrics_by_repetition",
            ["dataset", "model", "dropped_modality", "mechanism", "target_rate"], SEMI, False,
        ),
        (
            "semi_synthetic", "intnmf_sensitivity_metrics_by_repetition",
            ["dataset", "experiment", "model", "dropped_modality", "mechanism", "target_rate"],
            SEMI, False,
        ),
    ]

    # Summarize each run-level table and optionally compare with the archived summary
    for folder, filename, keys, metrics, obs in summaries:
        d = read(ROOT / f"results/{folder}/{filename}.csv")
        s = summarize(d, keys, metrics, obs)
        if obs:
            name = "bootstrap_metric_summary"
        elif filename == "intnmf_sensitivity_metrics_by_repetition":
            name = "intnmf_sensitivity_metric_summary"
        else:
            name = "recovery_metric_summary"

        if a.verify:
            compare(
                s, read(ROOT / f"results/{folder}/{name}.csv"), keys + ["metric"],
                ["recorded_records", "successful_records", "n", "mean", "median", "low", "high"]
            )
        output_file = out / f"{folder.replace('/', '_')}_{name}.csv"
        s.to_csv(output_file, index=False, float_format="%.17g")
        outputs[name] = len(s)

    # Observed-data fit coverage, absent replicates, and recorded run times
    d = read(ROOT / "results/observed_data/bootstrap_metrics_by_replicate.csv")

    # Fit counts per dataset, method, and strategy, plus the missing replicate numbers
    counts = []
    absent = []
    for (ds, model, style), g in d.groupby(["dataset", "model", "style"]):
        counts.append(
            dict(
                dataset=ds, model=model, style=style, records=len(g), expected_replicates=1000,
                absent_replicates=1000 - len(g),
                successful_replicates=int(
                    (g.status.isin(["success", "partial"]) & g.n_successful.gt(0)).sum()
                ),
                fully_successful_replicates=int(g.status.eq("success").sum()),
                partial_replicates=int(g.status.eq("partial").sum()),
                failed_replicates=int(g.n_successful.eq(0).sum()),
                attempted_fits=int(g.n_attempted.sum()),
                successful_fits=int(g.n_successful.sum()),
                failed_fits=int(g.n_failed.sum()),
            )
        )
        for boot in sorted(set(range(1, 1001)) - set(g.boot)):
            absent.append(dict(dataset=ds, model=model, style=style, boot=boot))

    c = pd.DataFrame(counts)
    if a.verify:
        compare(
            c, read(ROOT / "results/observed_data/bootstrap_fit_success_counts.csv"),
            ["dataset", "model", "style"], [x for x in c if x not in ["dataset", "model", "style"]]
        )
    c.to_csv(out / "observed_fit_counts.csv", index=False)

    missing = pd.DataFrame(absent)
    if a.verify:
        compare(
            missing, read(ROOT / "results/observed_data/missing_bootstrap_records.csv"),
            ["dataset", "model", "style", "boot"], []
        )
    missing.to_csv(out / "observed_absent_records.csv", index=False)

    # Run-time summary over all recorded fits
    times = d.copy()
    times["status"] = "success"
    times["n_successful"] = 1
    runtime = summarize(times, ["dataset", "model", "style"], ["elapsed_seconds"], True)
    if a.verify:
        compare(
            runtime, read(ROOT / "results/observed_data/bootstrap_recorded_fit_time_summary.csv"),
            ["dataset", "model", "style", "metric"],
            ["recorded_records", "successful_records", "n", "mean", "median", "low", "high"]
        )
    runtime.to_csv(out / "recorded_runtime_summary.csv", index=False, float_format="%.17g")

    # Semi-synthetic fit coverage
    d = read(ROOT / "results/semi_synthetic/recovery_metrics_by_repetition.csv")

    # Fit counts per condition
    rows = []
    keys = ["dataset", "model", "dropped_modality", "mechanism", "target_rate"]
    for key, g in d.groupby(keys):
        rows.append(
            dict(zip(keys, key))
            | dict(
                records=len(g), expected_repetitions=100,
                successful=int(g.status.eq("success").sum()),
                failed=int(g.status.ne("success").sum()), absent=100 - len(g),
            )
        )

    c = pd.DataFrame(rows)
    if a.verify:
        compare(
            c, read(ROOT / "results/semi_synthetic/recovery_fit_success_counts.csv"), keys,
            [x for x in c if x not in keys]
        )
    c.to_csv(out / "semisynthetic_fit_counts.csv", index=False)

    # Record what was produced and whether it was verified
    outputs["observed_fit_count_conditions"] = len(counts)
    outputs["semisynthetic_fit_count_conditions"] = len(rows)
    outputs["absent_records"] = len(absent)
    (out / "verification.json").write_text(
        json.dumps(
            {
                "summary_rows": outputs, "frozen_summaries_match": a.verify,
                "patient_inputs_used": False, "model_fitting": False,
            },
            indent=2,
        )
        + "\n"
    )
    print("Summaries match the saved tables." if a.verify else "Summaries written.")


if __name__ == "__main__":
    main()
