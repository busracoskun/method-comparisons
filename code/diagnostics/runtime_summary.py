"""Summarize recorded fitting times for the observed-data bootstraps and semi-synthetic tasks."""

from pathlib import Path
import numpy as np
import pandas as pd
import argparse

# Output folder (never inside results/)
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "config/default.json").is_file())
p = argparse.ArgumentParser()
p.add_argument("--output", type=Path, default=ROOT / "build/runtime")
args = p.parse_args()
OUT = args.output.resolve()
assert OUT != ROOT / "results" and ROOT / "results" not in OUT.parents
OUT.mkdir(parents=True, exist_ok=True)
SRC = ROOT


def read(p):
    return pd.read_csv(p, keep_default_na=False, na_values=[""])


# Recorded run times, one row per fit record
obs = read(SRC / "results/observed_data/bootstrap_recorded_fit_times.csv")
sem = read(SRC / "results/semi_synthetic/recovery_metrics_by_repetition.csv")
assert not obs.duplicated(["dataset", "model", "style", "boot"]).any()
assert not sem.duplicated(["dataset", "model", "task_id", "dropped_modality"]).any()

methods = (
    "mofa gfa icluster fabia rgcca sgcca block_pls o2pls gpcca cia mcia jive "
    "gipca bidifac bsf slide intnmf cmf tica nemo"
).split()
styles = ["CC_SF", "MI_SF", "NA"]
rows = []


def stats(v):
    v = pd.to_numeric(v, errors="coerce")
    v = v[np.isfinite(v) & (v > 0)]
    return dict(
        n_timed=len(v), median_seconds=v.median(), q25_seconds=v.quantile(0.25),
        q75_seconds=v.quantile(0.75), p025_seconds=v.quantile(0.025),
        p975_seconds=v.quantile(0.975), mean_seconds=v.mean()
    )


# Observed data: per-fit times for fully successful records, by method and strategy
for ds in ["adni", "tcga"]:
    for m in methods:
        for st in styles:
            g = obs[(obs.dataset == ds) & (obs.model == m) & (obs["style"] == st)]
            if g.empty:
                rows.append(
                    dict(
                        dataset=ds, model=m, strategy=st,
                        coverage="not represented in saved workflow", n_recorded=0, n_absent=np.nan
                    )
                )
                continue

            expected = 5 if st == "MI_SF" else 1
            eligible = (
                (g.status == "success")
                & (g.n_attempted == expected)
                & (g.n_successful == expected)
                & (g.n_failed == 0)
            )
            v = g.loc[eligible, "elapsed_seconds"] / expected
            rows.append(
                dict(
                    dataset=ds,
                    model=m,
                    strategy=st,
                    coverage=(
                        "available"
                        if v.notna().any()
                        else (
                            "no successful fits"
                            if not (g.status == "success").any()
                            else "timing missing"
                        )
                    ),
                    n_recorded=len(g),
                    n_absent=1000 - len(g),
                    n_success=int((g.status == "success").sum()),
                    n_partial=int((g.status == "partial").sum()),
                    n_failed=int((g.status == "failed").sum()),
                    n_missing_time=int(g.elapsed_seconds.isna().sum()),
                    fits_per_record=expected,
                    **stats(v)
                )
            )

summary = pd.DataFrame(rows)
summary.to_csv(OUT / "observed_runtime_summary.csv", index=False)

# Keep failure/partial time distributions in their original per-record units
fail = []
for key, g in obs[obs.status != "success"].groupby(["dataset", "model", "style", "status"]):
    fail.append(
        {
            **dict(zip(["dataset", "model", "strategy", "status"], key)),
            "n_records": len(g),
            "n_attempted_total": g.n_attempted.sum(),
            "n_successful_total": g.n_successful.sum(),
            "n_failed_total": g.n_failed.sum(),
            **stats(g.elapsed_seconds),
        }
    )
pd.DataFrame(fail).to_csv(OUT / "failed_partial_record_times.csv", index=False)

# Recovery timings: stratify by modality/rate, pool the three equally scheduled mechanisms
keys = ["dataset", "model", "dropped_modality", "target_rate"]
srows = []
for key, g in sem.groupby(keys):
    good = g[g.status == "success"]
    srows.append(
        dict(
            zip(keys, key), n_recorded=len(g), n_success=len(good),
            n_failed=int((g.status != "success").sum()), **stats(good.fit_elapsed_seconds)
        )
    )
ss = pd.DataFrame(srows)
ss.to_csv(OUT / "semisynthetic_runtime_summary.csv", index=False)

# Per mechanism details avoid hiding condition-specific runtime
detail = []
for key, g in sem.groupby(keys + ["mechanism"]):
    detail.append(
        dict(
            zip(keys + ["mechanism"], key), n_recorded=len(g),
            n_success=int((g.status == "success").sum()),
            **stats(g.loc[g.status == "success", "fit_elapsed_seconds"])
        )
    )
pd.DataFrame(detail).to_csv(OUT / "semisynthetic_runtime_by_mechanism.csv", index=False)

# Deduplicate imputation only within engine: Python and R used independent pipelines
sem["engine"] = np.where(sem.model.isin(["mofa", "cmf"]), "python", "r")
ik = ["dataset", "engine", "task_id", "dropped_modality"]
assert (sem.groupby(ik).imputation_elapsed_seconds.nunique(dropna=True) <= 1).all()
imp = sem.drop_duplicates(ik)
imp[ik + ["mechanism", "target_rate", "imputation_elapsed_seconds"]].to_csv(
    OUT / "unique_scenario_imputation_times.csv", index=False
)

ir = []
for key, g in imp.groupby(["dataset", "engine", "dropped_modality", "target_rate"]):
    ir.append(
        dict(
            zip(["dataset", "engine", "dropped_modality", "target_rate"], key),
            **stats(g.imputation_elapsed_seconds)
        )
    )
pd.DataFrame(ir).to_csv(OUT / "scenario_imputation_summary.csv", index=False)
