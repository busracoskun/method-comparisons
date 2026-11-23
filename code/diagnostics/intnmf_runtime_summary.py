"""Summarize intNMF fitting times from the sensitivity analysis (no fitting).

These runs are kept separate: the repaired TCGA intNMF masks are not pooled with the main
semi-synthetic experiment.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd

# Output folder (under build/)
ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT
p = argparse.ArgumentParser()
p.add_argument("--output", type=Path, default=ROOT / "build/intnmf_runtime")
a = p.parse_args()
OUT = a.output.resolve()
if not OUT.is_relative_to(ROOT / "build"):
    p.error("Choose output under build/")
OUT.mkdir(parents=True, exist_ok=True)


def read(p):
    return pd.read_csv(p, keep_default_na=False, na_values=[""])


def stats(v):
    v = pd.to_numeric(v, errors="coerce")
    v = v[np.isfinite(v) & (v > 0)]
    return dict(
        n_timed=len(v), median_seconds=v.median(), q25_seconds=v.quantile(0.25),
        q75_seconds=v.quantile(0.75), p025_seconds=v.quantile(0.025),
        p975_seconds=v.quantile(0.975), mean_seconds=v.mean()
    )


# Run-time summary per experiment and condition
sp = SRC / "results/semi_synthetic/intnmf_sensitivity_metrics_by_repetition.csv"
sens = read(sp)
sr = []
keys = ["dataset", "experiment", "model", "dropped_modality", "target_rate", "mechanism"]
for key, group in sens.groupby(keys):
    successful = group.loc[group.status == "success", "fit_elapsed_seconds"]
    sr.append(
        {
            **dict(zip(keys, key)),
            "n_recorded": len(group),
            "n_success": int((group.status == "success").sum()),
            "n_failed": int((group.status != "success").sum()),
            **stats(successful),
        }
    )
pd.DataFrame(sr).to_csv(OUT / "intnmf_sensitivity_runtime_summary.csv", index=False)
