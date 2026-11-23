#!/usr/bin/env python3
"""Annotate the full-cohort metric estimates for publication; the raw source CSV is unchanged.

Each metric is labelled as recorded, not available, or undefined (missingness clustering for
CC_SF), together with the number of finite saved runs behind it.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from summarize_results import read, OBS

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "config/default.json").is_file())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--output", type=Path, default=ROOT / "build/full_cohort_metrics_for_publication.csv"
    )
    a = p.parse_args()

    # Historical full-cohort metrics and their finite-run counts
    d = read(ROOT / "results/observed_data/full_cohort_historical_metrics_wide.csv")
    den = pd.concat(
        [
            read(ROOT / f"results/observed_data/{ds}_full_cohort_finite_fit_counts_and_means.csv")
            for ds in ["adni", "tcga"]
        ],
        ignore_index=True,
    )

    # One row per dataset, method, strategy, and metric
    rows = []
    for r in d.to_dict("records"):
        for metric in OBS:
            x = pd.to_numeric(r.get(metric), errors="coerce")
            undefined = r["style"] == "CC_SF" and metric.endswith("_miss")
            q = den.loc[
                den.dataset.eq(r["dataset"])
                & den.model.eq(r["model"])
                & den["style"].eq(r["style"])
                & den.metric.eq(metric)
            ]
            rows.append(
                dict(
                    dataset=r["dataset"], model=r["model"], style=r["style"], metric=metric,
                    value=np.nan if undefined else x,
                    availability=(
                        "undefined_complete_case_missingness"
                        if undefined
                        else "recorded" if np.isfinite(x) else "not_available"
                    ),
                    finite_saved_metric_runs=(
                        0 if undefined else int(q.iloc[0].n) if len(q) == 1 else np.nan
                    ),
                    denominator_source=(
                        "not_recoverable_from_aggregate"
                        if q.empty and not undefined
                        else "raw_saved_scalar_runs" if not undefined else "undefined"
                    ),
                    historical_source_retained=(
                        "results/observed_data/full_cohort_historical_metrics_wide.csv"
                    ),
                )
            )

    a.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.output, index=False, float_format="%.17g")
    print(f"Annotated {len(rows)} full-cohort metric cells; raw historical source unchanged.")


if __name__ == "__main__":
    main()
