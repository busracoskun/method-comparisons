#!/usr/bin/env python3
"""Select the tICA MCAR repetition shown in Figure 7 for each dataset and dropped modality.

The same repetition is used at both rates: the one whose RV ranks are closest to the median
(sum of absolute percentile-rank distances). Selection uses only the numerical results, before
any embedding is inspected. The generated combined selection is written under ``build/`` and
does not replace the published per-modality selection records under ``results/``.
"""

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from project_paths import authorized_input_path, check_output_tree  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adni-results",
        type=Path,
        default=ROOT / "results/semi_synthetic/adni_semi_synthetic_results.csv",
    )
    parser.add_argument(
        "--tcga-results",
        type=Path,
        default=ROOT / "results/semi_synthetic/tcga_semi_synthetic_results.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "build/patient_visualizations/semisynthetic_UMAP_selection.csv",
        help="Generated combined selection CSV under build/.",
    )
    args = parser.parse_args()

    output = args.output.resolve()
    if not output.is_relative_to(ROOT / "build"):
        parser.error("Choose an output path under build/.")
    output_dir = check_output_tree(output.parent)

    import pandas as pd

    inputs = {
        "ADNI": (authorized_input_path(args.adni_results), ["gene", "meta"]),
        "TCGA": (authorized_input_path(args.tcga_results), ["methylation", "rppa"]),
    }
    required = {
        "model", "mechanism", "dropped_modality", "status", "repetition", "target_rate", "rv"
    }
    rows = []
    for dataset, (source, modalities) in inputs.items():
        if not source.is_file():
            parser.error(f"Missing semi-synthetic result file: {source}")
        data = pd.read_csv(source)
        missing = sorted(required - set(data.columns))
        if missing:
            parser.error(f"{source} is missing required columns: {', '.join(missing)}")

        for modality in modalities:
            group = data[
                (data.model == "tica")
                & (data.mechanism == "MCAR")
                & (data.dropped_modality == modality)
                & (data.status == "success")
            ]
            pivoted = group.pivot(
                index="repetition", columns="target_rate", values="rv"
            ).dropna()
            if pivoted.empty:
                parser.error(f"No complete tICA MCAR repetitions for {dataset} {modality}")

            # Preserve the historical rank-distance criterion and stable tie ordering exactly.
            score = (pivoted.rank(pct=True) - 0.5).abs().sum(axis=1)
            repetition = score.sort_values(kind="stable").index[0]
            selected = group[group.repetition == repetition].copy()
            selected["selection_distance"] = score.loc[repetition]
            selected["plot_dataset"] = dataset
            rows.append(selected)

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(rows).to_csv(output, index=False)
    print(f"Wrote {len(rows)} dataset/modality selections to {output}")


if __name__ == "__main__":
    main()
