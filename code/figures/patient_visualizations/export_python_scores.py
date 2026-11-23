#!/usr/bin/env python3
"""Export TCGA MOFA full-run scores (NA and MI_SF imputation 1) with patient IDs.

This helper requires the private patient-level fitted-model artifacts; aggregate result tables
cannot be substituted. Pass the directory that contains the ``mofa/`` fit subdirectory with
``--fits-dir``. Generated patient-level CSVs are written under ``private_outputs/``.
"""

import argparse
from pathlib import Path
import pickle
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from project_paths import authorized_input_path, check_output_tree  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fits-dir",
        type=Path,
        required=True,
        help=(
            "Trusted private fit directory containing mofa/mofa_NA_0.pkl and "
            "mofa/mofa_MI_SF_0_imp1.pkl (for example, private_outputs/tcga/fits_python)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "private_outputs/patient_visualizations/python_scores",
        help="Patient-level output directory under private_outputs/.",
    )
    args = parser.parse_args()

    fits_dir = authorized_input_path(args.fits_dir)
    output_dir = check_output_tree(args.output_dir)
    if not output_dir.is_relative_to(ROOT / "private_outputs"):
        parser.error("Choose a patient-level output directory under private_outputs/.")

    sources = {
        "NA": fits_dir / "mofa/mofa_NA_0.pkl",
        "MI_SF": fits_dir / "mofa/mofa_MI_SF_0_imp1.pkl",
    }
    missing = [str(path) for path in sources.values() if not path.is_file()]
    if missing:
        parser.error("Missing required patient-level fitted model(s): " + ", ".join(missing))

    import numpy as np
    import pandas as pd

    output_dir.mkdir(parents=True, exist_ok=True)
    for strategy, source in sources.items():
        # Pickle is unsafe for untrusted input. --fits-dir must point to authorized artifacts.
        with source.open("rb") as f:
            fit = pickle.load(f)
        z = fit["Z"]
        assert (
            isinstance(z, pd.DataFrame)
            and z.index.is_unique
            and np.isfinite(z.to_numpy()).all()
        )
        assert z.index.astype(str).tolist() == list(map(str, fit["rows"]))
        z.rename_axis("patient_id").to_csv(
            output_dir / f"TCGA_mofa_{strategy}_latest_Z.csv"
        )

    print(f"Exported {len(sources)} patient-level score files to {output_dir}")


if __name__ == "__main__":
    main()
