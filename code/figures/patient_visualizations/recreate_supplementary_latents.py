#!/usr/bin/env python3
"""Optional: recreate MCIA/GFA latent scores for the selected MCAR repetitions (Figures S1, S2).

Requires the repetition selection (mcia_gfa_umap_selection.csv) derived from the metric tables.
Existing score CSVs are reused. Runs recreate_supplementary_latents.R with the configured R
installation; it does not run run_benchmark.py or the full benchmark.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "config/default.json").is_file())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--selection", type=Path,
        default=ROOT / "results/semi_synthetic/mcia_gfa_umap_selection.csv"
    )
    p.add_argument(
        "--latents", type=Path, default=ROOT / "private_outputs/supplementary/latent_scores"
    )
    p.add_argument("--results", type=Path, default=ROOT / "private_outputs")
    a = p.parse_args()
    import pandas as pd

    # Check the selection and output location
    if not a.selection.exists():
        p.error("Generate representative-repetition selection from real metric tables first")
    d = pd.read_csv(a.selection)
    if d.empty:
        p.error("No representative repetitions selected; metric inputs are missing")
    if not set(d.model).issubset({"mcia", "gfa"}):
        p.error("Only predeclared MCIA/GFA methods supported")
    if not set(d.dataset).issubset({"adni", "tcga"}):
        p.error("Invalid dataset")
    if not a.latents.resolve().is_relative_to(ROOT / "private_outputs"):
        p.error("Write generated latent scores under private_outputs/.")

    # Configuration and environment for the R helper
    config = json.loads((ROOT / "config/default.json").read_text())
    host = ROOT / "config/host.local.json"
    if host.exists():
        config.update(json.loads(host.read_text()))
    methods = [m for m in config["semisynthetic_methods"] if m not in ["mofa", "cmf"]]
    env = os.environ.copy()
    env.update(
        CLEAN_ROOT=str(ROOT), PROJECT_ROOT=str(ROOT), R_METHODS=",".join(methods),
        SEMISYNTH_METHODS=",".join(methods), TCGA_INTNMF_REPAIR="0", SAVE_SEMISYNTH_DATA="0",
        RETRY_FAILED="0", SEMISYNTH_REPETITIONS="100", SEMISYNTH_N_FACTORS="3", SEMISYNTH_BETA="2",
        MISSFOREST_NTREE="100", MISSFOREST_MAXITER="10",
        BENCH_WORK_DIR=str((a.latents / "work").resolve()),
        SEMISYNTH_OUT_DIR=str((a.latents / "work/semisynthetic").resolve()),
        ADNI_INPUT_DIR=str(ROOT / config["adni_input_dir"]),
        TCGA_INPUT_DIR=str(ROOT / config["tcga_input_dir"]),
        CACHE_INPUT_DIR=str(ROOT / config["cache_input_dir"])
    )
    if config.get("r_home"):
        env["R_HOME"] = config["r_home"]
    if config.get("r_library"):
        env["R_LIBS_USER"] = config["r_library"]
    env["LD_LIBRARY_PATH"] = ":".join(
        config.get("library_paths", []) + [env.get("LD_LIBRARY_PATH", "")]
    )

    subprocess.run(
        [
            config["rscript"],
            "--vanilla",
            str(ROOT / "code/figures/patient_visualizations/recreate_supplementary_latents.R"),
            str(a.selection.resolve()),
            str(a.latents.resolve()),
            str(a.results.resolve()),
        ],
        env=env,
        cwd=ROOT,
        check=True,
    )


if __name__ == "__main__":
    main()
