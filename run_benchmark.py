#!/usr/bin/env python3
"""Run observed-data or semi-synthetic experiments using authorized local data."""

import argparse
import datetime
import json
import os
import subprocess
from pathlib import Path

from project_paths import (
    authorized_input_path,
    check_output_tree,
    input_file,
    load_config,
    local_path,
)

ROOT = Path(__file__).resolve().parent
R_METHODS = {
    "rgcca", "gpcca", "block_pls", "gipca", "bsf", "nemo",
    "bidifac", "icluster", "fabia", "mcia", "intnmf", "gfa",
    "sgcca", "o2pls", "cia", "tica", "slide", "jive",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pipeline",
        choices=["adni-observed", "tcga-observed", "adni-semisynthetic", "tcga-semisynthetic"],
    )
    parser.add_argument("--config", type=Path, help="Optional configuration overrides")
    parser.add_argument("--methods", help="Comma-separated names; order affects original RNG streams")
    parser.add_argument("--engine", choices=["all", "r", "python"], default="all")
    parser.add_argument("--mode", choices=[
        "bootstrap", "full", "prepare", "rescore", "ground_truth", "run", "combine",
    ])
    parser.add_argument("--count", type=int)
    parser.add_argument("--start", type=int)
    parser.add_argument("--task", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true", help="Check inputs without fitting")
    parser.add_argument(
        "--tcga-intnmf-repair", action="store_true",
        help="Use the documented TCGA-only intNMF repair seed index",
    )
    return parser, parser.parse_args()


def select_methods(parser, args, config, dataset, semi):
    key = "semisynthetic_methods" if semi else "observed_methods"
    methods = args.methods.split(",") if args.methods else config[key]

    if args.tcga_intnmf_repair and (dataset != "tcga" or not semi or methods != ["intnmf"]):
        parser.error("Repair profile requires tcga-semisynthetic --methods intnmf")
    if not methods or len(set(methods)) != len(methods) or set(methods) - R_METHODS - {"mofa", "cmf"}:
        parser.error("Unknown or duplicate methods")

    r_methods = [name for name in methods if name in R_METHODS]
    python_methods = [name for name in methods if name in {"mofa", "cmf"}]
    if args.engine == "r":
        python_methods = []
    if args.engine == "python":
        r_methods = []
    return methods, r_methods, python_methods


def prepare_environment(args, config, dataset, mode, r_methods, python_methods, count, start):
    output = local_path(args.output or config["output_dir"])
    work = check_output_tree(output / dataset)
    adni_input = authorized_input_path(config["adni_input_dir"])
    tcga_input = authorized_input_path(config["tcga_input_dir"])
    cache_input = (
        authorized_input_path(config["cache_input_dir"])
        if config.get("cache_input_dir") else None
    )

    env = os.environ.copy()
    # Remove legacy paths that might point to a different run or dataset.
    for key in (
        "BENCH_DIR", "DATA_DIR", "RESULTS_FULL_DIR", "FULL_MODEL_DIR",
        "PYTHON_METHODS_FILE", "FINAL_R_FILE", "FINAL_TCGA_FILE", "PREPROCESSED_DIR",
    ):
        env.pop(key, None)

    env.update({
        "CLEAN_ROOT": str(ROOT),
        "BENCH_WORK_DIR": str(work),
        "BENCH_DATASET": dataset,
        "ADNI_INPUT_DIR": str(adni_input),
        "TCGA_INPUT_DIR": str(tcga_input),
        "B_BOOT": str(count),
        "BOOT_START": str(start),
        "BOOT_CORES": str(config["cores"]),
        "LSB_DJOB_NUMPROC": str(config["cores"]),
        "SAVE_BOOT_MODELS": "1" if config["save_models"] else "0",
        "OBSERVED_MODE": mode,
        "R_METHODS": ",".join(r_methods),
        "PY_METHODS": ",".join(python_methods),
        "TCGA_INTNMF_REPAIR": "1" if args.tcga_intnmf_repair else "0",
        "SEMISYNTH_METHODS": ",".join(r_methods or ["icluster"]),
        "SEMISYNTH_MODE": mode,
        "SEMISYNTH_TASK_ID": str(args.task or 0),
        "PY_SEMISYNTH_MODE": mode,
        "PY_SEMISYNTH_TASK_ID": str(args.task or 0),
        "SEMISYNTH_OUT_DIR": str(work / "semisynthetic_r"),
        "PY_SEMISYNTH_OUT_DIR": str(work / "semisynthetic_python"),
        "PROJECT_ROOT": str(ROOT),
        "TMPDIR": str(ROOT / "private_outputs/tmp"),
        "PYTHONPATH": os.pathsep.join([str(ROOT), str(ROOT / "code/methods")]),
        "PYTHONDONTWRITEBYTECODE": "1",
    })
    if config.get("r_home"):
        env["R_HOME"] = config["r_home"]
    if config.get("r_library"):
        env["R_LIBS_USER"] = config["r_library"]

    env["CACHE_INPUT_DIR"] = str(cache_input) if cache_input else ""
    env["TMP"] = env["TEMP"] = env["TMPDIR"]
    env["JOBLIB_TEMP_FOLDER"] = str(ROOT / "private_outputs/tmp/joblib")
    env["XDG_CACHE_HOME"] = str(ROOT / "private_outputs/.cache")
    env["MPLCONFIGDIR"] = str(ROOT / "private_outputs/.cache/matplotlib")
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        config.get("library_paths", []) + [env.get("LD_LIBRARY_PATH", "")]
    )
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env.setdefault(name, "1")
    return work, env


def build_commands(dataset, mode, methods_r, methods_python, args, config, count, start):
    r = [config["rscript"], "--vanilla"]
    python = [config["python"]]
    commands = []
    semi = args.pipeline.endswith("semisynthetic")

    if mode == "rescore":
        commands.append(r + [str(ROOT / "code/observed_data/evaluate_saved_adni_fits.R")])
    elif semi:
        if methods_r:
            commands.append(r + [str(ROOT / "code/semi_synthetic/run_recovery_benchmark.R")])
        if methods_python:
            commands.append(python + [str(ROOT / "code/semi_synthetic/run_mofa_cmf_recovery.py")])
    else:
        if methods_r or mode == "prepare" or (dataset == "adni" and methods_python):
            # ADNI Python fitting uses variables and imputations prepared by R.
            override = {"OBSERVED_MODE": "prepare"} if not methods_r or mode == "prepare" else {}
            commands.append((r + [str(ROOT / "code/observed_data/run_observed_benchmark.R")], override))

        if methods_python and mode != "prepare":
            if dataset == "adni" and mode == "bootstrap":
                for boot in range(start, start + count):
                    commands.append((
                        python + [str(ROOT / "code/observed_data/run_adni_mofa_cmf_bootstrap.py")],
                        {"BOOT_ID": str(boot)},
                    ))
            else:
                script = (
                    "code/observed_data/run_adni_mofa_cmf_full_cohort.py"
                    if dataset == "adni" else "code/observed_data/run_tcga_mofa_cmf.py"
                )
                commands.append(python + [str(ROOT / script)])

    if dataset == "adni" and mode == "full" and config["save_models"]:
        commands.append(r + [str(ROOT / "code/observed_data/evaluate_saved_adni_fits.R")])
    return commands


def ensure_inputs(parser, dataset, env):
    names = (
        ["meta.csv", "rnaseq.csv", "clin.csv"]
        if dataset == "adni" else
        ["irsem_fs.csv", "imeth_fs.csv", "irppa_fs.csv", "clin_common.csv"]
    )
    folder = Path(env["ADNI_INPUT_DIR" if dataset == "adni" else "TCGA_INPUT_DIR"])
    paths = [input_file(folder, name) if dataset == "adni" else folder / name for name in names]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        parser.error("Missing input files: " + ", ".join(missing))


def main():
    parser, args = parse_args()
    config = load_config(args.config)
    dataset, analysis = args.pipeline.split("-", 1)
    semi = analysis == "semisynthetic"
    mode = args.mode or ("run" if semi else "bootstrap")

    allowed = ["ground_truth", "run", "combine"] if semi else ["bootstrap", "full", "prepare", "rescore"]
    if mode not in allowed:
        parser.error("Mode does not apply to this pipeline")
    if semi and mode == "run" and not (args.task and 1 <= args.task <= 600):
        parser.error("--task must be 1..600")
    if config.get("repetitions", 100) != 100:
        parser.error("The preserved design requires 100 repetitions; use task subsets for smaller runs")
    if mode == "rescore" and dataset != "adni":
        parser.error("The historical saved-fit rescoring correction is ADNI-specific")

    methods, methods_r, methods_python = select_methods(parser, args, config, dataset, semi)
    count = args.count if args.count is not None else config["bootstrap_count"]
    start = args.start if args.start is not None else config["bootstrap_start"]
    if count < 1 or start < 1 or config["cores"] < 1:
        parser.error("count, start and cores must be positive")

    try:
        work, env = prepare_environment(
            args, config, dataset, mode, methods_r, methods_python, count, start
        )
    except ValueError as error:
        parser.error(str(error))

    commands = build_commands(dataset, mode, methods_r, methods_python, args, config, count, start)
    for entry in commands:
        command = entry[0] if isinstance(entry, tuple) else entry
        print(command, flush=True)
    ensure_inputs(parser, dataset, env)
    if args.check:
        print("Configuration and input files OK. No fitting performed.")
        return

    work.mkdir(parents=True, exist_ok=True)
    for key in ("TMPDIR", "JOBLIB_TEMP_FOLDER", "XDG_CACHE_HOME", "MPLCONFIGDIR"):
        local_path(env[key]).mkdir(parents=True, exist_ok=True)
    logs = work / "logs"
    logs.mkdir(exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    record = {"arguments": {**vars(args), "config": str(args.config), "output": str(args.output or config["output_dir"])},
              "config": config, "methods": methods}
    (logs / f"{stamp}_config.json").write_text(json.dumps(record, default=str, indent=2))

    for i, entry in enumerate(commands):
        command, overrides = entry if isinstance(entry, tuple) else (entry, {})
        with (logs / f"{stamp}_{i}.log").open("w") as log:
            result = subprocess.run(
                command, cwd=work, env={**env, **overrides},
                stdout=log, stderr=subprocess.STDOUT,
            )
        if result.returncode:
            raise SystemExit(f"Command failed ({result.returncode}). See {log.name}")
    print(f"Completed {args.pipeline} {mode}. Outputs: {work}")


if __name__ == "__main__":
    main()
