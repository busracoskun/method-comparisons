"""Run a supplementary R helper in this folder with the configured R installation.

Usage: python standalone_r.py <script.R> [args...]
"""

import json
import os
from pathlib import Path
import subprocess

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "config/default.json").is_file())


def run_r(script, *args):
    # Configured R installation and libraries
    c = json.loads((ROOT / "config/default.json").read_text())
    if (ROOT / "config/host.local.json").exists():
        c.update(json.loads((ROOT / "config/host.local.json").read_text()))
    env = os.environ.copy()
    if c.get("r_home"):
        env["R_HOME"] = c["r_home"]
    if c.get("r_library"):
        env["R_LIBS_USER"] = c["r_library"]
    env["LD_LIBRARY_PATH"] = ":".join(c.get("library_paths", []) + [env.get("LD_LIBRARY_PATH", "")])

    result = subprocess.run(
        [
            c["rscript"],
            "--vanilla",
            str(Path(__file__).resolve().parent / script),
            *[str(a) for a in args],
        ],
        env=env,
        check=False,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(f"R helper {script} failed:\n{result.stdout}\n{result.stderr}")


if __name__ == "__main__":
    import sys

    run_r(sys.argv[1], *sys.argv[2:])
