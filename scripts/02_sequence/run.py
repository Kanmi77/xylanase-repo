#!/usr/bin/env python3
# Purpose: Run the workflow stage.

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    ("python", "compute_protparam.py"),
    ("bash", "run_mafft.sh"),
    ("bash", "run_fasttree.sh"),
    ("python", "add_temperature.py"),
    ("python", "summarise_conservation.py"),
    ("python", "map_conservation.py"),
    ("python", "package_pymol.py"),
]


def run_step(mode: str, script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[SEQUENCE/PHYLOGENY] Running {script_name}")
    if mode == "bash":
        subprocess.run(["bash", str(script_path)], check=True)
    else:
        subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for mode, step in STEPS:
        run_step(mode, step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
