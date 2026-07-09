#!/usr/bin/env python3
# Purpose: Run the workflow stage.

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "structure_manifest.py",
    "standardise_structures.py",
    "compute_features.py",
    "prepare_tmalign.py",
    "run_tmalign.py",
    "parse_tmalign.py",
    "parse_tmalign_full.py",
    "summarise_foldx.py",
    "salt_bridges.py",
    "compare_sources.py",
    "summarise_groups.py",
    "summarise_features.py",
]


def run_step(script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[STRUCTURE ANALYSIS] Running {script_name}")
    subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for step in STEPS:
        run_step(step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
