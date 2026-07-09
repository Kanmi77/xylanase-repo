#!/usr/bin/env python3
# Purpose: Run the workflow stage.

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "build_master.py",
    "merge_saltbridges.py",
    "integrate_scores.py",
    "stratify_sources.py",
    "prepare_homology.py",
    "annotate_candidates.py",
    "fill_metadata.py",
    "make_figures.py",
]


def run_step(script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[INTEGRATION] Running {script_name}")
    subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for step in STEPS:
        run_step(step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
