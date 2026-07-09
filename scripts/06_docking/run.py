#!/usr/bin/env python3
# Purpose: Run the workflow stage.

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "clean_mutants.py",
    "run_mutants.py",
    "run_vina.py",
    "parse_mutants.py",
    "prepare_wt.py",
    "compare_wt.py",
    "summarise_pdb60.py",
    "summarise_modeller.py",
    "verify_modeller.py",
    "verify_ligands.py",
    "verify_pdb.py",
]


def run_step(script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[DOCKING] Running {script_name}")
    subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for step in STEPS:
        run_step(step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
