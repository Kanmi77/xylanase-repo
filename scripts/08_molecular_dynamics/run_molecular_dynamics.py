#!/usr/bin/env python3
"""
Run the molecular-dynamics summary stage.

This stage covers:
- extraction of generated MD metrics
- primary WT/mutant MD tables
- DSSP-supported secondary-structure summaries
- thesis-ready MD metric tables
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    ("bash", "01_collect_md_result_files.sh"),
    ("python", "10_extract_md_metrics_from_original_files.py"),
    ("python", "11_make_primary_md_result_tables.py"),
    ("python", "18_make_thesis_md_metric_tables.py"),
    ("python", "19_extract_all_md_generated_metrics.py"),
    ("python", "09_make_md_24_simulation_summary_with_dssp.py"),
    ("python", "10_make_md_24_summary_standalone_mkdssp.py"),
]


def run_step(mode: str, script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[MOLECULAR DYNAMICS] Running {script_name}")
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
