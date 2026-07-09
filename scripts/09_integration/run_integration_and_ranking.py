#!/usr/bin/env python3
"""
Run the integration and candidate-ranking stage.

This stage covers:
- thesis master-table integration
- salt-bridge merge
- integrated scoring
- PDB vs homology source-stratified ranking
- metadata annotation
- final visual outputs
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "01_build_thesis_master_table.py",
    "04_compute_and_merge_salt_bridges_all_structures.py",
    "01_integrate_scoring.py",
    "21_stratified_pdb_vs_homology_integration.py",
    "23_prepare_homology_full_docking_candidate_level.py",
    "24_annotate_integrated_candidates_with_master_metadata.py",
    "25_fill_integrated_candidate_metadata_fallbacks.py",
    "02_visualize_integrated_results.py",
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
