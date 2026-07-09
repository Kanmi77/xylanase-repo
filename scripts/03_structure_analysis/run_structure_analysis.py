#!/usr/bin/env python3
"""
Run the structural analysis stage.

This stage covers:
- combined PDB/MODELLER structure manifest construction
- structure standardisation
- structural feature extraction
- TM-align preparation, execution and parsing
- source-stratified structure summaries
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "01_build_combined_structure_manifest.py",
    "02_standardize_combined_structures.py",
    "03_compute_combined_structural_features.py",
    "04_prepare_tmalign_jobs.py",
    "05_run_tmalign.py",
    "06_extract_tmalign_best_reference_per_model.py",
    "parse_tmalign_outputs_full.py",
    "summarise_structural_features_foldx.py",
    "calculate_salt_bridges_from_foldx_structures.py",
    "analyse_pdb_and_modeller_results.py",
    "analyse_structural_features_separately.py",
    "analyse_structural_features_without_foldx.py",
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
