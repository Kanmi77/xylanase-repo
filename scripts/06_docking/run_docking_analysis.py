#!/usr/bin/env python3
"""
Run the docking-analysis stage.

This stage covers:
- mutant receptor preparation
- PDBQT preparation
- AutoDock Vina score parsing
- WT-mutant docking delta calculation
- PDB and MODELLER docking summary tables
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "08_clean_pdb60_mutants_prepare_docking.py",
    "09_make_pdbqt_fallback_run_pdb60_mutant_docking.py",
    "10_make_vina_compatible_pdbqt_and_run.py",
    "11_parse_pdb60_mutant_vina_scores.py",
    "12_prepare_matching_wt_docking_for_pdb60_mutants.py",
    "13_parse_pdb60_wt_and_calculate_docking_delta.py",
    "summarise_pdb60_docking_results.py",
    "summarise_modeller_docking_results.py",
    "make_verified_modeller_docking_tables.py",
    "make_verified_ligand_specific_docking_tables.py",
    "make_verified_pdb_ligand_specific_docking_tables.py",
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
