#!/usr/bin/env python3
"""
Run the mutation-screening preparation and parsing stage.

This stage covers:
- PDB mutation/docking reconciliation
- PDB60 consensus mutation preparation
- FoldX mutation parsing
- candidate-level mutation table generation
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "03_reconcile_pdb_mutation_docking.py",
    "05_prepare_pdb60_consensus_foldx_mutations.py",
    "06_parse_pdb60_foldx_prepare_mutant_docking.py",
    "16_make_pdb60_final_ranked_no_tiers.py",
]


def run_step(script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[MUTATION SCREENING] Running {script_name}")
    subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for step in STEPS:
        run_step(step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
