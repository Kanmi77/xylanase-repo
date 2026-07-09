#!/usr/bin/env python3
"""
Run the sequence analysis and phylogeny stage.

This stage covers:
- ProtParam-style sequence feature calculation
- MAFFT alignment
- FastTree phylogeny
- conservation and thermal-position summaries
- mapping conserved positions to representative structures
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    ("python", "04_compute_protparam_features.py"),
    ("bash", "02_run_mafft_alignments.sh"),
    ("bash", "03_run_fasttree.sh"),
    ("python", "07_add_brenda_raw_temperature_to_master_692.py"),
    ("python", "summarise_alignment_conservation.py"),
    ("python", "map_conserved_positions_to_structure.py"),
    ("python", "package_pymol_structure_mapping_for_git.py"),
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
