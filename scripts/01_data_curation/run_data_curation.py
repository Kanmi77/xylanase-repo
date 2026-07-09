#!/usr/bin/env python3
"""
Run the data curation stage for the xylanase thermostability workflow.

This stage covers:
- UniProt xylanase retrieval
- GH10/GH11 curation
- CAZy, BRENDA, PDB and RefSeq metadata integration
- bacterial/fungal/GH-family group export
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "01_fetch_uniprot_xylanase_tsv.py",
    "02_curate_uniprot_to_master_csv.py",
    "08_refetch_uniprot_add_cazy_and_merge.py",
    "08_merge_brenda_into_master.py",
    "08_merge_experimental_metadata.py",
    "09_build_master_all_curated.py",
    "11_export_gh10_gh11_fungi_bacteria_files.py",
    "04_build_pdb_inventory.py",
    "05_download_pdb_structures.py",
    "06_build_refseq_inventory.py",
    "07_fetch_refseq_fasta.py",
]


def run_step(script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[DATA CURATION] Running {script_name}")
    subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for step in STEPS:
        run_step(step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
