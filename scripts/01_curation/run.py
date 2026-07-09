#!/usr/bin/env python3
# Purpose: Run the workflow stage.

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "fetch_uniprot.py",
    "curate_master.py",
    "merge_cazy.py",
    "merge_brenda.py",
    "merge_metadata.py",
    "build_master.py",
    "export_groups.py",
    "pdb_inventory.py",
    "download_pdb.py",
    "refseq_inventory.py",
    "fetch_refseq.py",
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
