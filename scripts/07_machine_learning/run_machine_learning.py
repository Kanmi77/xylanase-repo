#!/usr/bin/env python3
"""
Run the machine-learning stage.

This stage covers:
- experimental Topt/Tm label preparation
- 60C-aligned label correction
- sequence/structure mapping
- sequence/structure classifier training
- feature importance extraction
- WT and mutation-prioritisation panels
"""

from pathlib import Path
import subprocess
import sys


BASE_DIR = Path(__file__).resolve().parents[2]
STEP_DIR = Path(__file__).resolve().parent

STEPS = [
    "01_prepare_tmxyl_toptxyl_experimental_labels.py",
    "02_train_experimental_sequence_classifier.py",
    "03_match_external_temperatures_to_master.py",
    "04_merge_60c_aligned_internal_external_and_retrain.py",
    "05_build_sequence_structure_mapping.py",
    "06_train_sequence_structure_classifier.py",
    "extract_selected_extratrees_feature_importance.py",
    "00_review_stability_ml_inputs.py",
    "02_build_wt_stability_panel.py",
    "05_train_wt_stability_ml.py",
    "07_build_mutation_prioritisation_panel.py",
]


def run_step(script_name: str) -> None:
    script_path = STEP_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing workflow step: {script_path}")

    print(f"\n[MACHINE LEARNING] Running {script_name}")
    subprocess.run([sys.executable, str(script_path)], check=True)


def main() -> int:
    for step in STEPS:
        run_step(step)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
