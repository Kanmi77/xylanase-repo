#!/usr/bin/env python3
# Purpose: Package PyMOL mapping files.

from pathlib import Path
import shutil
import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"
EVIDENCE_DIR = Path.home() / "thesis-evidence-pack"

SOURCE_DIR = PROJECT_DIR / "results/thermal_sequence_analysis/structure_mapping"

PACKAGE_DIR = (
    EVIDENCE_DIR
    / "12_thermal_sequence_analysis/results/structure_mapping"
)

PDB_DIR = PACKAGE_DIR / "pdb"
PML_DIR = PACKAGE_DIR / "pymol_scripts"
TABLE_DIR = PACKAGE_DIR / "tables"

SUMMARY_FILE = SOURCE_DIR / "structure_mapping_summary.csv"

PDB_DIR.mkdir(parents=True, exist_ok=True)
PML_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def copy_mapping_tables():
    for csv_file in SOURCE_DIR.glob("*.csv"):
        shutil.copy2(csv_file, TABLE_DIR / csv_file.name)


def residue_selection(name, positions):
    positions = sorted(set(int(pos) for pos in positions if pd.notna(pos)))

    if not positions:
        return f"select {name}, none\n"

    position_text = "+".join(str(pos) for pos in positions)

    return f"select {name}, protein_model and resi {position_text}\n"


def clean_positions(table, column):
    if table.empty or column not in table.columns:
        return []

    positions = []

    for value in table[column].dropna():
        try:
            position = int(value)
        except ValueError:
            continue

        if position > 0:
            positions.append(position)

    return sorted(set(positions))


def write_windows_friendly_pml(dataset, pdb_filename, conserved, contrast):
    pml_file = PML_DIR / f"{dataset}_3d_mapping_view.pml"

    conserved_positions = clean_positions(
        conserved.head(150),
        "reference_residue_position",
    )

    contrast_positions = clean_positions(
        contrast,
        "reference_residue_position",
    )

    if not contrast.empty and "high_confidence" in contrast.columns:
        high_confidence_positions = clean_positions(
            contrast[contrast["high_confidence"] == True],
            "reference_residue_position",
        )
    else:
        high_confidence_positions = []

    with pml_file.open("w") as handle:
        handle.write(f"load ../pdb/{pdb_filename}, protein_model\n")
        handle.write("hide everything\n")
        handle.write("show cartoon, protein_model\n")
        handle.write("color gray80, protein_model\n")
        handle.write("set cartoon_transparency, 0.10, protein_model\n")
        handle.write("bg_color white\n")
        handle.write("set antialias, 2\n")
        handle.write("set sphere_scale, 0.35\n")
        handle.write("set stick_radius, 0.18\n")
        handle.write("set label_size, 16\n")
        handle.write("set label_color, black\n")
        handle.write("\n")

        handle.write(residue_selection("conserved_90", conserved_positions))
        handle.write("color cyan, conserved_90\n")
        handle.write("show sticks, conserved_90\n")
        handle.write("\n")

        handle.write(residue_selection("thermal_contrast", contrast_positions))
        handle.write("color orange, thermal_contrast\n")
        handle.write("show spheres, thermal_contrast\n")
        handle.write("show sticks, thermal_contrast\n")
        handle.write("\n")

        handle.write(
            residue_selection(
                "high_confidence_thermal_contrast",
                high_confidence_positions,
            )
        )
        handle.write("color red, high_confidence_thermal_contrast\n")
        handle.write("show spheres, high_confidence_thermal_contrast\n")
        handle.write("show sticks, high_confidence_thermal_contrast\n")
        handle.write(
            "label high_confidence_thermal_contrast and name CA, "
            "'%s%s' % (resn, resi)\n"
        )
        handle.write("\n")

        handle.write("orient protein_model\n")
        handle.write("zoom protein_model, 5\n")

    return pml_file


def main():
    if not SUMMARY_FILE.exists():
        raise SystemExit(
            f"Missing summary file: {SUMMARY_FILE}\n"
            "Run map_conservation.py first."
        )

    summary = pd.read_csv(SUMMARY_FILE)

    if summary.empty:
        raise SystemExit("The structure mapping summary is empty.")

    package_rows = []

    copy_mapping_tables()

    for _, row in summary.iterrows():
        dataset = row["dataset"]
        structure_file = Path(row["structure_file"])

        if not structure_file.exists():
            print(f"Structure file missing, skipped: {structure_file}")
            continue

        pdb_filename = f"{dataset}_{structure_file.name}"
        target_pdb = PDB_DIR / pdb_filename
        shutil.copy2(structure_file, target_pdb)

        conserved_file = (
            SOURCE_DIR
            / f"{dataset}_conserved_positions_mapped_to_structure.csv"
        )

        contrast_file = (
            SOURCE_DIR
            / f"{dataset}_thermal_contrast_positions_mapped_to_structure.csv"
        )

        conserved = (
            pd.read_csv(conserved_file)
            if conserved_file.exists()
            else pd.DataFrame()
        )

        contrast = (
            pd.read_csv(contrast_file)
            if contrast_file.exists()
            else pd.DataFrame()
        )

        pml_file = write_windows_friendly_pml(
            dataset=dataset,
            pdb_filename=pdb_filename,
            conserved=conserved,
            contrast=contrast,
        )

        package_rows.append(
            {
                "dataset": dataset,
                "reference_accession": row.get("reference_accession", ""),
                "reference_thermal_class": row.get("reference_thermal_class", ""),
                "pdb_file": str(target_pdb.relative_to(PACKAGE_DIR)),
                "pymol_script": str(pml_file.relative_to(PACKAGE_DIR)),
                "mapped_conserved_positions": row.get("mapped_conserved_positions", ""),
                "mapped_thermal_contrast_positions": row.get(
                    "mapped_thermal_contrast_positions",
                    "",
                ),
                "mapped_high_confidence_contrast_positions": row.get(
                    "mapped_high_confidence_contrast_positions",
                    "",
                ),
            }
        )

    package_summary = pd.DataFrame(package_rows)
    package_summary.to_csv(
        PACKAGE_DIR / "pymol_package_summary.csv",
        index=False,
    )

    readme = PACKAGE_DIR / "README.md"
    readme.write_text(
        """# PyMOL structure mapping package

This folder contains Windows-friendly PyMOL files for viewing conserved and thermal-contrast positions on representative xylanase structures.

## Folder contents

- `pdb/` contains the structure files.
- `pymol_scripts/` contains PyMOL `.pml` scripts.
- `tables/` contains mapped residue-position tables.
- `pymol_package_summary.csv` summarises the datasets and files.

## How to open on Windows

1. Download or clone the `thesis-evidence-pack` repository on Windows.
2. Open PyMOL.
3. Go to `File > Run Script`.
4. Select one `.pml` file from:

   `12_thermal_sequence_analysis/results/structure_mapping/pymol_scripts/`

Recommended first file:

   `fungal_gh11_thermal_annotated_3d_mapping_view.pml`

## Colour code

- Grey: full protein structure
- Cyan: 90% conserved positions
- Orange: thermal contrast positions
- Red: high-confidence thermal contrast positions

## Important interpretation note

The residue positions are mapped onto one representative structure per group. They should be interpreted as structure-mapped candidate positions, not as experimentally validated thermostability residues.
"""
    )

    print("PyMOL package created.")
    print(f"Package directory: {PACKAGE_DIR}")
    print("\nPackage summary:")
    print(package_summary.to_string(index=False))


if __name__ == "__main__":
    main()
