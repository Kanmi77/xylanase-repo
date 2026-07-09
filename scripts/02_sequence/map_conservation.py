#!/usr/bin/env python3
# Purpose: Map conserved positions to structures.

from pathlib import Path
import re

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"

ALIGNMENT_DIR = PROJECT_DIR / "results/thermal_sequence_analysis/alignments"
CONSERVATION_DIR = PROJECT_DIR / "results/thermal_sequence_analysis/conservation"
OUTPUT_DIR = PROJECT_DIR / "results/thermal_sequence_analysis/structure_mapping"

STRUCTURE_SEARCH_DIRS = [
    PROJECT_DIR / "data",
    PROJECT_DIR / "results",
]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATASETS_TO_MAP = [
    "bacterial_gh10_thermal_annotated",
    "bacterial_gh11_thermal_annotated",
    "fungal_gh10_thermal_annotated",
    "fungal_gh11_thermal_annotated",
]

MAX_CONSERVED_TO_SHOW = 150


def read_fasta(path):
    records = []
    header = None
    sequence_parts = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_parts)))

                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if header is not None:
        records.append((header, "".join(sequence_parts)))

    return records


def parse_header(header):
    details = {
        "accession": "",
        "organism_type": "unknown",
        "gh_family": "unknown",
        "thermal_class": "unknown",
    }

    parts = str(header).split("|")

    if len(parts) > 0:
        details["accession"] = parts[0]

    if len(parts) > 1:
        details["organism_type"] = parts[1]

    if len(parts) > 2:
        details["gh_family"] = parts[2]

    for part in parts:
        if part.startswith("thermal="):
            details["thermal_class"] = part.replace("thermal=", "")

    return details


def find_structure_file(accession):
    patterns = [
        f"{accession}.pdb",
        f"{accession}_*.pdb",
        f"*{accession}*.pdb",
    ]

    for root in STRUCTURE_SEARCH_DIRS:
        if not root.exists():
            continue

        for pattern in patterns:
            matches = sorted(root.rglob(pattern))

            if matches:
                return matches[0]

    return None


def choose_reference_with_structure(records):
    labelled_priority = ["thermophilic", "mesophilic_lower", "unknown"]

    for thermal_class in labelled_priority:
        for header, sequence in records:
            details = parse_header(header)

            if details["thermal_class"] != thermal_class:
                continue

            structure_file = find_structure_file(details["accession"])

            if structure_file is not None:
                return header, sequence, structure_file

    return None, None, None


def build_alignment_to_reference_map(reference_sequence):
    mapping = {}
    residue_position = 0

    for alignment_position, residue in enumerate(reference_sequence, start=1):
        if residue == "-":
            mapping[alignment_position] = {
                "reference_residue_position": pd.NA,
                "reference_residue": "-",
            }
        else:
            residue_position += 1
            mapping[alignment_position] = {
                "reference_residue_position": residue_position,
                "reference_residue": residue,
            }

    return mapping


def clean_residue_positions(values):
    positions = []

    for value in values:
        if pd.isna(value):
            continue

        try:
            position = int(value)
        except ValueError:
            continue

        if position > 0:
            positions.append(position)

    return sorted(set(positions))


def get_column(table, possible_names):
    for name in possible_names:
        if name in table.columns:
            return name

    raise KeyError(f"None of these columns were found: {possible_names}")


def map_conserved_positions(dataset_name, reference_map):
    conserved_file = CONSERVATION_DIR / f"{dataset_name}_conserved_positions_90pct.csv"

    if not conserved_file.exists():
        return pd.DataFrame()

    conserved = pd.read_csv(conserved_file)

    major_col = get_column(
        conserved,
        ["overall_major_residue", "overall_major_aa"],
    )

    rows = []

    for _, row in conserved.iterrows():
        alignment_position = int(row["alignment_position"])
        reference_info = reference_map.get(alignment_position)

        if reference_info is None:
            continue

        rows.append(
            {
                "alignment_position": alignment_position,
                "conserved_residue": row[major_col],
                "conservation_fraction": row["overall_major_fraction"],
                "non_gap_count": row["non_gap_count"],
                "reference_residue_position": reference_info["reference_residue_position"],
                "reference_residue": reference_info["reference_residue"],
            }
        )

    mapped = pd.DataFrame(rows)
    mapped = mapped.dropna(subset=["reference_residue_position"]).copy()

    if not mapped.empty:
        mapped["reference_residue_position"] = mapped[
            "reference_residue_position"
        ].astype(int)

    return mapped


def map_thermal_contrast_positions(dataset_name, reference_map):
    contrast_file = CONSERVATION_DIR / f"{dataset_name}_thermal_contrast_positions.csv"

    if not contrast_file.exists():
        return pd.DataFrame()

    contrast = pd.read_csv(contrast_file)

    if contrast.empty:
        return pd.DataFrame()

    therm_col = get_column(
        contrast,
        ["thermophilic_major_residue", "thermophilic_major_aa"],
    )

    meso_col = get_column(
        contrast,
        ["mesophilic_lower_major_residue", "mesophilic_lower_major_aa"],
    )

    rows = []

    for _, row in contrast.iterrows():
        alignment_position = int(row["alignment_position"])
        reference_info = reference_map.get(alignment_position)

        if reference_info is None:
            continue

        therm_fraction = float(row["thermophilic_major_fraction"])
        meso_fraction = float(row["mesophilic_lower_major_fraction"])

        high_confidence = (
            row["thermophilic_non_gap"] >= 5
            and row["mesophilic_lower_non_gap"] >= 5
            and therm_fraction >= 0.70
            and meso_fraction >= 0.70
        )

        rows.append(
            {
                "alignment_position": alignment_position,
                "thermophilic_major_residue": row[therm_col],
                "thermophilic_major_fraction": therm_fraction,
                "thermophilic_non_gap": row["thermophilic_non_gap"],
                "mesophilic_lower_major_residue": row[meso_col],
                "mesophilic_lower_major_fraction": meso_fraction,
                "mesophilic_lower_non_gap": row["mesophilic_lower_non_gap"],
                "contrast_strength": (therm_fraction + meso_fraction) / 2,
                "high_confidence": high_confidence,
                "reference_residue_position": reference_info["reference_residue_position"],
                "reference_residue": reference_info["reference_residue"],
            }
        )

    mapped = pd.DataFrame(rows)
    mapped = mapped.dropna(subset=["reference_residue_position"]).copy()

    if not mapped.empty:
        mapped["reference_residue_position"] = mapped[
            "reference_residue_position"
        ].astype(int)

        mapped = mapped.sort_values(
            ["high_confidence", "contrast_strength"],
            ascending=[False, False],
        )

    return mapped


def pymol_residue_selection(name, positions):
    if not positions:
        return f"select {name}, none\n"

    chunks = []
    sorted_positions = sorted(set(positions))

    for position in sorted_positions:
        chunks.append(str(position))

    residue_text = "+".join(chunks)

    return f"select {name}, protein_model and resi {residue_text}\n"


def write_pymol_script(
    dataset_name,
    structure_file,
    reference_details,
    conserved_positions,
    contrast_positions,
    high_confidence_positions,
):
    pml_file = OUTPUT_DIR / f"{dataset_name}_conserved_positions_3d_view.pml"
    png_file = OUTPUT_DIR / f"{dataset_name}_conserved_positions_3d_view.png"

    with pml_file.open("w") as handle:
        handle.write(f"load {structure_file}, protein_model\n")
        handle.write("hide everything\n")
        handle.write("show cartoon, protein_model\n")
        handle.write("color gray80, protein_model\n")
        handle.write("set cartoon_transparency, 0.15, protein_model\n")
        handle.write("bg_color white\n")
        handle.write("set ray_opaque_background, off\n")
        handle.write("set antialias, 2\n")
        handle.write("set sphere_scale, 0.35\n")
        handle.write("set stick_radius, 0.18\n")
        handle.write("set label_size, 16\n")
        handle.write("set label_color, black\n")
        handle.write("\n")

        handle.write(pymol_residue_selection("conserved_90", conserved_positions))
        handle.write("color cyan, conserved_90\n")
        handle.write("show sticks, conserved_90\n")
        handle.write("\n")

        handle.write(pymol_residue_selection("thermal_contrast", contrast_positions))
        handle.write("color orange, thermal_contrast\n")
        handle.write("show spheres, thermal_contrast\n")
        handle.write("show sticks, thermal_contrast\n")
        handle.write("\n")

        handle.write(
            pymol_residue_selection(
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
        handle.write("ray 2400, 1800\n")
        handle.write(f"png {png_file}, dpi=600\n")

    return pml_file, png_file


def process_dataset(dataset_name):
    alignment_file = ALIGNMENT_DIR / f"{dataset_name}_aligned.fasta"

    if not alignment_file.exists():
        print(f"Missing alignment: {alignment_file}")
        return None

    records = read_fasta(alignment_file)

    if not records:
        print(f"No records in alignment: {alignment_file}")
        return None

    reference_header, reference_sequence, structure_file = choose_reference_with_structure(records)

    if structure_file is None:
        print(f"No structure file found for any sequence in {dataset_name}")
        return None

    reference_details = parse_header(reference_header)
    reference_map = build_alignment_to_reference_map(reference_sequence)

    conserved = map_conserved_positions(dataset_name, reference_map)
    contrast = map_thermal_contrast_positions(dataset_name, reference_map)

    conserved_output = OUTPUT_DIR / f"{dataset_name}_conserved_positions_mapped_to_structure.csv"
    contrast_output = OUTPUT_DIR / f"{dataset_name}_thermal_contrast_positions_mapped_to_structure.csv"

    conserved.to_csv(conserved_output, index=False)
    contrast.to_csv(contrast_output, index=False)

    conserved_positions = clean_residue_positions(
        conserved.sort_values(
            "conservation_fraction",
            ascending=False,
        )["reference_residue_position"].head(MAX_CONSERVED_TO_SHOW)
    )

    contrast_positions = clean_residue_positions(
        contrast["reference_residue_position"]
        if not contrast.empty
        else []
    )

    high_confidence_positions = clean_residue_positions(
        contrast.loc[
            contrast["high_confidence"] == True,
            "reference_residue_position",
        ]
        if not contrast.empty
        else []
    )

    pml_file, png_file = write_pymol_script(
        dataset_name=dataset_name,
        structure_file=structure_file,
        reference_details=reference_details,
        conserved_positions=conserved_positions,
        contrast_positions=contrast_positions,
        high_confidence_positions=high_confidence_positions,
    )

    return {
        "dataset": dataset_name,
        "reference_accession": reference_details["accession"],
        "reference_thermal_class": reference_details["thermal_class"],
        "structure_file": str(structure_file),
        "mapped_conserved_positions": len(conserved),
        "mapped_thermal_contrast_positions": len(contrast),
        "mapped_high_confidence_contrast_positions": len(high_confidence_positions),
        "pymol_script": str(pml_file),
        "pymol_png_output": str(png_file),
    }


def main():
    summaries = []

    for dataset_name in DATASETS_TO_MAP:
        print(f"Mapping dataset: {dataset_name}")
        summary = process_dataset(dataset_name)

        if summary is not None:
            summaries.append(summary)
            print(f"  Reference: {summary['reference_accession']}")
            print(f"  Structure: {summary['structure_file']}")
            print(f"  PyMOL script: {summary['pymol_script']}")
        print()

    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(
        OUTPUT_DIR / "structure_mapping_summary.csv",
        index=False,
    )

    print("Structure mapping summary:")
    if summary_table.empty:
        print("No mappings were created. Check that structure PDB files exist.")
    else:
        print(summary_table.to_string(index=False))

    print(f"\nOutput directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
