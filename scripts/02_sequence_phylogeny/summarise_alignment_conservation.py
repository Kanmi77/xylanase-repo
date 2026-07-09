#!/usr/bin/env python3

from pathlib import Path
from collections import Counter

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"

ALIGNMENT_DIR = PROJECT_DIR / "results/thermal_sequence_analysis/alignments"
CONSERVATION_DIR = PROJECT_DIR / "results/thermal_sequence_analysis/conservation"

CONSERVATION_DIR.mkdir(parents=True, exist_ok=True)

OVERALL_CONSERVATION_CUTOFF = 0.90
MIN_CLASS_COUNT = 5


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
    parts = header.split("|")

    details = {
        "header": header,
        "accession": parts[0] if len(parts) > 0 else "unknown",
        "organism_type": parts[1] if len(parts) > 1 else "unknown",
        "gh_family": parts[2] if len(parts) > 2 else "unknown",
        "thermal_class": "unknown",
    }

    for part in parts:
        if part.startswith("thermal="):
            details["thermal_class"] = part.replace("thermal=", "")

    return details


def majority_residue(residues):
    residues = [residue for residue in residues if residue != "-"]

    if not residues:
        return "-", 0, 0.0

    counts = Counter(residues)
    residue, count = counts.most_common(1)[0]

    return residue, count, count / len(residues)


def analyse_alignment(alignment_file):
    records = read_fasta(alignment_file)

    if not records:
        return None

    metadata = [parse_header(header) for header, _ in records]
    sequences = [sequence for _, sequence in records]

    alignment_length = len(sequences[0])
    dataset_name = alignment_file.stem.replace("_aligned", "")

    position_rows = []
    conserved_rows = []
    thermal_rows = []

    for position_index in range(alignment_length):
        column = [sequence[position_index] for sequence in sequences]

        overall_residue, overall_count, overall_fraction = majority_residue(column)
        non_gap_count = sum(1 for residue in column if residue != "-")

        position_row = {
            "alignment_position": position_index + 1,
            "non_gap_count": non_gap_count,
            "overall_major_residue": overall_residue,
            "overall_major_count": overall_count,
            "overall_major_fraction": overall_fraction,
        }

        position_rows.append(position_row)

        if (
            overall_residue != "-"
            and overall_fraction >= OVERALL_CONSERVATION_CUTOFF
        ):
            conserved_rows.append(position_row)

        class_summary = {}

        for thermal_class in ["thermophilic", "mesophilic_lower", "unknown"]:
            class_column = [
                sequences[index][position_index]
                for index, item in enumerate(metadata)
                if item["thermal_class"] == thermal_class
            ]

            class_residue, class_count, class_fraction = majority_residue(class_column)
            class_non_gap = sum(1 for residue in class_column if residue != "-")

            class_summary[thermal_class] = {
                "residue": class_residue,
                "count": class_count,
                "fraction": class_fraction,
                "non_gap": class_non_gap,
            }

        thermal_rows.append(
            {
                "alignment_position": position_index + 1,
                "thermophilic_non_gap": class_summary["thermophilic"]["non_gap"],
                "thermophilic_major_residue": class_summary["thermophilic"]["residue"],
                "thermophilic_major_fraction": class_summary["thermophilic"]["fraction"],
                "mesophilic_lower_non_gap": class_summary["mesophilic_lower"]["non_gap"],
                "mesophilic_lower_major_residue": class_summary["mesophilic_lower"]["residue"],
                "mesophilic_lower_major_fraction": class_summary["mesophilic_lower"]["fraction"],
                "unknown_non_gap": class_summary["unknown"]["non_gap"],
                "unknown_major_residue": class_summary["unknown"]["residue"],
                "unknown_major_fraction": class_summary["unknown"]["fraction"],
                "thermal_major_residue_differs": (
                    class_summary["thermophilic"]["residue"]
                    != class_summary["mesophilic_lower"]["residue"]
                ),
            }
        )

    position_table = pd.DataFrame(position_rows)
    conserved_table = pd.DataFrame(conserved_rows)
    thermal_table = pd.DataFrame(thermal_rows)

    thermal_contrast = thermal_table[
        (thermal_table["thermophilic_non_gap"] >= MIN_CLASS_COUNT)
        & (thermal_table["mesophilic_lower_non_gap"] >= MIN_CLASS_COUNT)
        & (thermal_table["thermal_major_residue_differs"])
    ].copy()

    thermal_counts = Counter(item["thermal_class"] for item in metadata)

    position_table.to_csv(
        CONSERVATION_DIR / f"{dataset_name}_position_conservation.csv",
        index=False,
    )

    conserved_table.to_csv(
        CONSERVATION_DIR / f"{dataset_name}_conserved_positions_90pct.csv",
        index=False,
    )

    thermal_table.to_csv(
        CONSERVATION_DIR / f"{dataset_name}_thermal_class_conservation.csv",
        index=False,
    )

    thermal_contrast.to_csv(
        CONSERVATION_DIR / f"{dataset_name}_thermal_contrast_positions.csv",
        index=False,
    )

    return {
        "dataset": dataset_name,
        "sequence_count": len(records),
        "alignment_length": alignment_length,
        "thermophilic_count": thermal_counts.get("thermophilic", 0),
        "mesophilic_lower_count": thermal_counts.get("mesophilic_lower", 0),
        "unknown_count": thermal_counts.get("unknown", 0),
        "conserved_positions_90pct": len(conserved_table),
        "thermal_contrast_positions": len(thermal_contrast),
    }


def main():
    summaries = []

    for alignment_file in sorted(ALIGNMENT_DIR.glob("*_aligned.fasta")):
        print(f"Analysing: {alignment_file.name}")
        summary = analyse_alignment(alignment_file)

        if summary is not None:
            summaries.append(summary)

    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(
        CONSERVATION_DIR / "conservation_summary.csv",
        index=False,
    )

    print("\nConservation summary:")
    print(summary_table.to_string(index=False))

    print(f"\nOutput directory: {CONSERVATION_DIR}")


if __name__ == "__main__":
    main()
