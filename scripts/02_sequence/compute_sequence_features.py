#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import argparse

import pandas as pd


AA_MASS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "Q": 146.15, "E": 147.13, "G": 75.07, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}

HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

HELIX = set("AEHKLMQR")
TURN = set("NPGSD")
SHEET = set("VIYFWLT")


def clean_sequence(sequence):
    return "".join(aa for aa in str(sequence).upper() if aa in AA_MASS)


def fallback_features(sequence):
    sequence = clean_sequence(sequence)
    length = len(sequence)

    if length == 0:
        return {}

    counts = Counter(sequence)

    molecular_weight = sum(AA_MASS[aa] * count for aa, count in counts.items())
    gravy = sum(HYDROPATHY[aa] * count for aa, count in counts.items()) / length
    aromaticity = sum(counts[aa] for aa in "FWY") / length

    features = {
        "sequence_length": length,
        "molecular_weight": molecular_weight,
        "aromaticity": aromaticity,
        "instability_index": "",
        "isoelectric_point": "",
        "gravy": gravy,
        "helix_fraction": sum(counts[aa] for aa in HELIX) / length,
        "turn_fraction": sum(counts[aa] for aa in TURN) / length,
        "sheet_fraction": sum(counts[aa] for aa in SHEET) / length,
    }

    for aa in sorted(AA_MASS):
        features[f"aa_fraction_{aa}"] = counts[aa] / length

    return features


def biopython_features(sequence):
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError:
        return fallback_features(sequence)

    sequence = clean_sequence(sequence)

    if not sequence:
        return {}

    analysis = ProteinAnalysis(sequence)
    secondary = analysis.secondary_structure_fraction()
    features = {
        "sequence_length": len(sequence),
        "molecular_weight": analysis.molecular_weight(),
        "aromaticity": analysis.aromaticity(),
        "instability_index": analysis.instability_index(),
        "isoelectric_point": analysis.isoelectric_point(),
        "gravy": analysis.gravy(),
        "helix_fraction": secondary[0],
        "turn_fraction": secondary[1],
        "sheet_fraction": secondary[2],
    }

    counts = Counter(sequence)
    for aa in sorted(AA_MASS):
        features[f"aa_fraction_{aa}"] = counts[aa] / len(sequence)

    return features


def main():
    parser = argparse.ArgumentParser(
        description="Compute ProtParam-style sequence features."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    master = pd.read_csv(args.input, low_memory=False)
    rows = []

    for _, row in master.iterrows():
        features = biopython_features(row["sequence"])
        features["uniprot_accession"] = row["uniprot_accession"]
        features["gh_family"] = row.get("gh_family", "")
        features["organism_type"] = row.get("organism_type", "")
        rows.append(features)

    result = pd.DataFrame(rows)
    front_cols = ["uniprot_accession", "organism_type", "gh_family"]
    result = result[front_cols + [col for col in result.columns if col not in front_cols]]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    print(f"Sequence-feature rows: {len(result)}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
