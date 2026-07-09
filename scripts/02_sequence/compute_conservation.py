#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import argparse

import pandas as pd
import yaml


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_fasta(path):
    records = []
    header = None
    sequence_parts = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
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


def summarise_alignment(path, threshold):
    records = read_fasta(path)
    sequences = [sequence for _, sequence in records]

    if not sequences:
        raise ValueError(f"No sequences found in {path}")

    lengths = sorted(set(len(sequence) for sequence in sequences))

    if len(lengths) != 1:
        raise ValueError(f"Unequal alignment lengths in {path}: {lengths}")

    sequence_count = len(sequences)
    alignment_length = lengths[0]
    conserved_positions = 0

    for index in range(alignment_length):
        column = [sequence[index] for sequence in sequences]
        residues = [residue for residue in column if residue != "-"]

        if not residues:
            continue

        _, dominant_count = Counter(residues).most_common(1)[0]
        dominant_fraction = dominant_count / sequence_count

        if dominant_fraction >= threshold:
            conserved_positions += 1

    return {
        "group": path.stem.replace("_aligned", ""),
        "sequences": sequence_count,
        "alignment_length": alignment_length,
        "conserved_positions_90": conserved_positions,
        "conserved_percent": conserved_positions / alignment_length * 100,
        "conservation_threshold": threshold,
        "source_alignment": str(path),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compute conserved alignment positions."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    threshold = float(config["sequence_analysis"]["conservation_threshold"])

    input_dir = Path(args.input_dir)
    rows = []

    for alignment in sorted(input_dir.glob("*_aligned.fasta")):
        rows.append(summarise_alignment(alignment, threshold))

    if not rows:
        raise SystemExit(f"No aligned FASTA files found in {input_dir}")

    result = pd.DataFrame(rows)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    print(result.to_string(index=False))
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
