#!/usr/bin/env python3

from pathlib import Path
import argparse
import re

import pandas as pd
import yaml


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def safe_name(value):
    value = str(value).strip().lower()
    value = re.sub(r"[^a-zA-Z0-9]+", "_", value)
    return value.strip("_")


def wrap_sequence(sequence, width=80):
    sequence = str(sequence)
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def main():
    parser = argparse.ArgumentParser(
        description="Split curated sequences into group FASTA files."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    families = config["enzyme_search"]["required_families"]
    organism_groups = list(config["enzyme_search"]["organism_groups"].keys())

    master = pd.read_csv(args.input, low_memory=False)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    written = 0

    for organism_group in organism_groups:
        for family in families:
            subset = master[
                (master["organism_type"] == organism_group)
                & (master["gh_family"] == family)
            ]

            if subset.empty:
                continue

            fasta_path = output_dir / f"{safe_name(organism_group)}_{family}.fasta"

            with open(fasta_path, "w", encoding="utf-8") as handle:
                for _, row in subset.iterrows():
                    accession = row["uniprot_accession"]
                    sequence = row["sequence"]
                    organism = row.get("organism", "")
                    handle.write(f">{accession}|{organism_group}|{family}|{organism}\n")
                    handle.write(wrap_sequence(sequence) + "\n")

            written += 1
            print(f"Wrote {fasta_path}: {len(subset)} sequences")

    if written == 0:
        raise SystemExit("No group FASTA files were written.")


if __name__ == "__main__":
    main()
