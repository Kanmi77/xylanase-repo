#!/usr/bin/env python3

from pathlib import Path
import argparse
import re

import pandas as pd
import yaml


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def split_pdb_ids(value):
    if pd.isna(value):
        return []

    text = str(value).strip()

    if not text:
        return []

    candidates = re.split(r"[;,|\s]+", text)
    pdb_ids = []

    for item in candidates:
        item = item.strip().upper()

        if re.fullmatch(r"[0-9][A-Z0-9]{3}", item):
            pdb_ids.append(item)

    return sorted(set(pdb_ids))


def main():
    parser = argparse.ArgumentParser(
        description="Build a structure inventory from the curated master table."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    master = pd.read_csv(args.input, low_memory=False)

    pdb_column = None
    for candidate in ["pdb", "xref_pdb", "pdb_ids"]:
        if candidate in master.columns:
            pdb_column = candidate
            break

    rows = []

    for _, row in master.iterrows():
        accession = row.get("uniprot_accession", "")
        organism_type = row.get("organism_type", "")
        gh_family = row.get("gh_family", "")
        organism = row.get("organism", "")

        pdb_ids = split_pdb_ids(row.get(pdb_column, "")) if pdb_column else []

        for pdb_id in pdb_ids:
            rows.append(
                {
                    "uniprot_accession": accession,
                    "organism_type": organism_type,
                    "gh_family": gh_family,
                    "organism": organism,
                    "structure_source": "pdb",
                    "pdb_id": pdb_id,
                    "structure_file": f"data/structures/pdb/{pdb_id}.pdb",
                }
            )

    inventory = pd.DataFrame(rows)

    if inventory.empty:
        inventory = pd.DataFrame(
            columns=[
                "uniprot_accession",
                "organism_type",
                "gh_family",
                "organism",
                "structure_source",
                "pdb_id",
                "structure_file",
            ]
        )

    inventory = inventory.drop_duplicates(
        subset=["uniprot_accession", "pdb_id"]
    ).sort_values(["organism_type", "gh_family", "uniprot_accession", "pdb_id"])

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(output, index=False)

    print(f"Structure inventory rows: {len(inventory)}")
    print(f"Unique PDB IDs: {inventory['pdb_id'].nunique() if len(inventory) else 0}")
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
