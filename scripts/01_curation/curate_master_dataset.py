#!/usr/bin/env python3

from pathlib import Path
import argparse
import re

import pandas as pd
import yaml


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def clean_column_name(name):
    return (
        name.strip()
        .lower()
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("[", "")
        .replace("]", "")
        .replace("-", "_")
    )


def standardise_columns(table):
    table = table.copy()
    table.columns = [clean_column_name(column) for column in table.columns]
    return table


def find_column(table, candidates):
    for candidate in candidates:
        if candidate in table.columns:
            return candidate
    return None


def detect_family(row, required_families):
    combined_text = " ".join(str(value) for value in row.values).upper()

    for family in required_families:
        family_upper = family.upper()
        family_number = re.sub(r"\D", "", family_upper)

        patterns = [
            family_upper,
            f"GH {family_number}",
            f"FAMILY {family_number}",
            f"GLYCOSIDE HYDROLASE {family_number}",
            f"GLYCOSIDE HYDROLASE FAMILY {family_number}",
        ]

        if any(pattern in combined_text for pattern in patterns):
            return family

    return ""


def detect_organism_group(lineage, organism_groups):
    lineage_text = str(lineage).lower()

    for group_name, group_config in organism_groups.items():
        for term in group_config.get("taxonomy_terms", []):
            if term.lower() in lineage_text:
                return group_name

    return "other"


def is_fragment(row):
    combined_text = " ".join(str(value) for value in row.values).lower()
    return "fragment" in combined_text


def main():
    parser = argparse.ArgumentParser(
        description="Curate raw UniProt records into the master enzyme dataset."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    curation = config["curation"]
    required_families = config["enzyme_search"]["required_families"]
    organism_groups = config["enzyme_search"]["organism_groups"]

    raw = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")
    raw = standardise_columns(raw)

    accession_col = find_column(raw, ["entry", "accession"])
    sequence_col = find_column(raw, ["sequence"])
    length_col = find_column(raw, ["length"])
    protein_col = find_column(raw, ["protein_names", "protein_name"])
    organism_col = find_column(raw, ["organism", "organism_name"])
    lineage_col = find_column(raw, ["lineage", "taxonomic_lineage"])

    if accession_col is None or sequence_col is None:
        raise ValueError("The UniProt table must contain accession and sequence columns.")

    curated = raw.copy()
    curated["uniprot_accession"] = curated[accession_col]
    curated["sequence"] = curated[sequence_col].str.replace(" ", "", regex=False)

    if length_col:
        curated["sequence_length"] = pd.to_numeric(curated[length_col], errors="coerce")
    else:
        curated["sequence_length"] = curated["sequence"].str.len()

    if protein_col:
        curated["protein_name"] = curated[protein_col]

    if organism_col:
        curated["organism"] = curated[organism_col]

    if lineage_col:
        curated["lineage"] = curated[lineage_col]
    else:
        curated["lineage"] = ""

    curated["gh_family"] = curated.apply(
        lambda row: detect_family(row, required_families),
        axis=1,
    )

    curated["organism_type"] = curated["lineage"].apply(
        lambda lineage: detect_organism_group(lineage, organism_groups)
    )

    curated = curated[curated["gh_family"].isin(required_families)]

    if not curation.get("keep_unmatched_organisms", False):
        configured_groups = set(organism_groups.keys())
        curated = curated[curated["organism_type"].isin(configured_groups)]

    if curation.get("require_sequence", True):
        curated = curated[curated["sequence"].str.len() > 0]

    if curation.get("remove_fragments", True):
        curated = curated[~curated.apply(is_fragment, axis=1)]

    curated = curated[
        curated["sequence_length"].between(
            curation["min_sequence_length"],
            curation["max_sequence_length"],
            inclusive="both",
        )
    ]

    curated = curated.drop_duplicates(subset=["uniprot_accession"])
    curated = curated.sort_values(["organism_type", "gh_family", "uniprot_accession"])

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    curated.to_csv(output_path, index=False)

    print(f"Curated records: {len(curated)}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
