#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"
OUTPUT_DIR = PROJECT_DIR / "results/structural_results/separate_structural_features"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


MASTER_FILE = (
    PROJECT_DIR
    / "data/curated/xylanase_master_gh10_gh11_frozen_thesis_with_refseq_experimental_metadata.csv"
)

PDB_FEATURE_FILE = (
    PROJECT_DIR
    / "data/curated/xylanase_structured_subset_with_foldx_norm.csv"
)

MODELLER_FEATURE_FILE = (
    PROJECT_DIR
    / "results/foldx/modeller/structural_features_with_length.csv"
)


ACCESSION_COLUMNS = [
    "uniprot_accession",
    "accession",
    "Entry",
    "entry",
    "query_accession",
    "target_accession",
    "protein_accession",
]


GROUP_COLUMNS = [
    "organism_type",
    "gh_family",
]


RAW_FEATURE_COLUMNS = [
    "chain_length",
    "hbond_proxy_count",
    "hbond_count",
    "salt_bridge_count",
    "disulfide_count",
    "sasa_total",
    "total_sasa",
    "foldx_wt_total_energy",
    "foldx_energy_per_residue",
]


NORMALISED_FEATURE_COLUMNS = [
    "hbond_proxy_per_residue",
    "hbond_per_residue",
    "salt_bridge_per_residue",
    "disulfide_per_residue",
    "sasa_per_residue",
    "foldx_energy_per_residue",
]


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    table = pd.read_csv(path)
    table.columns = [str(column).strip() for column in table.columns]

    return table


def find_column(table, possible_names):
    lower_map = {column.lower(): column for column in table.columns}

    for name in possible_names:
        if name.lower() in lower_map:
            return lower_map[name.lower()]

    return None


def infer_accession_from_text(value):
    if pd.isna(value):
        return pd.NA

    text = str(value)

    patterns = [
        r"\b[A-NR-Z][0-9][A-Z0-9]{3}[0-9]\b",
        r"\b[A-Z0-9]{6,10}\b",
    ]

    for pattern in patterns:
        match = re.search(pattern, text)

        if match:
            return match.group(0)

    return pd.NA


def add_accession_column(table, table_name):
    table = table.copy()

    accession_column = find_column(table, ACCESSION_COLUMNS)

    if accession_column is not None:
        table["accession_for_merge"] = table[accession_column].astype(str)
        return table

    text_columns = [
        column
        for column in table.columns
        if table[column].dtype == "object"
    ]

    for column in text_columns:
        inferred = table[column].apply(infer_accession_from_text)

        if inferred.notna().sum() > 0:
            table["accession_for_merge"] = inferred
            print(
                f"{table_name}: accession inferred from column '{column}' "
                f"for {inferred.notna().sum()} rows."
            )
            return table

    table["accession_for_merge"] = pd.NA
    print(f"{table_name}: no accession column could be identified.")

    return table


def prepare_master_metadata(master):
    master = add_accession_column(master, "master")

    metadata_columns = [
        "accession_for_merge",
        "organism_type",
        "gh_family",
        "organism",
        "protein_name",
    ]

    metadata_columns = [
        column for column in metadata_columns
        if column in master.columns
    ]

    metadata = master[metadata_columns].drop_duplicates("accession_for_merge")

    return metadata


def merge_metadata(feature_table, metadata):
    feature_table = feature_table.copy()

    if "organism_type" in feature_table.columns and "gh_family" in feature_table.columns:
        return feature_table

    merged = feature_table.merge(
        metadata,
        on="accession_for_merge",
        how="left",
        suffixes=("", "_master"),
    )

    return merged


def add_normalised_features(table):
    table = table.copy()

    length_column = find_column(table, ["chain_length", "sequence_length", "length"])

    if length_column is None:
        return table

    table[length_column] = pd.to_numeric(table[length_column], errors="coerce")

    feature_map = {
        "hbond_proxy_count": "hbond_proxy_per_residue",
        "hbond_count": "hbond_per_residue",
        "salt_bridge_count": "salt_bridge_per_residue",
        "disulfide_count": "disulfide_per_residue",
        "sasa_total": "sasa_per_residue",
        "total_sasa": "sasa_per_residue",
    }

    for raw_column, normalised_column in feature_map.items():
        if raw_column not in table.columns:
            continue

        table[raw_column] = pd.to_numeric(table[raw_column], errors="coerce")
        table[normalised_column] = table[raw_column] / table[length_column]

    if (
        "foldx_wt_total_energy" in table.columns
        and "foldx_energy_per_residue" not in table.columns
    ):
        table["foldx_wt_total_energy"] = pd.to_numeric(
            table["foldx_wt_total_energy"],
            errors="coerce",
        )

        table["foldx_energy_per_residue"] = (
            table["foldx_wt_total_energy"] / table[length_column]
        )

    return table


def numeric_summary(table, feature_columns):
    available = [
        column for column in feature_columns
        if column in table.columns
    ]

    if not available:
        return pd.DataFrame()

    working = table.copy()

    for column in available:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    summary = (
        working[available]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .T
        .reset_index()
        .rename(columns={"index": "feature"})
    )

    return summary


def grouped_summary(table, group_columns, feature_columns):
    available_groups = [
        column for column in group_columns
        if column in table.columns
    ]

    available_features = [
        column for column in feature_columns
        if column in table.columns
    ]

    if not available_groups or not available_features:
        return pd.DataFrame()

    working = table.copy()

    for column in available_features:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    summary = (
        working.groupby(available_groups, dropna=False)[available_features]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )

    summary.columns = [
        "_".join([str(part) for part in column if str(part) != ""]).strip("_")
        for column in summary.columns
    ]

    return summary


def count_table(table, group_columns):
    available_groups = [
        column for column in group_columns
        if column in table.columns
    ]

    if not available_groups:
        return pd.DataFrame()

    return (
        table.groupby(available_groups, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(available_groups)
    )


def save_branch_results(name, table):
    branch_dir = OUTPUT_DIR / name
    branch_dir.mkdir(parents=True, exist_ok=True)

    raw_features = [
        column for column in RAW_FEATURE_COLUMNS
        if column in table.columns
    ]

    normalised_features = [
        column for column in NORMALISED_FEATURE_COLUMNS
        if column in table.columns
    ]

    all_features = list(dict.fromkeys(raw_features + normalised_features))

    table.to_csv(branch_dir / f"{name}_structural_features_with_metadata.csv", index=False)

    count_table(table, GROUP_COLUMNS).to_csv(
        branch_dir / f"{name}_counts_by_organism_and_family.csv",
        index=False,
    )

    numeric_summary(table, all_features).to_csv(
        branch_dir / f"{name}_overall_feature_summary.csv",
        index=False,
    )

    grouped_summary(table, GROUP_COLUMNS, all_features).to_csv(
        branch_dir / f"{name}_feature_summary_by_organism_and_family.csv",
        index=False,
    )

    grouped_summary(table, ["gh_family"], all_features).to_csv(
        branch_dir / f"{name}_feature_summary_by_gh_family.csv",
        index=False,
    )

    grouped_summary(table, ["organism_type"], all_features).to_csv(
        branch_dir / f"{name}_feature_summary_by_organism_type.csv",
        index=False,
    )

    return {
        "branch": name,
        "rows": len(table),
        "unique_accessions": table["accession_for_merge"].dropna().nunique()
        if "accession_for_merge" in table.columns
        else pd.NA,
        "features_summarised": ", ".join(all_features),
    }


def main():
    master = read_csv(MASTER_FILE)
    metadata = prepare_master_metadata(master)

    pdb = read_csv(PDB_FEATURE_FILE)
    modeller = read_csv(MODELLER_FEATURE_FILE)

    pdb = add_accession_column(pdb, "pdb")
    modeller = add_accession_column(modeller, "modeller")

    pdb = merge_metadata(pdb, metadata)
    modeller = merge_metadata(modeller, metadata)

    pdb = add_normalised_features(pdb)
    modeller = add_normalised_features(modeller)

    branch_summaries = [
        save_branch_results("pdb", pdb),
        save_branch_results("modeller", modeller),
    ]

    branch_summary = pd.DataFrame(branch_summaries)
    branch_summary.to_csv(OUTPUT_DIR / "separate_structural_feature_branch_summary.csv", index=False)

    readme = OUTPUT_DIR / "README.md"
    readme.write_text(
        """# Separate structural-feature analysis

This folder reports structural features separately for the PDB-supported branch and the MODELLER-derived branch.

## Why separate?

PDB-linked structures and MODELLER-derived structures do not have the same evidence level. Therefore, they are summarised separately and should not be presented as one uniform downstream structural dataset.

## Structural features

The extracted features include chain length, hydrogen-bond proxy count, salt-bridge count, disulfide count, total solvent-accessible surface area, FoldX total energy and FoldX energy per residue where available.

Normalised per-residue features are calculated where the required raw feature and chain length are available.

## Output folders

- `pdb/` contains PDB-derived feature summaries.
- `modeller/` contains MODELLER-derived feature summaries.
"""
    )

    print("\nSeparate structural-feature branch summary:")
    print(branch_summary.to_string(index=False))

    print(f"\nOutput directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
