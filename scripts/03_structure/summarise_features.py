#!/usr/bin/env python3
# Purpose: Summarise structural features.

from pathlib import Path
import re

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"

OUTPUT_DIR = (
    PROJECT_DIR
    / "results/structural_results/structural_features_without_foldx"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


MASTER_FILE = (
    PROJECT_DIR
    / "data/curated/master_metadata.csv"
)

PDB_FEATURE_FILE = (
    PROJECT_DIR
    / "data/curated/xylanase_structured_subset_with_foldx_norm.csv"
)

MODELLER_CANDIDATE_FILES = [
    PROJECT_DIR
    / "results/structure_features/combined_structural_features.csv",
    PROJECT_DIR
    / "results/structural_features_foldx_summary/structural_features_with_calculated_normalised_values.csv",
    PROJECT_DIR
    / "results/foldx/modeller/structural_features_full.csv",
    PROJECT_DIR
    / "results/foldx/modeller/structural_features_with_length.csv",
]


ACCESSION_COLUMNS = [
    "uniprot_accession",
    "accession",
    "entry",
    "Entry",
    "query_accession",
    "target_accession",
    "protein_accession",
    "model_accession",
]

GROUP_COLUMNS = [
    "organism_type",
    "gh_family",
]

FEATURE_ALIASES = {
    "chain_length": [
        "chain_length",
        "sequence_length",
        "length",
        "residue_count",
    ],
    "hbond_proxy_count": [
        "hbond_proxy_count",
        "hbond_count",
        "hydrogen_bond_count",
        "hydrogen_bonds",
        "h_bond_count",
    ],
    "salt_bridge_count": [
        "salt_bridge_count",
        "salt_bridges",
        "salt_bridge_proxy_count",
    ],
    "disulfide_count": [
        "disulfide_count",
        "disulfide_bond_count",
        "disulfide_bonds",
    ],
    "sasa_total": [
        "sasa_total",
        "total_sasa",
        "sasa",
        "solvent_accessible_surface_area",
    ],
}

FOLDX_KEYWORDS = [
    "foldx",
    "energy",
    "stability",
    "ddg",
]


def read_csv(path):
    table = pd.read_csv(path)
    table.columns = [str(column).strip() for column in table.columns]
    return table


def find_column(table, candidates):
    lower_map = {column.lower(): column for column in table.columns}

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    return None


def infer_accession_from_text(value):
    if pd.isna(value):
        return pd.NA

    text = str(value)

    patterns = [
        r"\b[A-NR-Z][0-9][A-Z0-9]{3}[0-9]\b",
        r"\b[A-Z][A-Z0-9]{5,9}\b",
    ]

    for pattern in patterns:
        match = re.search(pattern, text)

        if match:
            return match.group(0)

    return pd.NA


def add_accession_column(table, table_name):
    table = table.copy()

    accession_column = find_column(table, ACCESSION_COLUMNS)

    if accession_column:
        table["accession_for_merge"] = table[accession_column].astype(str)
        return table

    object_columns = [
        column for column in table.columns
        if table[column].dtype == "object"
    ]

    best_column = None
    best_count = 0
    best_values = None

    for column in object_columns:
        inferred = table[column].apply(infer_accession_from_text)
        count = inferred.notna().sum()

        if count > best_count:
            best_column = column
            best_count = count
            best_values = inferred

    if best_column and best_count > 0:
        table["accession_for_merge"] = best_values
        print(
            f"{table_name}: accession inferred from '{best_column}' "
            f"for {best_count} rows."
        )
    else:
        table["accession_for_merge"] = pd.NA
        print(f"{table_name}: no accession column could be identified.")

    return table


def prepare_master_metadata():
    master = read_csv(MASTER_FILE)
    master = add_accession_column(master, "master")

    wanted_columns = [
        "accession_for_merge",
        "organism_type",
        "gh_family",
        "organism",
        "protein_name",
    ]

    wanted_columns = [
        column for column in wanted_columns
        if column in master.columns
    ]

    metadata = master[wanted_columns].drop_duplicates("accession_for_merge")

    return metadata


def merge_metadata(table, metadata):
    table = table.copy()

    if "accession_for_merge" not in table.columns:
        table = add_accession_column(table, "feature_table")

    table = table.merge(
        metadata,
        on="accession_for_merge",
        how="left",
        suffixes=("", "_master"),
    )

    for column in ["organism_type", "gh_family"]:
        master_column = f"{column}_master"

        if master_column in table.columns:
            if column in table.columns:
                table[column] = table[column].fillna(table[master_column])
            else:
                table[column] = table[master_column]

    return table


def standardise_feature_columns(table):
    table = table.copy()

    for standard_name, aliases in FEATURE_ALIASES.items():
        source_column = find_column(table, aliases)

        if source_column is not None:
            table[standard_name] = pd.to_numeric(
                table[source_column],
                errors="coerce",
            )

    return table


def remove_foldx_columns(table):
    keep_columns = []

    for column in table.columns:
        lower = column.lower()

        if any(keyword in lower for keyword in FOLDX_KEYWORDS):
            continue

        keep_columns.append(column)

    return table[keep_columns].copy()


def add_per_residue_features(table):
    table = table.copy()

    if "chain_length" not in table.columns:
        return table

    table["chain_length"] = pd.to_numeric(
        table["chain_length"],
        errors="coerce",
    )

    normalisation_map = {
        "hbond_proxy_count": "hbond_proxy_per_residue",
        "salt_bridge_count": "salt_bridge_per_residue",
        "disulfide_count": "disulfide_per_residue",
        "sasa_total": "sasa_per_residue",
    }

    for raw_column, normalised_column in normalisation_map.items():
        if raw_column not in table.columns:
            continue

        table[raw_column] = pd.to_numeric(table[raw_column], errors="coerce")

        table[normalised_column] = table[raw_column] / table["chain_length"]

    return table


def valid_feature_rows(table):
    table = table.copy()

    if "chain_length" in table.columns:
        table = table[
            pd.to_numeric(table["chain_length"], errors="coerce").fillna(0) > 0
        ].copy()

    return table


def feature_columns_available(table):
    possible_features = [
        "chain_length",
        "hbond_proxy_count",
        "salt_bridge_count",
        "disulfide_count",
        "sasa_total",
        "hbond_proxy_per_residue",
        "salt_bridge_per_residue",
        "disulfide_per_residue",
        "sasa_per_residue",
    ]

    return [
        column for column in possible_features
        if column in table.columns
    ]


def score_modeller_candidate(path):
    if not path.exists():
        return None

    table = read_csv(path)

    source_column = find_column(
        table,
        ["structure_source", "source", "model_source"],
    )

    if source_column:
        modeller_rows = table[
            table[source_column].astype(str).str.lower().str.contains(
                "modeller",
                na=False,
            )
        ].copy()

        if not modeller_rows.empty:
            table = modeller_rows

    table = standardise_feature_columns(table)
    table = remove_foldx_columns(table)

    available_features = feature_columns_available(table)

    non_length_features = [
        column for column in available_features
        if column != "chain_length"
    ]

    score = len(non_length_features) * 10 + len(available_features)

    return {
        "path": path,
        "table": table,
        "rows": len(table),
        "available_features": available_features,
        "score": score,
    }


def load_best_modeller_table():
    candidates = []

    for path in MODELLER_CANDIDATE_FILES:
        candidate = score_modeller_candidate(path)

        if candidate:
            candidates.append(candidate)

    if not candidates:
        raise FileNotFoundError("No MODELLER structural-feature file found.")

    candidates = sorted(
        candidates,
        key=lambda item: (item["score"], item["rows"]),
        reverse=True,
    )

    best = candidates[0]

    source_summary = pd.DataFrame(
        [
            {
                "candidate_file": str(item["path"]),
                "rows": item["rows"],
                "score": item["score"],
                "features_found": ", ".join(item["available_features"]),
            }
            for item in candidates
        ]
    )

    source_summary.to_csv(
        OUTPUT_DIR / "modeller_feature_source_summary.csv",
        index=False,
    )

    print("Selected MODELLER feature source:")
    print(best["path"])
    print("Features found:", ", ".join(best["available_features"]))

    return best["table"], best["path"]


def numeric_summary(table, feature_columns):
    if not feature_columns:
        return pd.DataFrame()

    working = table.copy()

    for column in feature_columns:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    summary = (
        working[feature_columns]
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

    if not available_groups or not feature_columns:
        return pd.DataFrame()

    working = table.copy()

    for column in feature_columns:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    summary = (
        working.groupby(available_groups, dropna=False)[feature_columns]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )

    summary.columns = [
        "_".join([str(part) for part in column if str(part) != ""]).strip("_")
        for column in summary.columns
    ]

    return summary


def count_by_group(table, group_columns):
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


def save_branch(name, table, source_file):
    branch_dir = OUTPUT_DIR / name
    branch_dir.mkdir(parents=True, exist_ok=True)

    table = remove_foldx_columns(table)
    table = standardise_feature_columns(table)
    table = add_per_residue_features(table)
    table = valid_feature_rows(table)

    features = feature_columns_available(table)

    table.to_csv(
        branch_dir / f"{name}_clean_structural_features.csv",
        index=False,
    )

    count_by_group(table, ["organism_type", "gh_family"]).to_csv(
        branch_dir / f"{name}_counts_by_organism_and_family.csv",
        index=False,
    )

    numeric_summary(table, features).to_csv(
        branch_dir / f"{name}_overall_feature_summary.csv",
        index=False,
    )

    grouped_summary(table, ["gh_family"], features).to_csv(
        branch_dir / f"{name}_feature_summary_by_gh_family.csv",
        index=False,
    )

    grouped_summary(table, ["organism_type"], features).to_csv(
        branch_dir / f"{name}_feature_summary_by_organism_type.csv",
        index=False,
    )

    grouped_summary(table, ["organism_type", "gh_family"], features).to_csv(
        branch_dir / f"{name}_feature_summary_by_organism_and_family.csv",
        index=False,
    )

    return {
        "branch": name,
        "source_file": str(source_file),
        "valid_rows_used": len(table),
        "unique_accessions": table["accession_for_merge"].dropna().nunique()
        if "accession_for_merge" in table.columns
        else pd.NA,
        "features_summarised": ", ".join(features),
    }


def main():
    metadata = prepare_master_metadata()

    pdb = read_csv(PDB_FEATURE_FILE)
    pdb = add_accession_column(pdb, "pdb")
    pdb = merge_metadata(pdb, metadata)

    modeller, modeller_source = load_best_modeller_table()
    modeller = add_accession_column(modeller, "modeller")
    modeller = merge_metadata(modeller, metadata)

    summaries = [
        save_branch("pdb", pdb, PDB_FEATURE_FILE),
        save_branch("modeller", modeller, modeller_source),
    ]

    summary = pd.DataFrame(summaries)
    summary.to_csv(
        OUTPUT_DIR / "structural_feature_source_summary.csv",
        index=False,
    )

    readme = OUTPUT_DIR / "README.md"
    readme.write_text(
        """# Structural features without FoldX

This folder contains source-separated structural-feature summaries for PDB-supported and MODELLER-derived xylanase structures.

FoldX energy columns are intentionally excluded from this result section because FoldX stability outputs are reported separately in the FoldX section.

## Reported structural features

- chain length
- hydrogen-bond proxy count
- salt-bridge count
- disulfide count
- total solvent-accessible surface area
- per-residue normalised values where possible

## Interpretation

PDB and MODELLER features are kept separate because they represent different evidence levels. PDB-linked structures provide experimentally supported structural evidence. MODELLER-derived structures expand structural coverage but require source-aware interpretation.
"""
    )

    print("\nClean structural-feature summary:")
    print(summary.to_string(index=False))
    print(f"\nOutput directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
