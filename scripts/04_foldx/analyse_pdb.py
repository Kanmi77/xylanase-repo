#!/usr/bin/env python3
# Purpose: Analyse PDB FoldX results.

from pathlib import Path
import re

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"
OUTPUT_DIR = PROJECT_DIR / "results/foldx_pdb_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FILE = (
    PROJECT_DIR
    / "data/curated/xylanase_master_gh10_gh11_frozen_thesis_with_refseq_experimental_metadata.csv"
)

PDB_WT_FILE = (
    PROJECT_DIR
    / "data/curated/xylanase_structured_subset_with_foldx_norm.csv"
)

MUTATION_CANDIDATE_FILES = [
    PROJECT_DIR / "results/stability_ml_goal/datasets/mutation_prioritisation_panel.csv",
    PROJECT_DIR / "results/foldx_clean/tier2_ddg_ranked_industrial_annotated.csv",
    PROJECT_DIR / "results/foldx_clean/foldx_mutation_ddg_clean.csv",
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

WT_TOTAL_ENERGY_COLUMNS = [
    "foldx_wt_total_energy",
    "foldx_total_energy",
    "total_energy",
    "stability",
    "foldx_stability",
]

WT_ENERGY_PER_RESIDUE_COLUMNS = [
    "foldx_energy_per_residue",
    "energy_per_residue",
    "foldx_wt_energy_per_residue",
]

CHAIN_LENGTH_COLUMNS = [
    "chain_length",
    "sequence_length",
    "length",
]

DDG_COLUMNS = [
    "foldx_ddg",
    "ddg",
    "delta_delta_g",
    "foldx_delta_delta_g",
    "total_ddg",
]


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

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

    best_column = None
    best_count = 0
    best_values = None

    object_columns = [
        column for column in table.columns
        if table[column].dtype == "object"
    ]

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

    columns = [
        "accession_for_merge",
        "organism_type",
        "gh_family",
        "organism",
        "protein_name",
    ]

    columns = [
        column for column in columns
        if column in master.columns
    ]

    return master[columns].drop_duplicates("accession_for_merge")


def merge_metadata(table, metadata):
    table = table.copy()

    if "accession_for_merge" not in table.columns:
        table = add_accession_column(table, "table")

    table = table.merge(
        metadata,
        on="accession_for_merge",
        how="left",
        suffixes=("", "_master"),
    )

    for column in ["organism_type", "gh_family", "organism", "protein_name"]:
        master_column = f"{column}_master"

        if master_column in table.columns:
            if column in table.columns:
                table[column] = table[column].fillna(table[master_column])
            else:
                table[column] = table[master_column]

    return table


def standardise_pdb_wt_table(table):
    table = table.copy()

    total_energy_col = find_column(table, WT_TOTAL_ENERGY_COLUMNS)
    per_residue_col = find_column(table, WT_ENERGY_PER_RESIDUE_COLUMNS)
    chain_length_col = find_column(table, CHAIN_LENGTH_COLUMNS)

    if total_energy_col is None:
        raise KeyError("No FoldX WT total energy column was found.")

    table["foldx_wt_total_energy"] = pd.to_numeric(
        table[total_energy_col],
        errors="coerce",
    )

    if chain_length_col:
        table["chain_length"] = pd.to_numeric(
            table[chain_length_col],
            errors="coerce",
        )

    if per_residue_col:
        table["foldx_energy_per_residue"] = pd.to_numeric(
            table[per_residue_col],
            errors="coerce",
        )
    elif "chain_length" in table.columns:
        table["foldx_energy_per_residue"] = (
            table["foldx_wt_total_energy"] / table["chain_length"]
        )

    table["structure_source"] = "pdb"

    return table


def numeric_summary(table, columns):
    available = [
        column for column in columns
        if column in table.columns
    ]

    if not available:
        return pd.DataFrame()

    working = table.copy()

    for column in available:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    return (
        working[available]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .T
        .reset_index()
        .rename(columns={"index": "metric"})
    )


def grouped_numeric_summary(table, group_columns, numeric_columns):
    groups = [
        column for column in group_columns
        if column in table.columns
    ]

    numbers = [
        column for column in numeric_columns
        if column in table.columns
    ]

    if not groups or not numbers:
        return pd.DataFrame()

    working = table.copy()

    for column in numbers:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    summary = (
        working.groupby(groups, dropna=False)[numbers]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )

    summary.columns = [
        "_".join([str(part) for part in column if str(part) != ""]).strip("_")
        for column in summary.columns
    ]

    return summary


def count_by_group(table, group_columns):
    groups = [
        column for column in group_columns
        if column in table.columns
    ]

    if not groups:
        return pd.DataFrame()

    return (
        table.groupby(groups, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(groups)
    )


def is_pdb_mutation_row(table):
    possible_branch_columns = [
        "foldx_structure_branch",
        "foldx_source_branch",
        "evidence_source_branch",
        "foldx_branches",
        "structure_source",
        "source",
    ]

    mask = pd.Series(False, index=table.index)

    for column in possible_branch_columns:
        if column not in table.columns:
            continue

        text = table[column].astype(str).str.lower()

        mask = mask | text.str.contains("pdb", na=False)
        mask = mask | text.str.contains("experimental", na=False)

    return mask


def choose_pdb_mutation_table():
    source_rows = []
    best_table = None
    best_file = None
    best_ddg_col = None
    best_count = 0

    for path in MUTATION_CANDIDATE_FILES:
        if not path.exists():
            continue

        table = read_csv(path)
        ddg_col = find_column(table, DDG_COLUMNS)

        if ddg_col is None:
            for column in table.columns:
                lower = column.lower()

                if "ddg" in lower or "delta" in lower:
                    ddg_col = column
                    break

        if ddg_col is None:
            source_rows.append(
                {
                    "candidate_file": str(path),
                    "rows": len(table),
                    "ddg_column": "",
                    "pdb_rows": 0,
                    "selected": False,
                }
            )
            continue

        pdb_mask = is_pdb_mutation_row(table)
        pdb_rows = table[pdb_mask].copy()
        numeric_count = pd.to_numeric(
            pdb_rows[ddg_col],
            errors="coerce",
        ).notna().sum()

        selected = numeric_count > best_count

        if selected:
            best_table = pdb_rows
            best_file = path
            best_ddg_col = ddg_col
            best_count = numeric_count

        source_rows.append(
            {
                "candidate_file": str(path),
                "rows": len(table),
                "ddg_column": ddg_col,
                "pdb_rows": len(pdb_rows),
                "numeric_pdb_ddg_rows": numeric_count,
                "selected": selected,
            }
        )

    source_summary = pd.DataFrame(source_rows)
    source_summary.to_csv(OUTPUT_DIR / "foldx_pdb_mutation_source_summary.csv", index=False)

    if best_table is None or best_count == 0:
        return None, None, None

    return best_table, best_file, best_ddg_col


def standardise_mutation_table(table, ddg_col):
    table = table.copy()

    table["foldx_ddg"] = pd.to_numeric(table[ddg_col], errors="coerce")
    table = table.dropna(subset=["foldx_ddg"]).copy()

    if "accession_for_merge" not in table.columns:
        table = add_accession_column(table, "pdb_mutation")

    if "mutation" not in table.columns:
        mutation_col = find_column(
            table,
            [
                "mutation_code",
                "mutation_name",
                "mutation_string",
                "individual_list_mutation",
                "mutant",
                "substitution",
            ],
        )

        if mutation_col:
            table["mutation"] = table[mutation_col].astype(str)
        else:
            table["mutation"] = ""

    table["foldx_effect_class"] = "neutral_zero"
    table.loc[table["foldx_ddg"] < 0, "foldx_effect_class"] = "stabilising"
    table.loc[table["foldx_ddg"] > 0, "foldx_effect_class"] = "destabilising"

    return table


def main():
    metadata = prepare_master_metadata()

    pdb_wt = read_csv(PDB_WT_FILE)
    pdb_wt = add_accession_column(pdb_wt, "pdb_wt")
    pdb_wt = merge_metadata(pdb_wt, metadata)
    pdb_wt = standardise_pdb_wt_table(pdb_wt)

    wt_valid = pdb_wt.dropna(
        subset=[
            "foldx_wt_total_energy",
            "foldx_energy_per_residue",
        ],
        how="all",
    ).copy()

    wt_valid.to_csv(
        OUTPUT_DIR / "foldx_pdb_wt_stability_clean.csv",
        index=False,
    )

    wt_numeric_columns = [
        "foldx_wt_total_energy",
        "foldx_energy_per_residue",
        "chain_length",
    ]

    numeric_summary(
        wt_valid,
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_pdb_wt_overall_summary.csv",
        index=False,
    )

    count_by_group(
        wt_valid,
        ["organism_type", "gh_family"],
    ).to_csv(
        OUTPUT_DIR / "foldx_pdb_wt_counts_by_organism_and_family.csv",
        index=False,
    )

    grouped_numeric_summary(
        wt_valid,
        ["gh_family"],
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_pdb_wt_summary_by_gh_family.csv",
        index=False,
    )

    grouped_numeric_summary(
        wt_valid,
        ["organism_type"],
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_pdb_wt_summary_by_organism_type.csv",
        index=False,
    )

    grouped_numeric_summary(
        wt_valid,
        ["organism_type", "gh_family"],
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_pdb_wt_summary_by_organism_and_family.csv",
        index=False,
    )

    wt_valid.sort_values(
        "foldx_energy_per_residue",
        ascending=True,
    ).head(30).to_csv(
        OUTPUT_DIR / "top30_pdb_lowest_foldx_energy_per_residue.csv",
        index=False,
    )

    wt_valid.sort_values(
        "foldx_energy_per_residue",
        ascending=False,
    ).head(30).to_csv(
        OUTPUT_DIR / "top30_pdb_highest_foldx_energy_per_residue.csv",
        index=False,
    )

    mutation_table, mutation_file, ddg_col = choose_pdb_mutation_table()

    if mutation_table is not None:
        mutations = standardise_mutation_table(mutation_table, ddg_col)
        mutations = merge_metadata(mutations, metadata)

        mutations.to_csv(
            OUTPUT_DIR / "foldx_pdb_mutation_ddg_clean.csv",
            index=False,
        )

        numeric_summary(
            mutations,
            ["foldx_ddg"],
        ).to_csv(
            OUTPUT_DIR / "foldx_pdb_mutation_ddg_overall_summary.csv",
            index=False,
        )

        effect_counts = (
            mutations["foldx_effect_class"]
            .value_counts(dropna=False)
            .rename_axis("foldx_effect_class")
            .reset_index(name="count")
        )

        effect_counts["percent"] = (
            effect_counts["count"] / effect_counts["count"].sum() * 100
        )

        effect_counts.to_csv(
            OUTPUT_DIR / "foldx_pdb_mutation_effect_counts.csv",
            index=False,
        )

        grouped_numeric_summary(
            mutations,
            ["organism_type", "gh_family"],
            ["foldx_ddg"],
        ).to_csv(
            OUTPUT_DIR / "foldx_pdb_mutation_ddg_by_organism_and_family.csv",
            index=False,
        )

        mutations.sort_values("foldx_ddg", ascending=True).head(30).to_csv(
            OUTPUT_DIR / "top30_pdb_stabilising_foldx_mutations.csv",
            index=False,
        )

        mutations.sort_values("foldx_ddg", ascending=False).head(30).to_csv(
            OUTPUT_DIR / "top30_pdb_destabilising_foldx_mutations.csv",
            index=False,
        )

    readme = OUTPUT_DIR / "README.md"
    readme.write_text(
        """# PDB-only FoldX results

This folder contains FoldX summaries for the PDB-supported structural branch only.

## Wild-type stability

The wild-type FoldX result is forced to use:

`data/curated/xylanase_structured_subset_with_foldx_norm.csv`

## Mutation screening

PDB mutation rows are included only if a mutation table contains rows marked as PDB or experimental structure branch.

## Interpretation note

This section should be reported separately from the MODELLER FoldX result because PDB-supported structures and MODELLER-derived structures represent different structural evidence levels.
"""
    )

    report_lines = []
    report_lines.append("# PDB-only FoldX results\n")
    report_lines.append(f"- PDB WT source file: `{PDB_WT_FILE}`")
    report_lines.append(f"- Valid PDB WT rows: {len(wt_valid)}")

    report_lines.append("\n## PDB WT overall summary\n")
    report_lines.append(
        numeric_summary(wt_valid, wt_numeric_columns).to_markdown(index=False)
    )

    if mutation_table is not None:
        report_lines.append("\n## PDB mutation ΔΔG summary\n")
        report_lines.append(f"- Mutation source file: `{mutation_file}`")
        report_lines.append(f"- ΔΔG column used: `{ddg_col}`")
        report_lines.append(f"- PDB mutation rows: {len(mutations)}")
        report_lines.append(
            numeric_summary(mutations, ["foldx_ddg"]).to_markdown(index=False)
        )
        report_lines.append("\n## PDB mutation effect counts\n")
        report_lines.append(effect_counts.to_markdown(index=False))
    else:
        report_lines.append("\n## PDB mutation ΔΔG summary\n")
        report_lines.append(
            "- No PDB-specific mutation ΔΔG rows were detected automatically."
        )

    (OUTPUT_DIR / "PDB_FOLDX_RESULTS_SUMMARY.md").write_text(
        "\n".join(report_lines)
    )

    print("\nPDB FoldX WT summary:")
    print(numeric_summary(wt_valid, wt_numeric_columns).to_string(index=False))

    if mutation_table is not None:
        print("\nPDB FoldX mutation ΔΔG summary:")
        print(numeric_summary(mutations, ["foldx_ddg"]).to_string(index=False))

        print("\nPDB FoldX mutation effect counts:")
        print(effect_counts.to_string(index=False))
    else:
        print("\nNo PDB-specific mutation ΔΔG rows were detected.")

    print(f"\nOutput directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
