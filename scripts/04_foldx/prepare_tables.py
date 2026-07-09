#!/usr/bin/env python3
# Purpose: Prepare FoldX result tables.

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"

MODELLER_FILE = PROJECT_DIR / "results/foldx_results/foldx_mutation_ddg_clean.csv"
PDB_FILE = PROJECT_DIR / "results/foldx_pdb_results/foldx_pdb_mutation_ddg_clean.csv"

OUTPUT_DIR = PROJECT_DIR / "results/foldx_mutation_thesis_tables"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


KEEP_COLUMNS = [
    "protein",
    "uniprot_accession",
    "mutation",
    "candidate_id",
    "foldx_ddg",
    "foldx_effect_class",
    "foldx_effect",
    "foldx_structure_branch",
    "foldx_source_file",
    "organism",
    "organism_type",
    "gh_family",
    "pdb_id",
    "structure_id",
    "mutant_model",
    "docking_delta_mean",
    "docking_delta_best",
    "docking_delta_worst",
    "docking_improved_fraction",
    "docking_retained_fraction",
    "activity_proxy_status",
    "recommendation_class",
    "evidence_candidate_tier",
    "evidence_source_branch",
]


def read_mutation_table(path, branch_label):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    table = pd.read_csv(path)
    table.columns = [str(column).strip() for column in table.columns]

    table["result_branch"] = branch_label
    table["foldx_ddg"] = pd.to_numeric(table["foldx_ddg"], errors="coerce")
    table = table.dropna(subset=["foldx_ddg"]).copy()

    if "foldx_effect_class" not in table.columns:
        table["foldx_effect_class"] = "neutral_zero"
        table.loc[table["foldx_ddg"] < 0, "foldx_effect_class"] = "stabilising"
        table.loc[table["foldx_ddg"] > 0, "foldx_effect_class"] = "destabilising"

    return table


def numeric_summary(table, group_columns):
    summary = (
        table.groupby(group_columns, dropna=False)["foldx_ddg"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )

    summary["range"] = summary["max"] - summary["min"]

    return summary


def effect_counts(table, group_columns):
    summary = (
        table.groupby(group_columns + ["foldx_effect_class"], dropna=False)
        .size()
        .reset_index(name="count")
    )

    total = (
        summary.groupby(group_columns, dropna=False)["count"]
        .transform("sum")
    )

    summary["percent"] = summary["count"] / total * 100

    return summary


def clean_candidate_table(table):
    columns = ["result_branch"]

    for column in KEEP_COLUMNS:
        if column in table.columns:
            columns.append(column)

    clean = table[columns].copy()

    clean = clean.sort_values(
        ["result_branch", "foldx_ddg"],
        ascending=[True, True],
    )

    return clean


def make_top_table(table, branch_label, stabilising=True, n=15):
    branch_table = table[table["result_branch"] == branch_label].copy()

    if stabilising:
        branch_table = branch_table[branch_table["foldx_ddg"] < 0].copy()
        branch_table = branch_table.sort_values("foldx_ddg", ascending=True)
    else:
        branch_table = branch_table[branch_table["foldx_ddg"] > 0].copy()
        branch_table = branch_table.sort_values("foldx_ddg", ascending=False)

    clean = clean_candidate_table(branch_table).head(n)

    return clean


def main():
    modeller = read_mutation_table(MODELLER_FILE, "MODELLER-derived")
    pdb = read_mutation_table(PDB_FILE, "PDB-supported")

    combined = pd.concat([modeller, pdb], ignore_index=True)

    combined_clean = clean_candidate_table(combined)
    combined_clean.to_csv(
        OUTPUT_DIR / "table_00_all_foldx_mutation_rows_clean.csv",
        index=False,
    )

    overall_summary = numeric_summary(combined, ["result_branch"])
    overall_summary.to_csv(
        OUTPUT_DIR / "table_01_foldx_mutation_ddg_summary_by_branch.csv",
        index=False,
    )

    branch_effects = effect_counts(combined, ["result_branch"])
    branch_effects.to_csv(
        OUTPUT_DIR / "table_02_foldx_mutation_effect_counts_by_branch.csv",
        index=False,
    )

    group_summary = numeric_summary(
        combined,
        ["result_branch", "organism_type", "gh_family"],
    )
    group_summary.to_csv(
        OUTPUT_DIR / "table_03_foldx_mutation_ddg_by_branch_organism_family.csv",
        index=False,
    )

    group_effects = effect_counts(
        combined,
        ["result_branch", "organism_type", "gh_family"],
    )
    group_effects.to_csv(
        OUTPUT_DIR / "table_04_foldx_mutation_effect_counts_by_branch_organism_family.csv",
        index=False,
    )

    stabilising_summary = numeric_summary(
        combined[combined["foldx_effect_class"] == "stabilising"],
        ["result_branch"],
    )
    stabilising_summary.to_csv(
        OUTPUT_DIR / "table_05_stabilising_ddg_summary_by_branch.csv",
        index=False,
    )

    destabilising_summary = numeric_summary(
        combined[combined["foldx_effect_class"] == "destabilising"],
        ["result_branch"],
    )
    destabilising_summary.to_csv(
        OUTPUT_DIR / "table_06_destabilising_ddg_summary_by_branch.csv",
        index=False,
    )

    stabilising_group_summary = numeric_summary(
        combined[combined["foldx_effect_class"] == "stabilising"],
        ["result_branch", "organism_type", "gh_family"],
    )
    stabilising_group_summary.to_csv(
        OUTPUT_DIR / "table_07_stabilising_ddg_by_branch_organism_family.csv",
        index=False,
    )

    destabilising_group_summary = numeric_summary(
        combined[combined["foldx_effect_class"] == "destabilising"],
        ["result_branch", "organism_type", "gh_family"],
    )
    destabilising_group_summary.to_csv(
        OUTPUT_DIR / "table_08_destabilising_ddg_by_branch_organism_family.csv",
        index=False,
    )

    top_modeller = make_top_table(
        combined,
        "MODELLER-derived",
        stabilising=True,
        n=15,
    )
    top_modeller.to_csv(
        OUTPUT_DIR / "table_09_top15_modeller_stabilising_mutations.csv",
        index=False,
    )

    top_pdb = make_top_table(
        combined,
        "PDB-supported",
        stabilising=True,
        n=15,
    )
    top_pdb.to_csv(
        OUTPUT_DIR / "table_10_top15_pdb_stabilising_mutations.csv",
        index=False,
    )

    top_destabilising_modeller = make_top_table(
        combined,
        "MODELLER-derived",
        stabilising=False,
        n=15,
    )
    top_destabilising_modeller.to_csv(
        OUTPUT_DIR / "table_11_top15_modeller_destabilising_mutations.csv",
        index=False,
    )

    top_destabilising_pdb = make_top_table(
        combined,
        "PDB-supported",
        stabilising=False,
        n=15,
    )
    top_destabilising_pdb.to_csv(
        OUTPUT_DIR / "table_12_top15_pdb_destabilising_mutations.csv",
        index=False,
    )

    combined_top = (
        combined[combined["foldx_ddg"] < 0]
        .sort_values("foldx_ddg", ascending=True)
    )

    clean_candidate_table(combined_top).head(20).to_csv(
        OUTPUT_DIR / "table_13_top20_combined_stabilising_mutations.csv",
        index=False,
    )

    if "recommendation_class" in combined.columns:
        recommendation_summary = (
            combined.groupby(
                ["result_branch", "recommendation_class", "foldx_effect_class"],
                dropna=False,
            )
            .size()
            .reset_index(name="count")
        )

        recommendation_summary.to_csv(
            OUTPUT_DIR / "table_14_mutation_recommendation_class_summary.csv",
            index=False,
        )

    if "activity_proxy_status" in combined.columns:
        activity_summary = (
            combined.groupby(
                ["result_branch", "activity_proxy_status", "foldx_effect_class"],
                dropna=False,
            )
            .size()
            .reset_index(name="count")
        )

        activity_summary.to_csv(
            OUTPUT_DIR / "table_15_mutation_activity_proxy_summary.csv",
            index=False,
        )

    readme = OUTPUT_DIR / "README.md"
    readme.write_text(
        """# FoldX mutation thesis tables

These tables summarise the FoldX mutation-screening results for thesis reporting.

## Important note

The mutation outputs are reported as row-level FoldX mutation panels, not deduplicated unique mutation candidates.

## Numeric summaries

Numeric ΔΔG summaries include:

- count
- mean
- median
- standard deviation
- minimum
- maximum
- range

## Table guide

- table_01: ΔΔG summary by structure branch
- table_02: stabilising/destabilising counts by branch
- table_03: ΔΔG summary by branch, organism type and GH family
- table_04: stabilising/destabilising counts by branch, organism type and GH family
- table_05: stabilising-only ΔΔG summary by branch
- table_06: destabilising-only ΔΔG summary by branch
- table_07: stabilising-only ΔΔG summary by branch, organism type and GH family
- table_08: destabilising-only ΔΔG summary by branch, organism type and GH family
- table_09: top 15 MODELLER-supported stabilising mutations
- table_10: top 15 PDB-supported stabilising mutations
- table_11: top 15 MODELLER-supported destabilising mutations
- table_12: top 15 PDB-supported destabilising mutations
- table_13: top 20 stabilising mutations across both branches
- table_14: recommendation-class summary, if available
- table_15: activity-proxy summary, if available
"""
    )

    print()
    print("FoldX mutation thesis tables created with SD, median and range.")
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    print("Overall ΔΔG summary:")
    print(overall_summary.to_string(index=False))
    print()
    print("Effect counts:")
    print(branch_effects.to_string(index=False))
    print()
    print("Stabilising-only summary:")
    print(stabilising_summary.to_string(index=False))
    print()
    print("Destabilising-only summary:")
    print(destabilising_summary.to_string(index=False))


if __name__ == "__main__":
    main()
