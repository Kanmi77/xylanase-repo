#!/usr/bin/env python3
# Purpose: Compare PDB and modelled structures.

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"
OUTPUT_DIR = PROJECT_DIR / "results/structural_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


MASTER_FILES = [
    PROJECT_DIR / "data/curated/master_metadata.csv",
    PROJECT_DIR / "data/curated/xylanase_master_all_curated_with_brenda.csv",
    PROJECT_DIR / "data/curated/xylanase_master_all_curated.csv",
]

PDB_FILES = [
    PROJECT_DIR / "data/curated/xylanase_structured_subset_with_foldx_norm.csv",
    PROJECT_DIR / "data/curated/xylanase_structured_subset_with_foldx.csv",
    PROJECT_DIR / "data/curated/xylanase_structured_subset.csv",
    PROJECT_DIR / "results/structures/structure_manifest.csv",
]

MODELLER_FILES = [
    PROJECT_DIR / "results/foldx/modeller/structural_features_with_length.csv",
    PROJECT_DIR / "results/foldx/modeller/structural_features_full.csv",
]

COMBINED_FILES = [
    PROJECT_DIR / "results/structure_features/combined_structural_features.csv",
    PROJECT_DIR / "results/structure_features/combined_structural_features_group_summary.csv",
    PROJECT_DIR / "results/structural_features_foldx_summary/structural_features_with_calculated_normalised_values.csv",
]

TMALIGN_FILES = [
    PROJECT_DIR / "results/structures/tmalign_best_reference_per_model.csv",
    PROJECT_DIR / "results/structures/tmalign_results.csv",
    PROJECT_DIR / "results/reports/tmalign_best_reference_overall_summary.csv",
]

REPORT_FILES = [
    PROJECT_DIR / "results/reports/structural_features_overall_summary.csv",
    PROJECT_DIR / "results/reports/structural_features_by_organism.csv",
    PROJECT_DIR / "results/reports/structural_features_by_gh_family.csv",
    PROJECT_DIR / "results/reports/structural_features_by_organism_gh.csv",
    PROJECT_DIR / "results/structural_features_foldx_summary/structural_feature_raw_summary.csv",
    PROJECT_DIR / "results/structural_features_foldx_summary/structural_feature_normalised_summary.csv",
]


def first_existing(paths):
    for path in paths:
        if path.exists():
            return path

    return None


def read_csv_if_exists(path):
    if path is None or not path.exists():
        return pd.DataFrame()

    return pd.read_csv(path)


def normalise_column_names(table):
    table = table.copy()
    table.columns = [str(column).strip() for column in table.columns]
    return table


def find_column(table, candidates):
    lower_map = {column.lower(): column for column in table.columns}

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    return None


def count_unique(table, candidates):
    column = find_column(table, candidates)

    if column is None:
        return pd.NA

    return table[column].dropna().astype(str).replace("", pd.NA).dropna().nunique()


def summarise_numeric(table, group_columns, numeric_columns, output_file):
    available_groups = [
        column for column in group_columns
        if column in table.columns
    ]

    available_numeric = [
        column for column in numeric_columns
        if column in table.columns
    ]

    if not available_groups or not available_numeric:
        return pd.DataFrame()

    working = table.copy()

    for column in available_numeric:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    summary = (
        working.groupby(available_groups, dropna=False)[available_numeric]
        .agg(["count", "mean", "median", "min", "max"])
        .reset_index()
    )

    summary.columns = [
        "_".join([str(item) for item in column if str(item) != ""]).strip("_")
        for column in summary.columns
    ]

    summary.to_csv(output_file, index=False)
    return summary


def summarise_counts(table, group_columns, output_file):
    available_groups = [
        column for column in group_columns
        if column in table.columns
    ]

    if not available_groups:
        return pd.DataFrame()

    summary = (
        table.groupby(available_groups, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(available_groups)
    )

    summary.to_csv(output_file, index=False)
    return summary


def main():
    master_path = first_existing(MASTER_FILES)
    pdb_path = first_existing(PDB_FILES)
    modeller_path = first_existing(MODELLER_FILES)
    combined_path = first_existing(COMBINED_FILES)
    tmalign_path = first_existing(TMALIGN_FILES)

    master = normalise_column_names(read_csv_if_exists(master_path))
    pdb = normalise_column_names(read_csv_if_exists(pdb_path))
    modeller = normalise_column_names(read_csv_if_exists(modeller_path))
    combined = normalise_column_names(read_csv_if_exists(combined_path))
    tmalign = normalise_column_names(read_csv_if_exists(tmalign_path))

    overview_rows = []

    overview_rows.append(
        {
            "section": "master_dataset",
            "file_used": str(master_path) if master_path else "missing",
            "rows": len(master),
            "unique_accessions": count_unique(
                master,
                ["uniprot_accession", "accession", "Entry"],
            ),
            "unique_sequences": count_unique(
                master,
                ["sequence", "sequence_clean"],
            ),
        }
    )

    overview_rows.append(
        {
            "section": "pdb_supported_branch",
            "file_used": str(pdb_path) if pdb_path else "missing",
            "rows": len(pdb),
            "unique_pdb_ids": count_unique(pdb, ["pdb_id"]),
            "unique_accessions": count_unique(
                pdb,
                ["uniprot_accession", "accession", "Entry"],
            ),
            "parsed_chains": int(
                pd.to_numeric(
                    pdb.get("chain_length", pd.Series(dtype=float)),
                    errors="coerce",
                ).gt(0).sum()
            )
            if "chain_length" in pdb.columns
            else pd.NA,
        }
    )

    overview_rows.append(
        {
            "section": "modeller_branch",
            "file_used": str(modeller_path) if modeller_path else "missing",
            "rows": len(modeller),
            "unique_accessions": count_unique(
                modeller,
                ["uniprot_accession", "accession", "Entry", "query_accession"],
            ),
            "unique_models": count_unique(
                modeller,
                ["model_path", "pdb_path", "structure_file", "model_file"],
            ),
        }
    )

    overview_rows.append(
        {
            "section": "combined_structural_features",
            "file_used": str(combined_path) if combined_path else "missing",
            "rows": len(combined),
            "unique_accessions": count_unique(
                combined,
                ["uniprot_accession", "accession", "Entry", "query_accession"],
            ),
            "structure_sources": count_unique(
                combined,
                ["structure_source", "source", "model_source"],
            ),
        }
    )

    overview_rows.append(
        {
            "section": "tmalign_validation",
            "file_used": str(tmalign_path) if tmalign_path else "missing",
            "rows": len(tmalign),
            "unique_queries": count_unique(
                tmalign,
                ["query", "query_id", "query_accession", "model_id"],
            ),
        }
    )

    overview = pd.DataFrame(overview_rows)
    overview.to_csv(OUTPUT_DIR / "structural_dataset_overview.csv", index=False)

    # PDB branch summaries
    if not pdb.empty:
        summarise_counts(
            pdb,
            ["organism_type", "gh_family"],
            OUTPUT_DIR / "pdb_counts_by_organism_and_family.csv",
        )

        summarise_counts(
            pdb,
            ["file_type"],
            OUTPUT_DIR / "pdb_counts_by_file_type.csv",
        )

        pdb_numeric = [
            "chain_length",
            "hbond_count",
            "salt_bridge_count",
            "disulfide_count",
            "sasa_total",
            "total_sasa",
            "foldx_wt_total_energy",
            "foldx_energy_per_residue",
        ]

        summarise_numeric(
            pdb,
            ["organism_type", "gh_family"],
            pdb_numeric,
            OUTPUT_DIR / "pdb_structural_feature_summary_by_group.csv",
        )

    # MODELLER branch summaries
    if not modeller.empty:
        summarise_counts(
            modeller,
            ["organism_type", "gh_family"],
            OUTPUT_DIR / "modeller_counts_by_organism_and_family.csv",
        )

        modeller_numeric = [
            "sequence_length",
            "chain_length",
            "length",
            "hbond_count",
            "salt_bridge_count",
            "disulfide_count",
            "sasa_total",
            "total_sasa",
            "sasa",
            "foldx_wt_total_energy",
            "foldx_energy_per_residue",
        ]

        summarise_numeric(
            modeller,
            ["organism_type", "gh_family"],
            modeller_numeric,
            OUTPUT_DIR / "modeller_structural_feature_summary_by_group.csv",
        )

    # Combined branch summaries
    if not combined.empty:
        source_column = find_column(
            combined,
            ["structure_source", "source", "model_source"],
        )

        if source_column:
            summarise_counts(
                combined,
                [source_column],
                OUTPUT_DIR / "combined_counts_by_structure_source.csv",
            )

            summarise_counts(
                combined,
                [source_column, "organism_type", "gh_family"],
                OUTPUT_DIR / "combined_counts_by_source_organism_family.csv",
            )

        combined_numeric = [
            "sequence_length",
            "chain_length",
            "length",
            "hbond_count",
            "hbond_count_per_residue",
            "salt_bridge_count",
            "salt_bridge_count_per_residue",
            "disulfide_count",
            "sasa_total",
            "total_sasa",
            "sasa_per_residue",
            "foldx_wt_total_energy",
            "foldx_energy_per_residue",
        ]

        group_columns = ["organism_type", "gh_family"]

        if source_column:
            group_columns = [source_column, "organism_type", "gh_family"]

        summarise_numeric(
            combined,
            group_columns,
            combined_numeric,
            OUTPUT_DIR / "combined_structural_feature_summary_by_source_group.csv",
        )

    # TM-align summaries
    if not tmalign.empty:
        numeric_columns = []

        for column in tmalign.columns:
            lower = column.lower()

            if (
                "tm" in lower
                or "rmsd" in lower
                or "identity" in lower
                or "aligned" in lower
            ):
                numeric_columns.append(column)

        for column in numeric_columns:
            tmalign[column] = pd.to_numeric(tmalign[column], errors="coerce")

        tmalign_numeric_summary = (
            tmalign[numeric_columns]
            .agg(["count", "mean", "median", "min", "max"])
            .T
            .reset_index()
            .rename(columns={"index": "metric"})
            if numeric_columns
            else pd.DataFrame()
        )

        tmalign_numeric_summary.to_csv(
            OUTPUT_DIR / "tmalign_numeric_summary.csv",
            index=False,
        )

        quality_column = find_column(
            tmalign,
            ["quality_class", "tm_quality_class", "match_quality"],
        )

        if quality_column:
            summarise_counts(
                tmalign,
                [quality_column],
                OUTPUT_DIR / "tmalign_quality_class_counts.csv",
            )

    # Copy already prepared report tables if present
    copied_report_rows = []

    for report_file in REPORT_FILES:
        if report_file.exists():
            target = OUTPUT_DIR / report_file.name
            table = pd.read_csv(report_file)
            table.to_csv(target, index=False)

            copied_report_rows.append(
                {
                    "source_file": str(report_file),
                    "copied_to": str(target),
                    "rows": len(table),
                }
            )

    copied_reports = pd.DataFrame(copied_report_rows)
    copied_reports.to_csv(OUTPUT_DIR / "copied_existing_structural_reports.csv", index=False)

    # Markdown summary
    summary_md = OUTPUT_DIR / "STRUCTURAL_RESULTS_SUMMARY.md"

    with summary_md.open("w") as handle:
        handle.write("# Structural results summary\n\n")
        handle.write("## Files used\n\n")
        handle.write(overview.to_markdown(index=False))
        handle.write("\n\n")

        handle.write("## Interpretation guide\n\n")
        handle.write(
            "- The PDB-supported branch represents experimentally resolved structural records or PDB-linked structures.\n"
        )
        handle.write(
            "- The MODELLER branch represents homology-modelled structures used to increase structural coverage.\n"
        )
        handle.write(
            "- Combined structural-feature outputs should be used for source-stratified comparison, while PDB and MODELLER should also be reported separately.\n"
        )
        handle.write(
            "- TM-align outputs are used as structural validation evidence for model-reference similarity.\n"
        )

    print("\nStructural analysis files used:")
    print(overview.to_string(index=False))

    print(f"\nResult tables written to: {OUTPUT_DIR}")
    print(f"Summary written to: {summary_md}")


if __name__ == "__main__":
    main()
