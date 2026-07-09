#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


PROJECT_DIR = Path.home() / "xylanase-thesis"

OUTPUT_DIR = PROJECT_DIR / "results/foldx_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


MASTER_FILE = (
    PROJECT_DIR
    / "data/curated/xylanase_master_gh10_gh11_frozen_thesis_with_refseq_experimental_metadata.csv"
)

WT_CANDIDATE_FILES = [
    PROJECT_DIR / "data/curated/xylanase_structured_subset_with_foldx_norm.csv",
    PROJECT_DIR / "data/curated/xylanase_structured_subset_with_foldx.csv",
    PROJECT_DIR / "results/structural_features_foldx_summary/structural_features_with_calculated_normalised_values.csv",
    PROJECT_DIR / "results/foldx4_wt_stability.csv",
    PROJECT_DIR / "results/foldx/foldx4_wt_stability.csv",
    PROJECT_DIR / "results/foldx_clean/foldx4_wt_stability.csv",
]

MUTATION_SEARCH_ROOTS = [
    PROJECT_DIR / "results/foldx_clean",
    PROJECT_DIR / "results",
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

SOURCE_COLUMNS = [
    "structure_source",
    "source",
    "model_source",
    "foldx_source",
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
    "ddg",
    "delta_delta_g",
    "foldx_ddg",
    "foldx_delta_delta_g",
    "total_ddg",
    "difference",
    "energy_change",
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


def choose_wt_file():
    source_rows = []

    for path in WT_CANDIDATE_FILES:
        if not path.exists():
            continue

        table = read_csv(path)

        total_energy_col = find_column(table, WT_TOTAL_ENERGY_COLUMNS)
        per_residue_col = find_column(table, WT_ENERGY_PER_RESIDUE_COLUMNS)

        score = 0

        if total_energy_col:
            score += 10

        if per_residue_col:
            score += 10

        source_rows.append(
            {
                "candidate_file": str(path),
                "rows": len(table),
                "total_energy_column": total_energy_col or "",
                "energy_per_residue_column": per_residue_col or "",
                "score": score,
            }
        )

    source_summary = pd.DataFrame(source_rows)
    source_summary.to_csv(OUTPUT_DIR / "foldx_wt_source_summary.csv", index=False)

    if source_summary.empty:
        raise FileNotFoundError("No FoldX WT candidate files found.")

    source_summary = source_summary.sort_values(
        ["score", "rows"],
        ascending=[False, False],
    )

    selected_path = Path(source_summary.iloc[0]["candidate_file"])

    if int(source_summary.iloc[0]["score"]) == 0:
        raise ValueError(
            "FoldX WT files were found, but no recognised energy columns were detected."
        )

    print(f"Selected FoldX WT file: {selected_path}")

    return selected_path


def standardise_wt_table(table):
    table = table.copy()

    total_energy_col = find_column(table, WT_TOTAL_ENERGY_COLUMNS)
    per_residue_col = find_column(table, WT_ENERGY_PER_RESIDUE_COLUMNS)
    chain_length_col = find_column(table, CHAIN_LENGTH_COLUMNS)
    source_col = find_column(table, SOURCE_COLUMNS)

    if total_energy_col:
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
    elif (
        "foldx_wt_total_energy" in table.columns
        and "chain_length" in table.columns
    ):
        table["foldx_energy_per_residue"] = (
            table["foldx_wt_total_energy"] / table["chain_length"]
        )

    if source_col:
        table["structure_source"] = table[source_col].astype(str)
    elif "pdb_id" in table.columns:
        table["structure_source"] = "pdb"
    else:
        table["structure_source"] = "unknown"

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


def find_mutation_files():
    files = []

    for root in MUTATION_SEARCH_ROOTS:
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            name = path.name.lower()

            if path.suffix.lower() == ".csv" and (
                "ddg" in name
                or "mutation" in name
                or "buildmodel" in name
                or "foldx" in name
            ):
                files.append(path)

    return sorted(set(files))


def score_mutation_file(path):
    try:
        table = read_csv(path)
    except Exception:
        return None

    ddg_col = find_column(table, DDG_COLUMNS)

    if not ddg_col:
        # Try fuzzy detection.
        for column in table.columns:
            lower = column.lower()

            if "ddg" in lower or "delta" in lower:
                ddg_col = column
                break

    if not ddg_col:
        return None

    values = pd.to_numeric(table[ddg_col], errors="coerce")
    numeric_count = values.notna().sum()

    if numeric_count == 0:
        return None

    return {
        "path": path,
        "rows": len(table),
        "ddg_column": ddg_col,
        "numeric_ddg_count": numeric_count,
        "score": numeric_count,
    }


def choose_mutation_file():
    candidates = []

    for path in find_mutation_files():
        candidate = score_mutation_file(path)

        if candidate:
            candidates.append(candidate)

    source_summary = pd.DataFrame(
        [
            {
                "candidate_file": str(item["path"]),
                "rows": item["rows"],
                "ddg_column": item["ddg_column"],
                "numeric_ddg_count": item["numeric_ddg_count"],
                "score": item["score"],
            }
            for item in candidates
        ]
    )

    source_summary.to_csv(OUTPUT_DIR / "foldx_mutation_source_summary.csv", index=False)

    if source_summary.empty:
        print("No FoldX mutation ΔΔG file detected.")
        return None, None

    source_summary = source_summary.sort_values(
        ["score", "rows"],
        ascending=[False, False],
    )

    selected_path = Path(source_summary.iloc[0]["candidate_file"])
    selected_ddg_col = source_summary.iloc[0]["ddg_column"]

    print(f"Selected FoldX mutation file: {selected_path}")
    print(f"Selected ΔΔG column: {selected_ddg_col}")

    return selected_path, selected_ddg_col


def standardise_mutation_table(table, ddg_col):
    table = table.copy()

    table["foldx_ddg"] = pd.to_numeric(table[ddg_col], errors="coerce")

    table = table.dropna(subset=["foldx_ddg"]).copy()

    if "accession_for_merge" not in table.columns:
        table = add_accession_column(table, "mutation")

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

    # WT stability
    wt_file = choose_wt_file()
    wt = read_csv(wt_file)
    wt = add_accession_column(wt, "foldx_wt")
    wt = merge_metadata(wt, metadata)
    wt = standardise_wt_table(wt)

    wt_valid = wt.dropna(
        subset=[
            column
            for column in ["foldx_wt_total_energy", "foldx_energy_per_residue"]
            if column in wt.columns
        ],
        how="all",
    ).copy()

    wt_valid.to_csv(
        OUTPUT_DIR / "foldx_wt_stability_clean.csv",
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
        OUTPUT_DIR / "foldx_wt_overall_summary.csv",
        index=False,
    )

    count_by_group(
        wt_valid,
        ["structure_source"],
    ).to_csv(
        OUTPUT_DIR / "foldx_wt_counts_by_structure_source.csv",
        index=False,
    )

    grouped_numeric_summary(
        wt_valid,
        ["structure_source"],
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_wt_summary_by_structure_source.csv",
        index=False,
    )

    grouped_numeric_summary(
        wt_valid,
        ["organism_type", "gh_family"],
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_wt_summary_by_organism_and_family.csv",
        index=False,
    )

    grouped_numeric_summary(
        wt_valid,
        ["structure_source", "organism_type", "gh_family"],
        wt_numeric_columns,
    ).to_csv(
        OUTPUT_DIR / "foldx_wt_summary_by_source_organism_family.csv",
        index=False,
    )

    wt_sorted = wt_valid.sort_values(
        "foldx_energy_per_residue",
        ascending=True,
    )

    wt_sorted.head(30).to_csv(
        OUTPUT_DIR / "top30_lowest_foldx_energy_per_residue.csv",
        index=False,
    )

    wt_sorted.tail(30).to_csv(
        OUTPUT_DIR / "top30_highest_foldx_energy_per_residue.csv",
        index=False,
    )

    # Mutation ΔΔG
    mutation_file, ddg_col = choose_mutation_file()

    if mutation_file is not None:
        mutations = read_csv(mutation_file)
        mutations = standardise_mutation_table(mutations, ddg_col)
        mutations = merge_metadata(mutations, metadata)

        mutations.to_csv(
            OUTPUT_DIR / "foldx_mutation_ddg_clean.csv",
            index=False,
        )

        numeric_summary(
            mutations,
            ["foldx_ddg"],
        ).to_csv(
            OUTPUT_DIR / "foldx_mutation_ddg_overall_summary.csv",
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
            OUTPUT_DIR / "foldx_mutation_effect_counts.csv",
            index=False,
        )

        grouped_numeric_summary(
            mutations,
            ["organism_type", "gh_family"],
            ["foldx_ddg"],
        ).to_csv(
            OUTPUT_DIR / "foldx_mutation_ddg_by_organism_and_family.csv",
            index=False,
        )

        count_by_group(
            mutations,
            ["organism_type", "gh_family", "foldx_effect_class"],
        ).to_csv(
            OUTPUT_DIR / "foldx_mutation_effect_counts_by_group.csv",
            index=False,
        )

        mutations.sort_values("foldx_ddg", ascending=True).head(30).to_csv(
            OUTPUT_DIR / "top30_stabilising_foldx_mutations.csv",
            index=False,
        )

        mutations.sort_values("foldx_ddg", ascending=False).head(30).to_csv(
            OUTPUT_DIR / "top30_destabilising_foldx_mutations.csv",
            index=False,
        )

    readme = OUTPUT_DIR / "README.md"
    readme.write_text(
        """# FoldX stability and mutation-screening results

This folder contains FoldX result summaries for the xylanase thesis.

## Scope

FoldX results are reported separately from structural features, docking, molecular dynamics, and machine learning.

## Wild-type stability

The wild-type stability section summarises:

- FoldX total wild-type stability energy
- FoldX energy per residue
- source-aware summaries by structure source, organism type, and GH family

## Mutation screening

The mutation-screening section summarises FoldX BuildModel ΔΔG values where available.

Mutation classes:

- stabilising: ΔΔG < 0
- destabilising: ΔΔG > 0
- neutral_zero: ΔΔG = 0

## Interpretation note

FoldX values are computational energy estimates. They are used as a screening and prioritisation metric, not as direct experimental proof of thermostability.
"""
    )

    report_lines = []

    report_lines.append("# FoldX results summary\n")
    report_lines.append("## Wild-type FoldX stability\n")
    report_lines.append(f"- WT source file: `{wt_file}`")
    report_lines.append(f"- Valid WT rows: {len(wt_valid)}")

    report_lines.append("\n### WT overall summary\n")
    report_lines.append(
        numeric_summary(wt_valid, wt_numeric_columns).to_markdown(index=False)
    )

    if mutation_file is not None:
        report_lines.append("\n## FoldX mutation screening\n")
        report_lines.append(f"- Mutation source file: `{mutation_file}`")
        report_lines.append(f"- ΔΔG column used: `{ddg_col}`")
        report_lines.append(f"- Mutation rows: {len(mutations)}")

        report_lines.append("\n### Mutation ΔΔG summary\n")
        report_lines.append(
            numeric_summary(mutations, ["foldx_ddg"]).to_markdown(index=False)
        )

        report_lines.append("\n### Mutation effect counts\n")
        report_lines.append(effect_counts.to_markdown(index=False))
    else:
        report_lines.append("\n## FoldX mutation screening\n")
        report_lines.append("- No mutation ΔΔG table was detected automatically.")

    (OUTPUT_DIR / "FOLDX_RESULTS_SUMMARY.md").write_text(
        "\n".join(report_lines)
    )

    print("\nFoldX WT summary:")
    print(numeric_summary(wt_valid, wt_numeric_columns).to_string(index=False))

    if mutation_file is not None:
        print("\nFoldX mutation ΔΔG summary:")
        print(numeric_summary(mutations, ["foldx_ddg"]).to_string(index=False))

        print("\nFoldX mutation effect counts:")
        print(effect_counts.to_string(index=False))
    else:
        print("\nNo FoldX mutation ΔΔG table was detected automatically.")

    print(f"\nOutput directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
