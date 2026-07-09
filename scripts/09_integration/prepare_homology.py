#!/usr/bin/env python3
# Purpose: Prepare homology docking candidates.

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "results"
OUT_DIR = RESULTS_DIR / "integration"
OUT_CSV = OUT_DIR / "homology_candidate_level.csv"
OUT_REPORT = OUT_DIR / "homology_candidate_level_report.md"


def find_column(columns, options):
    for option in options:
        if option in columns:
            return option
    return None


def find_input_file():
    candidates = []

    for path in RESULTS_DIR.rglob("*.csv"):
        name = path.name.lower()

        if "homology" not in name or "docking" not in name:
            continue

        if "candidate_level" in name:
            continue

        candidates.append(path)

    for path in sorted(candidates):
        try:
            columns = set(pd.read_csv(path, nrows=5).columns)
        except Exception:
            continue

        has_identity = bool(
            {"candidate_id", "uniprot_accession", "protein"} & columns
        )
        has_delta = any("delta" in column.lower() for column in columns)

        if has_identity and has_delta:
            return path

    raise SystemExit("No suitable homology docking CSV found.")


def first_value(series):
    values = series.dropna()

    if values.empty:
        return ""

    return values.iloc[0]


def collect_ligands(frame):
    ligand_col = find_column(frame.columns, ["ligand", "ligand_name"])

    if ligand_col is None:
        return ""

    values = frame[ligand_col].dropna().astype(str).unique()
    return ";".join(sorted(values))


def summarise_numeric(frame, column):
    values = pd.to_numeric(frame[column], errors="coerce").dropna()

    if values.empty:
        return {
            f"{column}_mean": "",
            f"{column}_best": "",
            f"{column}_worst": "",
        }

    return {
        f"{column}_mean": values.mean(),
        f"{column}_best": values.min(),
        f"{column}_worst": values.max(),
    }


def build_candidate_table(frame):
    mutation_col = find_column(
        frame.columns,
        ["mutation", "foldx_mutation_code", "mutation_code"],
    )
    candidate_col = find_column(frame.columns, ["candidate_id", "protein"])

    group_cols = [
        column
        for column in [
            candidate_col,
            "protein",
            "uniprot_accession",
            mutation_col,
        ]
        if column is not None and column in frame.columns
    ]

    group_cols = list(dict.fromkeys(group_cols))

    if not group_cols:
        raise SystemExit("No candidate grouping columns found.")

    delta_col = find_column(
        frame.columns,
        [
            "delta_binding_mut_minus_wt",
            "mean_docking_change",
            "delta_affinity_mut_minus_wt",
        ],
    )

    top3_col = find_column(
        frame.columns,
        [
            "delta_top3_mut_minus_wt",
            "delta_top3_binding_mut_minus_wt",
        ],
    )

    foldx_col = find_column(
        frame.columns,
        ["foldx_ddg", "foldx_ddg_first", "ddg", "foldx_delta_delta_g"],
    )

    rows = []

    for keys, group in frame.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)

        row = dict(zip(group_cols, keys))
        row["n_ligand_rows"] = len(group)
        row["ligands"] = collect_ligands(group)

        if foldx_col is not None:
            row["foldx_ddg"] = first_value(
                pd.to_numeric(group[foldx_col], errors="coerce")
            )

        if delta_col is not None:
            row.update(summarise_numeric(group, delta_col))
            values = pd.to_numeric(group[delta_col], errors="coerce").dropna()
            row["docking_improved_any_ligand"] = bool((values <= 0).any())
            row["docking_improved_all_ligands"] = bool((values <= 0).all())

        if top3_col is not None:
            row.update(summarise_numeric(group, top3_col))

        rows.append(row)

    return pd.DataFrame(rows)


def write_report(input_file, table):
    lines = [
        "# Homology docking candidate summary",
        "",
        f"Input file: `{input_file}`",
        f"Output file: `{OUT_CSV}`",
        f"Candidate rows: {len(table)}",
        "",
    ]

    OUT_REPORT.write_text("\n".join(lines))


def main():
    input_file = find_input_file()
    frame = pd.read_csv(input_file)
    table = build_candidate_table(frame)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_CSV, index=False)
    write_report(input_file, table)

    print(f"Saved: {OUT_CSV}")
    print(f"Rows: {len(table)}")


if __name__ == "__main__":
    main()
