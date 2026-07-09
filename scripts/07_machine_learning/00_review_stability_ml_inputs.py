#!/usr/bin/env python3
"""
Review available input files for the stability ML goal.

This script searches existing data/results folders for files related to:
- experimental Tm/Topt labels
- sequence features
- structural features
- FoldX wild-type stability
- FoldX mutation ΔΔG
- docking
- MD validation
"""

from pathlib import Path

import pandas as pd


BASE_DIR = Path.home() / "xylanase-thesis"

SEARCH_ROOTS = [
    BASE_DIR / "data",
    BASE_DIR / "results",
]

OUTPUT_DIR = BASE_DIR / "results" / "stability_ml_goal" / "input_checks"

FILE_REVIEW = OUTPUT_DIR / "goal_input_file_candidates.csv"
COLUMN_REVIEW = OUTPUT_DIR / "goal_input_column_review.csv"
REPORT_FILE = OUTPUT_DIR / "goal_input_review_report.md"


KEYWORDS = [
    "tmxyl",
    "toptxyl",
    "temperature",
    "thermal",
    "stability",
    "foldx",
    "ddg",
    "mutation",
    "mutant",
    "docking",
    "vina",
    "xylobiose",
    "xylotriose",
    "md",
    "rmsd",
    "rmsf",
    "sasa",
    "structure",
    "structural",
    "modeller",
    "tmalign",
    "brenda",
    "sequence_features",
]


TABLE_EXTENSIONS = {".csv", ".tsv", ".txt", ".xlsx"}


def file_matches(path):
    """Check whether a file name/path is relevant."""
    text = str(path).lower()
    return any(keyword in text for keyword in KEYWORDS)


def read_table_preview(path):
    """Read a small preview of a table file."""
    suffix = path.suffix.lower()

    try:
        if suffix == ".csv":
            return pd.read_csv(path, nrows=5)
        if suffix == ".tsv":
            return pd.read_csv(path, sep="\t", nrows=5)
        if suffix == ".xlsx":
            return pd.read_excel(path, nrows=5)
        if suffix == ".txt":
            try:
                return pd.read_csv(path, sep="\t", nrows=5)
            except Exception:
                return pd.read_csv(path, nrows=5)
    except Exception:
        return None

    return None


def count_rows(path):
    """Count rows for readable table files."""
    suffix = path.suffix.lower()

    try:
        if suffix == ".csv":
            return len(pd.read_csv(path))
        if suffix == ".tsv":
            return len(pd.read_csv(path, sep="\t"))
        if suffix == ".xlsx":
            return len(pd.read_excel(path))
        if suffix == ".txt":
            try:
                return len(pd.read_csv(path, sep="\t"))
            except Exception:
                return len(pd.read_csv(path))
    except Exception:
        return None

    return None


def classify_file(path):
    """Assign a simple file category."""
    text = str(path).lower()

    if "tmxyl" in text or "toptxyl" in text or "experimental" in text:
        return "experimental_stability"
    if "foldx" in text and ("ddg" in text or "mutation" in text or "mutant" in text):
        return "mutation_foldx"
    if "foldx" in text and "stability" in text:
        return "wildtype_foldx"
    if "docking" in text or "vina" in text:
        return "docking_activity_proxy"
    if "md" in text or "rmsd" in text or "rmsf" in text:
        return "md_validation"
    if "structure" in text or "structural" in text or "modeller" in text or "tmalign" in text:
        return "structure_features"
    if "sequence" in text:
        return "sequence_features"

    return "other_relevant"


def main():
    """Run input review."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    file_rows = []
    column_rows = []

    for root in SEARCH_ROOTS:
        if not root.exists():
            continue

        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue

            if not file_matches(path):
                continue

            relative_path = path.relative_to(BASE_DIR)
            suffix = path.suffix.lower()
            size_mb = path.stat().st_size / (1024 * 1024)

            row_count = None
            columns = []

            if suffix in TABLE_EXTENSIONS:
                row_count = count_rows(path)
                preview = read_table_preview(path)

                if preview is not None:
                    columns = preview.columns.tolist()

                    column_rows.append(
                        {
                            "path": str(relative_path),
                            "category": classify_file(path),
                            "row_count": row_count,
                            "column_count": len(columns),
                            "columns": "; ".join(columns),
                        }
                    )

            file_rows.append(
                {
                    "path": str(relative_path),
                    "category": classify_file(path),
                    "suffix": suffix,
                    "size_mb": round(size_mb, 4),
                    "row_count": row_count,
                }
            )

    file_data = pd.DataFrame(file_rows)
    column_data = pd.DataFrame(column_rows)

    file_data.to_csv(FILE_REVIEW, index=False)
    column_data.to_csv(COLUMN_REVIEW, index=False)

    report_lines = [
        "# Stability ML Goal Input Review",
        "",
        f"Candidate files found: {len(file_data)}",
        f"Readable table files: {len(column_data)}",
        "",
        "## File categories",
        "",
        file_data["category"].value_counts(dropna=False).to_string()
        if len(file_data)
        else "No files found.",
        "",
        "## Top readable tables by row count",
        "",
    ]

    if len(column_data):
        report_lines.append(
            column_data.sort_values("row_count", ascending=False)
            .head(30)
            [["category", "row_count", "column_count", "path"]]
            .to_string(index=False)
        )

    REPORT_FILE.write_text("\n".join(report_lines) + "\n")

    print(f"Wrote: {FILE_REVIEW}")
    print(f"Wrote: {COLUMN_REVIEW}")
    print(f"Wrote: {REPORT_FILE}")
    print()
    print("File categories:")
    print(file_data["category"].value_counts(dropna=False))
    print()
    print("Top readable tables:")
    if len(column_data):
        print(
            column_data.sort_values("row_count", ascending=False)
            .head(20)
            [["category", "row_count", "column_count", "path"]]
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
