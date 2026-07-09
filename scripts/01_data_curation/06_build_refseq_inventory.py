#!/usr/bin/env python3
"""
Build a normalized RefSeq inventory from a curated enzyme master table.

Input:
    Curated master CSV containing at least:
        uniprot_accession, refseq_acc, organism_type, organism, gh_family, sequence

Outputs:
    1. RefSeq inventory CSV:
        uniprot_accession, organism_type, organism, gh_family, refseq_acc

    2. Missing RefSeq CSV:
        rows from the master table without usable RefSeq accessions

    3. Metadata JSON
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import sys
import time

import pandas as pd


REFSEQ_COLUMN_CANDIDATES = [
    "refseq_acc",
    "refseq_ids",
    "xref_refseq",
    "RefSeq",
    "refseq",
]


def load_yaml_config(path: str | None) -> dict:
    if not path:
        return {}

    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml") from exc

    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def deep_get(data: dict, keys: list[str], default=None):
    cur = data
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def resolve_path(path_value: str | None, project_root: str | None = None) -> Path | None:
    if not path_value:
        return None

    path = Path(os.path.expanduser(path_value))

    if path.is_absolute():
        return path

    if project_root:
        return Path(os.path.expanduser(project_root)) / path

    return Path.cwd() / path


class Logger:
    def __init__(self, log_path: Path | None):
        self.log_path = log_path
        if self.log_path:
            self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, msg: str) -> None:
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        if self.log_path:
            with self.log_path.open("a", encoding="utf-8") as handle:
                handle.write(line + "\n")


def clean_text(value) -> str:
    if value is None:
        return ""
    text = str(value)
    text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def normalise_accession(value) -> str:
    return clean_text(value).split("-")[0]


def is_missing(value: str) -> bool:
    text = clean_text(value).lower()
    return text in {"", "nan", "none", "na", "n/a", "null", "unknown"}


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    return None


def normalise_refseq_acc(value: str) -> str:
    value = clean_text(value)
    value = value.strip(";,")

    if is_missing(value):
        return ""

    # RefSeq accessions often look like WP_012345678.1, NP_..., XP_..., YP_...
    # Keep the accession/version but remove trailing punctuation.
    value = re.sub(r"[;,\s]+$", "", value)
    return value


def split_refseq_values(value: str) -> list[str]:
    value = clean_text(value)

    if is_missing(value):
        return []

    # Common UniProt separator is semicolon; allow comma and pipe too.
    parts = re.split(r"[;,|]\s*", value)

    out = []
    for part in parts:
        acc = normalise_refseq_acc(part)
        if acc:
            out.append(acc)

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for acc in out:
        if acc not in seen:
            unique.append(acc)
            seen.add(acc)

    return unique


def ensure_column(df: pd.DataFrame, col: str, default: str = "") -> None:
    if col not in df.columns:
        df[col] = default


def build_refseq_inventory(master: pd.DataFrame, logger: Logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    refseq_col = find_column(master, REFSEQ_COLUMN_CANDIDATES)

    if refseq_col is None:
        logger.write("No RefSeq column found in master table")
        refseq_col = "__missing_refseq__"
        master[refseq_col] = ""

    logger.write(f"Using RefSeq column: {refseq_col}")

    required_defaults = {
        "uniprot_accession": "",
        "organism_type": "",
        "organism": "",
        "gh_family": "",
        "sequence": "",
    }

    for col, default in required_defaults.items():
        ensure_column(master, col, default)

    rows = []
    has_refseq_mask = []

    for _, row in master.iterrows():
        accession = normalise_accession(row.get("uniprot_accession", ""))
        refseq_values = split_refseq_values(row.get(refseq_col, ""))

        has_refseq_mask.append(bool(refseq_values))

        if not accession or not refseq_values:
            continue

        for refseq_acc in refseq_values:
            rows.append({
                "uniprot_accession": accession,
                "organism_type": clean_text(row.get("organism_type", "")),
                "organism": clean_text(row.get("organism", "")),
                "gh_family": clean_text(row.get("gh_family", "")),
                "refseq_acc": refseq_acc,
            })

    inventory_columns = [
        "uniprot_accession",
        "organism_type",
        "organism",
        "gh_family",
        "refseq_acc",
    ]

    inventory = pd.DataFrame(rows, columns=inventory_columns)

    if not inventory.empty:
        inventory = inventory.drop_duplicates(
            subset=["uniprot_accession", "refseq_acc"],
            keep="first",
        ).copy()

    missing = master.loc[[not x for x in has_refseq_mask]].copy()

    return inventory, missing


def summarize_inventory(inv: pd.DataFrame, missing: pd.DataFrame, total_master_rows: int) -> dict:
    summary = {
        "master_rows": int(total_master_rows),
        "inventory_rows": int(len(inv)),
        "missing_refseq_rows": int(len(missing)),
    }

    if not inv.empty:
        summary.update({
            "unique_uniprot_accessions_with_refseq": int(inv["uniprot_accession"].nunique()),
            "unique_refseq_accessions": int(inv["refseq_acc"].nunique()),
            "organism_type_counts": dict(Counter(inv["organism_type"].astype(str))),
            "gh_family_counts": dict(Counter(inv["gh_family"].astype(str))),
        })
    else:
        summary.update({
            "unique_uniprot_accessions_with_refseq": 0,
            "unique_refseq_accessions": 0,
            "organism_type_counts": {},
            "gh_family_counts": {},
        })

    return summary


def write_metadata(output_path: Path, metadata: dict) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build normalized RefSeq inventory.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--input", required=True, help="Input curated master CSV.")
    parser.add_argument("--output", required=True, help="Output RefSeq inventory CSV.")
    parser.add_argument("--missing-output", required=True, help="Output CSV for master rows without RefSeq.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(args.input, project_root)
    output_path = resolve_path(args.output, project_root)
    missing_output_path = resolve_path(args.missing_output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "06_build_refseq_inventory.log"),
        project_root,
    )

    if input_path is None or output_path is None or missing_output_path is None:
        raise ValueError("--input, --output and --missing-output are required")

    logger = Logger(log_path)

    logger.write("Starting RefSeq inventory build")
    logger.write(f"Input master: {input_path}")
    logger.write(f"Output inventory: {output_path}")
    logger.write(f"Missing output: {missing_output_path}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input master file not found: {input_path}")

    master = pd.read_csv(input_path, dtype=str, keep_default_na=False).fillna("")

    inventory, missing = build_refseq_inventory(master, logger)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    missing_output_path.parent.mkdir(parents=True, exist_ok=True)

    inventory.to_csv(output_path, index=False)
    missing.to_csv(missing_output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "input": str(input_path),
        "output": str(output_path),
        "missing_output": str(missing_output_path),
        "summary": summarize_inventory(inventory, missing, len(master)),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Inventory rows: {len(inventory)}")
    logger.write(f"Missing RefSeq rows: {len(missing)}")
    logger.write(f"Unique UniProt with RefSeq: {inventory['uniprot_accession'].nunique() if not inventory.empty else 0}")
    logger.write(f"Unique RefSeq accessions: {inventory['refseq_acc'].nunique() if not inventory.empty else 0}")
    logger.write("Finished RefSeq inventory build successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
