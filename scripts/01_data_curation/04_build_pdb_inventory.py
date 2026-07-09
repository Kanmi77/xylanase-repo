#!/usr/bin/env python3
"""
Build or filter a normalized PDB inventory.

This script is generic and supports two modes:

1. Build mode:
   Extract PDB IDs from a master table containing columns such as:
   xref_pdb, pdb_ids, PDB.

2. Filter mode:
   If the master table does not contain PDB columns, provide an existing
   PDB inventory and filter it to the accessions present in the master table.

Standard output columns:
    uniprot_accession, organism_type, organism, gh_family, pdb_id
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


PDB_COLUMN_CANDIDATES = [
    "pdb_ids",
    "xref_pdb",
    "PDB",
    "pdb",
    "pdb_id",
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


def normalise_pdb_id(value) -> str:
    value = clean_text(value).upper()
    value = value.replace(" ", "")
    value = re.sub(r"[^A-Z0-9]", "", value)
    if len(value) == 4:
        return value
    return ""


def split_pdb_ids(value: str) -> list[str]:
    value = clean_text(value)
    if not value:
        return []

    parts = re.split(r"[;,|]\s*", value)
    ids = []
    for part in parts:
        pdb_id = normalise_pdb_id(part)
        if pdb_id:
            ids.append(pdb_id)

    return ids


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    return None


def ensure_column(df: pd.DataFrame, col: str, default: str = "") -> None:
    if col not in df.columns:
        df[col] = default


def build_inventory_from_master(master: pd.DataFrame, logger: Logger) -> pd.DataFrame:
    """Extract PDB IDs directly from a master table."""
    pdb_col = find_column(master, PDB_COLUMN_CANDIDATES)

    if pdb_col is None:
        logger.write("No PDB column found in master table")
        return pd.DataFrame()

    logger.write(f"Using PDB column from master: {pdb_col}")

    ensure_column(master, "uniprot_accession", "")
    ensure_column(master, "organism_type", "")
    ensure_column(master, "organism", "")
    ensure_column(master, "gh_family", "")

    rows = []

    for _, row in master.iterrows():
        accession = normalise_accession(row.get("uniprot_accession", ""))
        if not accession:
            continue

        pdb_ids = split_pdb_ids(row.get(pdb_col, ""))

        for pdb_id in pdb_ids:
            rows.append({
                "uniprot_accession": accession,
                "organism_type": clean_text(row.get("organism_type", "")),
                "organism": clean_text(row.get("organism", "")),
                "gh_family": clean_text(row.get("gh_family", "")),
                "pdb_id": pdb_id,
            })

    inv = pd.DataFrame(rows)
    if inv.empty:
        return inv

    inv = inv.drop_duplicates(subset=["uniprot_accession", "pdb_id"], keep="first")
    return inv


def filter_existing_inventory(
    master: pd.DataFrame,
    existing_inventory: pd.DataFrame,
    logger: Logger,
) -> pd.DataFrame:
    """Filter an existing PDB inventory to accessions in the master table."""
    if "uniprot_accession" not in master.columns:
        raise ValueError("Master table must contain uniprot_accession")

    if "uniprot_accession" not in existing_inventory.columns:
        raise ValueError("Existing PDB inventory must contain uniprot_accession")

    if "pdb_id" not in existing_inventory.columns:
        raise ValueError("Existing PDB inventory must contain pdb_id")

    master = master.copy()
    existing_inventory = existing_inventory.copy()

    master["uniprot_accession"] = master["uniprot_accession"].map(normalise_accession)
    existing_inventory["uniprot_accession"] = existing_inventory["uniprot_accession"].map(normalise_accession)
    existing_inventory["pdb_id"] = existing_inventory["pdb_id"].map(normalise_pdb_id)

    existing_inventory = existing_inventory[existing_inventory["pdb_id"].str.len() == 4].copy()

    keep_ids = set(master["uniprot_accession"])
    filtered = existing_inventory[existing_inventory["uniprot_accession"].isin(keep_ids)].copy()

    # Merge missing metadata from master
    meta_cols = ["uniprot_accession"]
    for c in ["organism_type", "organism", "gh_family"]:
        if c in master.columns:
            meta_cols.append(c)

    meta = master[meta_cols].drop_duplicates(subset=["uniprot_accession"], keep="first")

    # Avoid duplicate organism columns if existing inventory already has them
    drop_cols = [c for c in ["organism_type", "organism", "gh_family"] if c in filtered.columns]
    filtered = filtered.drop(columns=drop_cols, errors="ignore")
    filtered = filtered.merge(meta, on="uniprot_accession", how="left")

    for c in ["organism_type", "organism", "gh_family"]:
        ensure_column(filtered, c, "")

    filtered = filtered[
        ["uniprot_accession", "organism_type", "organism", "gh_family", "pdb_id"]
    ].copy()

    filtered = filtered.drop_duplicates(subset=["uniprot_accession", "pdb_id"], keep="first")

    logger.write(f"Filtered existing inventory from {len(existing_inventory)} rows to {len(filtered)} rows")
    return filtered


def summarize(inv: pd.DataFrame) -> dict:
    if inv.empty:
        return {
            "rows": 0,
            "unique_uniprot_accessions": 0,
            "unique_pdb_ids": 0,
        }

    return {
        "rows": int(len(inv)),
        "unique_uniprot_accessions": int(inv["uniprot_accession"].nunique()),
        "unique_pdb_ids": int(inv["pdb_id"].nunique()),
        "organism_type_counts": dict(Counter(inv["organism_type"].astype(str))),
        "gh_family_counts": dict(Counter(inv["gh_family"].astype(str))),
    }


def write_metadata(output_path: Path, metadata: dict) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build or filter normalized PDB inventory.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--input", required=True, help="Input master CSV.")
    parser.add_argument("--output", required=True, help="Output PDB inventory CSV.")
    parser.add_argument("--existing-inventory", default=None, help="Existing PDB inventory CSV to filter if master lacks PDB columns.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(args.input, project_root)
    output_path = resolve_path(args.output, project_root)
    existing_inventory_path = resolve_path(args.existing_inventory, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "04_build_pdb_inventory.log"),
        project_root,
    )

    if input_path is None or output_path is None:
        raise ValueError("--input and --output are required")

    logger = Logger(log_path)

    logger.write("Starting PDB inventory build/filter")
    logger.write(f"Input master: {input_path}")
    logger.write(f"Output inventory: {output_path}")
    logger.write(f"Existing inventory: {existing_inventory_path}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input master file not found: {input_path}")

    master = pd.read_csv(input_path, dtype=str, keep_default_na=False).fillna("")

    inv = build_inventory_from_master(master, logger)

    mode = "built_from_master"

    if inv.empty:
        if existing_inventory_path is None or not existing_inventory_path.exists():
            raise ValueError(
                "No PDB IDs found in master table and no valid --existing-inventory was supplied."
            )

        existing = pd.read_csv(existing_inventory_path, dtype=str, keep_default_na=False).fillna("")
        inv = filter_existing_inventory(master, existing, logger)
        mode = "filtered_existing_inventory"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    inv.to_csv(output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "mode": mode,
        "input": str(input_path),
        "existing_inventory": str(existing_inventory_path) if existing_inventory_path else None,
        "output": str(output_path),
        "summary": summarize(inv),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Mode: {mode}")
    logger.write(f"Rows: {len(inv)}")
    logger.write(f"Unique accessions: {inv['uniprot_accession'].nunique() if not inv.empty else 0}")
    logger.write(f"Unique PDB IDs: {inv['pdb_id'].nunique() if not inv.empty else 0}")
    logger.write("Finished PDB inventory build/filter successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
