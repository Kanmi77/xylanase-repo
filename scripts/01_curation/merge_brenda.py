#!/usr/bin/env python3
# Purpose: Merge BRENDA metadata.

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import sys
import time

import pandas as pd


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
    text = str(value).strip()
    if text.lower() in {"nan", "none", "na", "n/a", "null"}:
        return ""
    return text


def normalise_key(value) -> str:
    return clean_text(value).split("-")[0]


def is_present(value) -> bool:
    return bool(clean_text(value))


def get_brenda_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("brenda_")]


def prepare_source_mapping(source: pd.DataFrame, key: str, brenda_cols: list[str]) -> pd.DataFrame:
    if key not in source.columns:
        return pd.DataFrame()

    subset = source[[key] + brenda_cols].copy()
    subset[key] = subset[key].map(normalise_key)
    subset = subset[subset[key].astype(str).str.len() > 0].copy()

    # Keep the row with the most BRENDA information for duplicate keys
    subset["_brenda_nonempty_count"] = subset[brenda_cols].apply(
        lambda row: sum(is_present(v) for v in row),
        axis=1,
    )
    subset = subset.sort_values("_brenda_nonempty_count", ascending=False)
    subset = subset.drop_duplicates(subset=[key], keep="first")
    subset = subset.drop(columns=["_brenda_nonempty_count"])

    return subset


def fill_from_source(master: pd.DataFrame, source_map: pd.DataFrame, key: str, brenda_cols: list[str]) -> tuple[pd.DataFrame, int]:
    if source_map.empty or key not in master.columns:
        return master, 0

    master = master.copy()
    master[key] = master[key].map(normalise_key)

    rename_cols = {c: f"{c}__from_source" for c in brenda_cols}
    source_map_renamed = source_map.rename(columns=rename_cols)

    merged = master.merge(source_map_renamed, on=key, how="left")

    filled_count = 0

    for col in brenda_cols:
        src_col = f"{col}__from_source"

        if col not in merged.columns:
            merged[col] = ""

        if src_col not in merged.columns:
            continue

        before_missing = ~merged[col].map(is_present)
        source_present = merged[src_col].map(is_present)
        fill_mask = before_missing & source_present

        filled_count += int(fill_mask.sum())

        merged.loc[fill_mask, col] = merged.loc[fill_mask, src_col]
        merged = merged.drop(columns=[src_col])

    return merged, filled_count


def brenda_coverage(df: pd.DataFrame, brenda_cols: list[str]) -> dict:
    coverage = {}

    for col in brenda_cols:
        coverage[col] = int(df[col].map(is_present).sum()) if col in df.columns else 0

    main_fields = [
        "brenda_temperature_optimum",
        "brenda_temperature_range",
        "brenda_temperature_stability",
        "brenda_ph_optimum",
        "brenda_ph_range",
    ]

    coverage["rows_with_any_brenda"] = int(
        df[brenda_cols].apply(lambda row: any(is_present(v) for v in row), axis=1).sum()
    ) if brenda_cols else 0

    for field in main_fields:
        if field in df.columns:
            coverage[f"rows_with_{field}"] = int(df[field].map(is_present).sum())

    return coverage


def summarize(df: pd.DataFrame, brenda_cols: list[str]) -> dict:
    summary = {
        "rows": int(len(df)),
        "organism_type_counts": dict(Counter(df.get("organism_type", pd.Series(dtype=str)).astype(str))),
        "gh_family_counts": dict(Counter(df.get("gh_family", pd.Series(dtype=str)).astype(str))),
        "brenda_coverage": brenda_coverage(df, brenda_cols),
    }
    return summary


def write_metadata(output_path: Path, metadata: dict) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge BRENDA metadata into curated master table.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--master", required=True, help="Input curated master CSV.")
    parser.add_argument("--brenda-source", required=True, help="BRENDA-enriched source CSV.")
    parser.add_argument("--output", required=True, help="Output BRENDA-enriched master CSV.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    master_path = resolve_path(args.master, project_root)
    source_path = resolve_path(args.brenda_source, project_root)
    output_path = resolve_path(args.output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "08_merge_brenda_into_master.log"),
        project_root,
    )

    if master_path is None or source_path is None or output_path is None:
        raise ValueError("--master, --brenda-source and --output are required")

    logger = Logger(log_path)

    logger.write("Starting BRENDA merge")
    logger.write(f"Master: {master_path}")
    logger.write(f"BRENDA source: {source_path}")
    logger.write(f"Output: {output_path}")

    if not master_path.exists():
        raise FileNotFoundError(f"Master file not found: {master_path}")

    if not source_path.exists():
        raise FileNotFoundError(f"BRENDA source file not found: {source_path}")

    master = pd.read_csv(master_path, dtype=str, keep_default_na=False).fillna("")
    source = pd.read_csv(source_path, dtype=str, keep_default_na=False).fillna("")

    brenda_cols = get_brenda_columns(source)

    if not brenda_cols:
        raise ValueError(f"No BRENDA columns found in source file: {source_path}")

    logger.write(f"Master rows: {len(master)}")
    logger.write(f"Source rows: {len(source)}")
    logger.write(f"BRENDA columns: {brenda_cols}")

    # Ensure all BRENDA columns exist in master
    for col in brenda_cols:
        if col not in master.columns:
            master[col] = ""

    before_summary = summarize(master, brenda_cols)

    total_filled = 0
    fill_by_key = {}

    for key in ["uniprot_accession", "primary_id", "refseq_acc"]:
        if key in master.columns and key in source.columns:
            source_map = prepare_source_mapping(source, key, brenda_cols)
            master, filled = fill_from_source(master, source_map, key, brenda_cols)
            fill_by_key[key] = filled
            total_filled += filled
            logger.write(f"Filled BRENDA cells using {key}: {filled}")

    after_summary = summarize(master, brenda_cols)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "master": str(master_path),
        "brenda_source": str(source_path),
        "output": str(output_path),
        "brenda_columns": brenda_cols,
        "fill_by_key": fill_by_key,
        "total_filled_cells": int(total_filled),
        "before_summary": before_summary,
        "after_summary": after_summary,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Total filled BRENDA cells: {total_filled}")
    logger.write(f"Rows with any BRENDA after merge: {after_summary['brenda_coverage'].get('rows_with_any_brenda')}")
    logger.write("Finished BRENDA merge successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
