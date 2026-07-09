#!/usr/bin/env python3
# Purpose: Merge experimental metadata.

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


def parse_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def read_table(path: Path, sheet_name: str | None = None) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet_name or 0, dtype=str).fillna("")

    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")

    return pd.read_csv(path, dtype=str, keep_default_na=False).fillna("")


def select_metadata_columns(
    source: pd.DataFrame,
    merge_keys: list[str],
    metadata_columns: list[str],
    metadata_column_prefix: str | None,
) -> list[str]:
    if metadata_columns:
        selected = [c for c in metadata_columns if c in source.columns]
        missing = [c for c in metadata_columns if c not in source.columns]
        if missing:
            raise ValueError(f"Requested metadata columns not found in source: {missing}")
        return selected

    if metadata_column_prefix:
        selected = [c for c in source.columns if c.startswith(metadata_column_prefix)]
        if not selected:
            raise ValueError(f"No metadata columns found with prefix: {metadata_column_prefix}")
        return selected

    excluded = set(merge_keys)
    return [c for c in source.columns if c not in excluded]


def prepare_source_mapping(source: pd.DataFrame, key: str, metadata_cols: list[str]) -> pd.DataFrame:
    if key not in source.columns:
        return pd.DataFrame()

    subset = source[[key] + metadata_cols].copy()
    subset[key] = subset[key].map(normalise_key)
    subset = subset[subset[key].astype(str).str.len() > 0].copy()

    # Keep row with most metadata values when duplicate keys exist
    subset["_metadata_nonempty_count"] = subset[metadata_cols].apply(
        lambda row: sum(is_present(v) for v in row),
        axis=1,
    )

    subset = subset.sort_values("_metadata_nonempty_count", ascending=False)
    subset = subset.drop_duplicates(subset=[key], keep="first")
    subset = subset.drop(columns=["_metadata_nonempty_count"])

    return subset


def merge_by_key(
    master: pd.DataFrame,
    source_map: pd.DataFrame,
    key: str,
    metadata_cols: list[str],
    overwrite: bool,
) -> tuple[pd.DataFrame, int]:
    if source_map.empty or key not in master.columns:
        return master, 0

    master = master.copy()
    master[key] = master[key].map(normalise_key)

    rename_cols = {c: f"{c}__from_metadata_source" for c in metadata_cols}
    source_map_renamed = source_map.rename(columns=rename_cols)

    merged = master.merge(source_map_renamed, on=key, how="left")

    filled_count = 0

    for col in metadata_cols:
        src_col = f"{col}__from_metadata_source"

        if col not in merged.columns:
            merged[col] = ""

        if src_col not in merged.columns:
            continue

        source_present = merged[src_col].map(is_present)

        if overwrite:
            fill_mask = source_present
        else:
            current_missing = ~merged[col].map(is_present)
            fill_mask = current_missing & source_present

        filled_count += int(fill_mask.sum())

        merged.loc[fill_mask, col] = merged.loc[fill_mask, src_col]
        merged = merged.drop(columns=[src_col])

    return merged, filled_count


def metadata_coverage(df: pd.DataFrame, metadata_cols: list[str]) -> dict:
    coverage = {}

    for col in metadata_cols:
        coverage[col] = int(df[col].map(is_present).sum()) if col in df.columns else 0

    coverage["rows_with_any_experimental_metadata"] = int(
        df[metadata_cols].apply(lambda row: any(is_present(v) for v in row), axis=1).sum()
    ) if metadata_cols else 0

    return coverage


def add_row_level_source_label(
    df: pd.DataFrame,
    metadata_cols: list[str],
    source_name: str,
) -> pd.DataFrame:
    df = df.copy()

    has_any = df[metadata_cols].apply(lambda row: any(is_present(v) for v in row), axis=1) if metadata_cols else False

    if "experimental_metadata_source" not in df.columns:
        df["experimental_metadata_source"] = ""

    df.loc[has_any, "experimental_metadata_source"] = source_name

    return df


def summarize(df: pd.DataFrame, metadata_cols: list[str]) -> dict:
    return {
        "rows": int(len(df)),
        "organism_type_counts": dict(Counter(df.get("organism_type", pd.Series(dtype=str)).astype(str))),
        "gh_family_counts": dict(Counter(df.get("gh_family", pd.Series(dtype=str)).astype(str))),
        "metadata_coverage": metadata_coverage(df, metadata_cols),
    }


def write_metadata(output_path: Path, metadata: dict) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge experimental metadata into curated master table.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--master", required=True, help="Input curated master CSV.")
    parser.add_argument("--metadata-source", required=True, help="Experimental metadata source CSV/TSV/XLSX.")
    parser.add_argument("--output", required=True, help="Output merged master CSV.")
    parser.add_argument("--source-name", default="experimental_metadata", help="Name of metadata source, e.g. BRENDA or user_experiment.")
    parser.add_argument("--merge-keys", default="uniprot_accession,primary_id,refseq_acc", help="Comma-separated merge keys.")
    parser.add_argument("--metadata-columns", default=None, help="Comma-separated metadata columns to merge.")
    parser.add_argument("--metadata-column-prefix", default=None, help="Merge all source columns with this prefix.")
    parser.add_argument("--sheet-name", default=None, help="Excel sheet name if using XLSX.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing metadata values in master.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    master_path = resolve_path(args.master, project_root)
    source_path = resolve_path(args.metadata_source, project_root)
    output_path = resolve_path(args.output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "08_merge_experimental_metadata.log"),
        project_root,
    )

    if master_path is None or source_path is None or output_path is None:
        raise ValueError("--master, --metadata-source and --output are required")

    logger = Logger(log_path)

    logger.write("Starting experimental metadata merge")
    logger.write(f"Master: {master_path}")
    logger.write(f"Metadata source: {source_path}")
    logger.write(f"Output: {output_path}")
    logger.write(f"Source name: {args.source_name}")

    if not master_path.exists():
        raise FileNotFoundError(f"Master file not found: {master_path}")

    if not source_path.exists():
        raise FileNotFoundError(f"Metadata source file not found: {source_path}")

    master = read_table(master_path)
    source = read_table(source_path, sheet_name=args.sheet_name)

    merge_keys = parse_list(args.merge_keys)
    metadata_columns = parse_list(args.metadata_columns)
    metadata_cols = select_metadata_columns(
        source=source,
        merge_keys=merge_keys,
        metadata_columns=metadata_columns,
        metadata_column_prefix=args.metadata_column_prefix,
    )

    logger.write(f"Master rows: {len(master)}")
    logger.write(f"Metadata source rows: {len(source)}")
    logger.write(f"Merge keys: {merge_keys}")
    logger.write(f"Metadata columns: {metadata_cols}")

    for col in metadata_cols:
        if col not in master.columns:
            master[col] = ""

    before_summary = summarize(master, metadata_cols)

    fill_by_key = {}
    total_filled = 0

    for key in merge_keys:
        if key in master.columns and key in source.columns:
            source_map = prepare_source_mapping(source, key, metadata_cols)
            master, filled = merge_by_key(
                master=master,
                source_map=source_map,
                key=key,
                metadata_cols=metadata_cols,
                overwrite=args.overwrite,
            )
            fill_by_key[key] = int(filled)
            total_filled += int(filled)
            logger.write(f"Filled metadata cells using {key}: {filled}")
        else:
            fill_by_key[key] = 0
            logger.write(f"Skipped key {key}: not present in both tables")

    master = add_row_level_source_label(master, metadata_cols, args.source_name)

    after_summary = summarize(master, metadata_cols)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "master": str(master_path),
        "metadata_source": str(source_path),
        "output": str(output_path),
        "source_name": args.source_name,
        "merge_keys": merge_keys,
        "metadata_columns": metadata_cols,
        "overwrite": bool(args.overwrite),
        "fill_by_key": fill_by_key,
        "total_filled_cells": int(total_filled),
        "before_summary": before_summary,
        "after_summary": after_summary,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Total filled metadata cells: {total_filled}")
    logger.write(
        f"Rows with any experimental metadata: "
        f"{after_summary['metadata_coverage'].get('rows_with_any_experimental_metadata')}"
    )
    logger.write("Finished experimental metadata merge successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
