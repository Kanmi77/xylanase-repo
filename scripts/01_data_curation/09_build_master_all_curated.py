#!/usr/bin/env python3
"""
Build final curated enzyme master table from an enriched UniProt/CAZy table.

For the xylanase thesis, this creates the GH10/GH11 curated master table.
For other enzyme projects, accepted families and organism groups can be changed
through the config file or command-line options.

Example:
    python scripts/01_curation/09_build_master_all_curated.py \
        --config config/xylanase_config.yaml \
        --input data/curated/uniprot_xylanase_curated_with_cazy.csv \
        --output data/curated/xylanase_master_all_curated_universal.csv
"""

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
    return str(value).strip()


def parse_list(value, default=None) -> list[str]:
    if value is None:
        return default or []

    if isinstance(value, list):
        return [str(x).strip() for x in value if str(x).strip()]

    if isinstance(value, str):
        return [x.strip() for x in value.split(",") if x.strip()]

    return default or []


def ensure_column(df: pd.DataFrame, col: str, default: str = "") -> None:
    if col not in df.columns:
        df[col] = default


def build_final_master(
    df: pd.DataFrame,
    accepted_families: list[str],
    accepted_organism_groups: list[str],
    keep_unknown_family: bool = False,
    keep_unknown_organism: bool = False,
) -> pd.DataFrame:
    data = df.copy()

    required_defaults = {
        "source": "uniprot",
        "primary_id": "",
        "uniprot_accession": "",
        "refseq_acc": "",
        "organism_type": "unknown",
        "organism": "",
        "xref_cazy": "",
        "has_cazy_xref": False,
        "cazy_status": "",
        "gh_family": "unknown",
        "sequence": "",
    }

    for col, default in required_defaults.items():
        ensure_column(data, col, default)

    data["uniprot_accession"] = data["uniprot_accession"].fillna("").astype(str).str.strip()
    data["primary_id"] = data["primary_id"].fillna("").astype(str).str.strip()
    data["sequence"] = data["sequence"].fillna("").astype(str).str.replace(" ", "", regex=False)
    data["gh_family"] = data["gh_family"].fillna("unknown").astype(str).str.strip()
    data["organism_type"] = data["organism_type"].fillna("unknown").astype(str).str.strip()

    # Remove rows without accession or sequence
    data = data[data["uniprot_accession"].str.len() > 0].copy()
    data = data[data["sequence"].str.len() > 0].copy()

    # Deduplicate by accession
    data = data.drop_duplicates(subset=["uniprot_accession"], keep="first").copy()

    if accepted_families and not keep_unknown_family:
        data = data[data["gh_family"].isin(accepted_families)].copy()

    if accepted_organism_groups and not keep_unknown_organism:
        data = data[data["organism_type"].isin(accepted_organism_groups)].copy()

    # Standard final thesis columns first
    final_columns = [
        "source",
        "primary_id",
        "uniprot_accession",
        "refseq_acc",
        "organism_type",
        "organism",
        "xref_cazy",
        "has_cazy_xref",
        "cazy_status",
        "gh_family",
        "sequence",
    ]

    extra_columns = [c for c in data.columns if c not in final_columns]
    data = data[final_columns + extra_columns].copy()

    return data


def summarize(df: pd.DataFrame) -> dict:
    summary = {
        "rows": int(len(df)),
        "organism_type_counts": dict(Counter(df.get("organism_type", pd.Series(dtype=str)).astype(str))),
        "gh_family_counts": dict(Counter(df.get("gh_family", pd.Series(dtype=str)).astype(str))),
    }

    if "has_cazy_xref" in df.columns:
        summary["has_cazy_xref_counts"] = dict(Counter(df["has_cazy_xref"].astype(str)))

    if "has_pdb" in df.columns:
        summary["has_pdb_counts"] = dict(Counter(df["has_pdb"].astype(str)))

    if "has_refseq" in df.columns:
        summary["has_refseq_counts"] = dict(Counter(df["has_refseq"].astype(str)))

    return summary


def write_metadata(output_path: Path, metadata: dict) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build final curated enzyme master table.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--input", required=True, help="Input enriched curated CSV.")
    parser.add_argument("--output", required=True, help="Output final master CSV.")
    parser.add_argument("--log", default=None, help="Log file path.")
    parser.add_argument("--families", default=None, help="Comma-separated accepted families, e.g. GH10,GH11.")
    parser.add_argument("--organism-groups", default=None, help="Comma-separated accepted organism groups.")
    parser.add_argument("--keep-unknown-family", action="store_true", help="Keep rows with unknown family.")
    parser.add_argument("--keep-unknown-organism", action="store_true", help="Keep rows with unknown organism group.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(args.input, project_root)
    output_path = resolve_path(args.output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "09_build_master_all_curated.log"),
        project_root,
    )

    if input_path is None or output_path is None:
        raise ValueError("--input and --output are required")

    logger = Logger(log_path)

    logger.write("Starting final master curation")
    logger.write(f"Input: {input_path}")
    logger.write(f"Output: {output_path}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(input_path, dtype=str, keep_default_na=False).fillna("")

    accepted_families = parse_list(
        args.families,
        default=parse_list(deep_get(config, ["dataset", "families"], None)),
    )

    accepted_organism_groups = parse_list(
        args.organism_groups,
        default=parse_list(deep_get(config, ["curation", "accepted_organism_groups"], None)),
    )

    logger.write(f"Input rows: {len(df)}")
    logger.write(f"Accepted families: {accepted_families}")
    logger.write(f"Accepted organism groups: {accepted_organism_groups}")

    before_summary = summarize(df)

    final = build_final_master(
        df=df,
        accepted_families=accepted_families,
        accepted_organism_groups=accepted_organism_groups,
        keep_unknown_family=args.keep_unknown_family,
        keep_unknown_organism=args.keep_unknown_organism,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(output_path, index=False)

    after_summary = summarize(final)

    metadata = {
        "script": Path(__file__).name,
        "input": str(input_path),
        "output": str(output_path),
        "accepted_families": accepted_families,
        "accepted_organism_groups": accepted_organism_groups,
        "keep_unknown_family": args.keep_unknown_family,
        "keep_unknown_organism": args.keep_unknown_organism,
        "before_summary": before_summary,
        "after_summary": after_summary,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Final rows: {len(final)}")
    logger.write(f"Final organism counts: {after_summary['organism_type_counts']}")
    logger.write(f"Final GH-family counts: {after_summary['gh_family_counts']}")
    logger.write("Finished final master curation successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
