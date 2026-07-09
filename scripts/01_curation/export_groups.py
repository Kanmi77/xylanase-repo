#!/usr/bin/env python3
# Purpose: Export GH10/GH11 group files.

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import sys
import textwrap
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


def parse_list(value, default=None) -> list[str]:
    if value is None:
        return default or []

    if isinstance(value, list):
        return [str(x).strip() for x in value if str(x).strip()]

    if isinstance(value, str):
        return [x.strip() for x in value.split(",") if x.strip()]

    return default or []


def clean_text(value) -> str:
    if value is None:
        return ""
    text = str(value)
    text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def safe_filename(value: str) -> str:
    value = clean_text(value)
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def normalise_sequence(seq: str) -> str:
    seq = clean_text(seq).replace(" ", "")
    seq = re.sub(r"[^A-Za-z*.-]", "", seq)
    return seq.upper()


def make_fasta_header(row: pd.Series) -> str:
    accession = clean_text(row.get("uniprot_accession", "")) or clean_text(row.get("primary_id", "unknown"))
    family = clean_text(row.get("gh_family", "unknown"))
    organism_type = clean_text(row.get("organism_type", "unknown"))
    organism = clean_text(row.get("organism", ""))
    protein = clean_text(row.get("protein_name", ""))

    parts = [
        accession,
        f"family={family}",
        f"organism_type={organism_type}",
    ]

    if organism:
        parts.append(f"organism={organism}")

    if protein:
        parts.append(f"protein={protein}")

    return " | ".join(parts)


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for _, row in df.iterrows():
            seq = normalise_sequence(row.get("sequence", ""))
            if not seq:
                continue

            header = make_fasta_header(row)
            handle.write(f">{header}\n")
            handle.write("\n".join(textwrap.wrap(seq, width=80)))
            handle.write("\n")


def summarize_group(df: pd.DataFrame, organism_group: str, family: str, csv_path: Path, fasta_path: Path, txt_path: Path) -> dict:
    return {
        "organism_type": organism_group,
        "gh_family": family,
        "rows": int(len(df)),
        "csv": str(csv_path),
        "fasta": str(fasta_path),
        "txt": str(txt_path),
        "has_pdb_true": int((df.get("has_pdb", pd.Series(dtype=str)).astype(str).str.lower() == "true").sum()) if "has_pdb" in df.columns else None,
        "has_refseq_true": int((df.get("has_refseq", pd.Series(dtype=str)).astype(str).str.lower() == "true").sum()) if "has_refseq" in df.columns else None,
    }


def write_metadata(path: Path, metadata: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export grouped GH-family enzyme datasets.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--input", required=True, help="Input curated master CSV.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--log", default=None, help="Log file path.")
    parser.add_argument("--families", default="GH10,GH11", help="Comma-separated families.")
    parser.add_argument("--organism-groups", default=None, help="Comma-separated organism groups.")
    parser.add_argument("--prefix", default="", help="Optional filename prefix.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(args.input, project_root)
    output_dir = resolve_path(args.output_dir, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "11_export_gh10_gh11_groups.log"),
        project_root,
    )

    if input_path is None or output_dir is None:
        raise ValueError("--input and --output-dir are required")

    logger = Logger(log_path)

    logger.write("Starting grouped GH-family export")
    logger.write(f"Input: {input_path}")
    logger.write(f"Output directory: {output_dir}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(input_path, dtype=str, keep_default_na=False).fillna("")

    required_cols = ["organism_type", "gh_family", "sequence"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")

    families = parse_list(args.families)

    if args.organism_groups:
        organism_groups = parse_list(args.organism_groups)
    else:
        organism_groups = sorted(df["organism_type"].dropna().astype(str).unique().tolist())

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.write(f"Rows loaded: {len(df)}")
    logger.write(f"Families: {families}")
    logger.write(f"Organism groups: {organism_groups}")

    summary_rows = []

    for organism_group in organism_groups:
        for family in families:
            subset = df[
                (df["organism_type"].astype(str) == organism_group)
                & (df["gh_family"].astype(str) == family)
            ].copy()

            base_name = safe_filename(f"{args.prefix}{organism_group}_{family}")
            csv_path = output_dir / f"{base_name}.csv"
            fasta_path = output_dir / f"{base_name}.fasta"
            txt_path = output_dir / f"{base_name}.txt"

            subset.to_csv(csv_path, index=False)
            write_fasta(subset, fasta_path)
            write_fasta(subset, txt_path)

            summary = summarize_group(subset, organism_group, family, csv_path, fasta_path, txt_path)
            summary_rows.append(summary)

            logger.write(f"Exported {organism_group} {family}: {len(subset)} rows")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "group_export_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "input": str(input_path),
        "output_dir": str(output_dir),
        "input_rows": int(len(df)),
        "families": families,
        "organism_groups": organism_groups,
        "overall_organism_counts": dict(Counter(df["organism_type"].astype(str))),
        "overall_family_counts": dict(Counter(df["gh_family"].astype(str))),
        "summary_csv": str(summary_path),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_dir / "group_export_metadata.json"
    write_metadata(metadata_path, metadata)

    logger.write(f"Saved summary: {summary_path}")
    logger.write(f"Saved metadata: {metadata_path}")
    logger.write("Finished grouped GH-family export successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
