#!/usr/bin/env python3
"""
Compute sequence-derived ProtParam features from a sequence inventory.

Input:
    Sequence inventory CSV from:
        01_prepare_sequence_phylogeny_inputs.py

Required columns:
    uniprot_accession, organism_type, gh_family, clean_sequence

Outputs:
    - per-protein ProtParam feature table
    - group summary table
    - failed sequence table
    - metadata JSON
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


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


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


def standardize_sequence(seq: str) -> tuple[str, int, int]:
    seq = clean_text(seq).upper()
    seq = re.sub(r"\s+", "", seq)

    standard = []
    nonstandard_count = 0
    removed_count = 0

    for aa in seq:
        if aa in STANDARD_AA:
            standard.append(aa)
        elif aa.isalpha():
            nonstandard_count += 1
        else:
            removed_count += 1

    return "".join(standard), nonstandard_count, removed_count


def aa_composition(seq: str) -> dict:
    length = len(seq)
    counts = Counter(seq)

    out = {}

    for aa in sorted(STANDARD_AA):
        out[f"aa_count_{aa}"] = int(counts.get(aa, 0))
        out[f"aa_fraction_{aa}"] = float(counts.get(aa, 0) / length) if length else 0.0

    return out


def compute_features_for_sequence(seq: str) -> dict:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise ImportError("Biopython is required. Install with: pip install biopython") from exc

    analysis = ProteinAnalysis(seq)

    helix, turn, sheet = analysis.secondary_structure_fraction()

    features = {
        "sequence_length_protparam": int(len(seq)),
        "molecular_weight": float(analysis.molecular_weight()),
        "aromaticity": float(analysis.aromaticity()),
        "instability_index": float(analysis.instability_index()),
        "isoelectric_point": float(analysis.isoelectric_point()),
        "gravy": float(analysis.gravy()),
        "helix_fraction": float(helix),
        "turn_fraction": float(turn),
        "sheet_fraction": float(sheet),
    }

    features.update(aa_composition(seq))

    return features


def compute_features(df: pd.DataFrame, min_standard_length: int, logger: Logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = ["uniprot_accession", "organism_type", "gh_family"]

    for col in required:
        if col not in df.columns:
            raise ValueError(f"Input inventory is missing required column: {col}")

    sequence_col = "clean_sequence" if "clean_sequence" in df.columns else "sequence"

    if sequence_col not in df.columns:
        raise ValueError("Input inventory must contain either clean_sequence or sequence column")

    feature_rows = []
    failed_rows = []

    for _, row in df.iterrows():
        accession = clean_text(row.get("uniprot_accession", ""))
        raw_seq = clean_text(row.get(sequence_col, ""))

        standard_seq, nonstandard_count, removed_count = standardize_sequence(raw_seq)

        base = row.to_dict()
        base["protparam_sequence"] = standard_seq
        base["protparam_nonstandard_aa_removed"] = int(nonstandard_count)
        base["protparam_nonletter_removed"] = int(removed_count)

        if len(standard_seq) < min_standard_length:
            fail = base.copy()
            fail["protparam_status"] = "failed"
            fail["protparam_error"] = f"standard_sequence_shorter_than_{min_standard_length}"
            failed_rows.append(fail)
            continue

        try:
            features = compute_features_for_sequence(standard_seq)
            out = base.copy()
            out.update(features)
            out["protparam_status"] = "computed"
            out["protparam_error"] = ""
            feature_rows.append(out)

        except Exception as exc:
            fail = base.copy()
            fail["protparam_status"] = "failed"
            fail["protparam_error"] = str(exc)
            failed_rows.append(fail)
            logger.write(f"{accession}: failed ProtParam calculation: {exc}")

    features = pd.DataFrame(feature_rows)

    failed_columns = list(df.columns) + [
        "protparam_sequence",
        "protparam_nonstandard_aa_removed",
        "protparam_nonletter_removed",
        "protparam_status",
        "protparam_error",
    ]

    failed = pd.DataFrame(failed_rows, columns=failed_columns)

    return features, failed


def make_group_summary(features: pd.DataFrame) -> pd.DataFrame:
    if features.empty:
        return pd.DataFrame()

    numeric_cols = [
        "sequence_length_protparam",
        "molecular_weight",
        "aromaticity",
        "instability_index",
        "isoelectric_point",
        "gravy",
        "helix_fraction",
        "turn_fraction",
        "sheet_fraction",
    ]

    rows = []

    group_cols = ["organism_type", "gh_family"]

    for keys, group in features.groupby(group_cols):
        organism_type, gh_family = keys
        row = {
            "organism_type": organism_type,
            "gh_family": gh_family,
            "n": int(len(group)),
        }

        for col in numeric_cols:
            if col in group.columns:
                values = pd.to_numeric(group[col], errors="coerce")
                row[f"{col}_mean"] = float(values.mean())
                row[f"{col}_median"] = float(values.median())
                row[f"{col}_sd"] = float(values.std(ddof=1)) if len(values.dropna()) > 1 else 0.0
                row[f"{col}_min"] = float(values.min())
                row[f"{col}_max"] = float(values.max())

        rows.append(row)

    return pd.DataFrame(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute ProtParam sequence features.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--inventory", required=True, help="Input sequence inventory CSV.")
    parser.add_argument("--output", required=True, help="Output per-protein features CSV.")
    parser.add_argument("--summary", required=True, help="Output group summary CSV.")
    parser.add_argument("--failed", required=True, help="Output failed rows CSV.")
    parser.add_argument("--min-standard-length", type=int, default=30, help="Minimum standard amino-acid length.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    inventory_path = resolve_path(args.inventory, project_root)
    output_path = resolve_path(args.output, project_root)
    summary_path = resolve_path(args.summary, project_root)
    failed_path = resolve_path(args.failed, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "03_sequence_04_protparam.log"),
        project_root,
    )

    if inventory_path is None or output_path is None or summary_path is None or failed_path is None:
        raise ValueError("--inventory, --output, --summary and --failed are required")

    logger = Logger(log_path)

    logger.write("Starting ProtParam feature calculation")
    logger.write(f"Inventory: {inventory_path}")
    logger.write(f"Output: {output_path}")

    if not inventory_path.exists():
        raise FileNotFoundError(f"Input inventory not found: {inventory_path}")

    df = pd.read_csv(inventory_path, dtype=str, keep_default_na=False).fillna("")

    features, failed = compute_features(
        df=df,
        min_standard_length=args.min_standard_length,
        logger=logger,
    )

    summary = make_group_summary(features)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    failed_path.parent.mkdir(parents=True, exist_ok=True)

    features.to_csv(output_path, index=False)
    summary.to_csv(summary_path, index=False)
    failed.to_csv(failed_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "inventory": str(inventory_path),
        "output": str(output_path),
        "summary": str(summary_path),
        "failed": str(failed_path),
        "input_rows": int(len(df)),
        "computed_rows": int(len(features)),
        "failed_rows": int(len(failed)),
        "min_standard_length": int(args.min_standard_length),
        "organism_type_counts": dict(Counter(features.get("organism_type", pd.Series(dtype=str)).astype(str))),
        "gh_family_counts": dict(Counter(features.get("gh_family", pd.Series(dtype=str)).astype(str))),
        "calculation_backend": "Bio.SeqUtils.ProtParam.ProteinAnalysis",
        "method_description": "ProtParam-style descriptors calculated locally using Biopython ProteinAnalysis; sequences were not submitted to the ExPASy ProtParam web server.",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Input rows: {len(df)}")
    logger.write(f"Computed rows: {len(features)}")
    logger.write(f"Failed rows: {len(failed)}")
    logger.write(f"Output: {output_path}")
    logger.write(f"Summary: {summary_path}")
    logger.write("Finished ProtParam feature calculation")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
