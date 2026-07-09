#!/usr/bin/env python3
"""
Prepare TM-align pairwise job manifest.

Default thesis logic:
    - query structures: modelled structures, e.g. MODELLER
    - reference structures: experimental PDB-derived structures
    - compare only within the same GH family

The output is a job manifest that can be consumed by a TM-align runner.
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
    text = str(value).strip()
    if text.lower() in {"nan", "none", "na", "n/a", "null"}:
        return ""
    return text


def normalise_source(value) -> str:
    text = clean_text(value).lower()
    if "pdb" in text:
        return "pdb"
    if "modeller" in text:
        return "modeller"
    if "swiss" in text:
        return "swissmodel"
    if "alpha" in text:
        return "alphafold"
    return text or "unknown"


def pick_path(row: pd.Series) -> str:
    candidates = [
        "stage3_standardized_structure_path",
        "structure_path_for_features",
        "usable_structure_path",
        "standardized_structure_path",
        "foldx_input_pdb",
        "model_path",
        "structure_path",
    ]

    for col in candidates:
        if col in row.index:
            value = clean_text(row.get(col, ""))
            if value and Path(value).exists():
                return value

    for col in candidates:
        if col in row.index:
            value = clean_text(row.get(col, ""))
            if value:
                return value

    return ""


def select_rows(
    df: pd.DataFrame,
    query_source: str,
    reference_source: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()

    if "structure_source_normalized" not in df.columns:
        if "structure_source" in df.columns:
            df["structure_source_normalized"] = df["structure_source"].map(normalise_source)
        else:
            df["structure_source_normalized"] = "unknown"

    df["structure_source_normalized"] = df["structure_source_normalized"].map(normalise_source)
    df["tmalign_structure_path"] = df.apply(pick_path, axis=1)
    df["tmalign_structure_exists"] = df["tmalign_structure_path"].map(lambda p: Path(p).exists() if p else False)

    if "feature_status" in df.columns:
        df = df[df["feature_status"] == "computed"].copy()

    df = df[df["tmalign_structure_exists"] == True].copy()

    if query_source == "non_pdb":
        query = df[df["structure_source_normalized"] != "pdb"].copy()
    elif query_source == "all":
        query = df.copy()
    else:
        query = df[df["structure_source_normalized"] == query_source].copy()

    if reference_source == "all":
        reference = df.copy()
    else:
        reference = df[df["structure_source_normalized"] == reference_source].copy()

    return query, reference


def build_jobs(
    query: pd.DataFrame,
    reference: pd.DataFrame,
    match_gh_family: bool,
    match_organism_type: bool,
    include_same_accession: bool,
    max_references_per_query: int | None,
) -> pd.DataFrame:
    jobs = []

    query = query.reset_index(drop=True)
    reference = reference.reset_index(drop=True)

    for _, q in query.iterrows():
        refs = reference.copy()

        if match_gh_family:
            refs = refs[refs["gh_family"].astype(str) == str(q.get("gh_family", ""))].copy()

        if match_organism_type:
            refs = refs[refs["organism_type"].astype(str) == str(q.get("organism_type", ""))].copy()

        if not include_same_accession:
            refs = refs[refs["uniprot_accession"].astype(str) != str(q.get("uniprot_accession", ""))].copy()

        refs = refs.sort_values(["uniprot_accession", "structure_id"]).copy()

        if max_references_per_query is not None and max_references_per_query > 0:
            refs = refs.head(max_references_per_query).copy()

        for _, r in refs.iterrows():
            jobs.append({
                "job_id": f"tmalign_{len(jobs) + 1:08d}",
                "query_uniprot_accession": clean_text(q.get("uniprot_accession", "")),
                "query_structure_id": clean_text(q.get("structure_id", "")),
                "query_structure_source": clean_text(q.get("structure_source_normalized", "")),
                "query_organism": clean_text(q.get("organism", "")),
                "query_organism_type": clean_text(q.get("organism_type", "")),
                "query_gh_family": clean_text(q.get("gh_family", "")),
                "query_structure_path": clean_text(q.get("tmalign_structure_path", "")),
                "reference_uniprot_accession": clean_text(r.get("uniprot_accession", "")),
                "reference_structure_id": clean_text(r.get("structure_id", "")),
                "reference_structure_source": clean_text(r.get("structure_source_normalized", "")),
                "reference_organism": clean_text(r.get("organism", "")),
                "reference_organism_type": clean_text(r.get("organism_type", "")),
                "reference_gh_family": clean_text(r.get("gh_family", "")),
                "reference_structure_path": clean_text(r.get("tmalign_structure_path", "")),
                "same_gh_family": str(q.get("gh_family", "")) == str(r.get("gh_family", "")),
                "same_organism_type": str(q.get("organism_type", "")) == str(r.get("organism_type", "")),
                "same_uniprot_accession": str(q.get("uniprot_accession", "")) == str(r.get("uniprot_accession", "")),
            })

    return pd.DataFrame(jobs)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare TM-align pairwise job manifest.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--input", required=True, help="Input structure feature/manifest CSV.")
    parser.add_argument("--output", required=True, help="Output TM-align job manifest CSV.")
    parser.add_argument("--query-source", default="modeller", help="Query source: modeller, non_pdb, all, etc.")
    parser.add_argument("--reference-source", default="pdb", help="Reference source: pdb, all, etc.")
    parser.add_argument("--match-gh-family", action="store_true", help="Only compare within same GH family.")
    parser.add_argument("--match-organism-type", action="store_true", help="Only compare within same organism type.")
    parser.add_argument("--include-same-accession", action="store_true", help="Allow query/reference from same UniProt accession.")
    parser.add_argument("--max-references-per-query", type=int, default=0, help="Optional cap per query; 0 means no cap.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(args.input, project_root)
    output_path = resolve_path(args.output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "04_structure_04_prepare_tmalign_jobs.log"),
        project_root,
    )

    logger = Logger(log_path)

    logger.write("Starting TM-align job preparation")
    logger.write(f"Input: {input_path}")
    logger.write(f"Output: {output_path}")
    logger.write(f"Query source: {args.query_source}")
    logger.write(f"Reference source: {args.reference_source}")
    logger.write(f"Match GH family: {args.match_gh_family}")
    logger.write(f"Match organism type: {args.match_organism_type}")

    if input_path is None or not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(input_path, dtype=str, keep_default_na=False).fillna("")

    query, reference = select_rows(
        df=df,
        query_source=args.query_source,
        reference_source=args.reference_source,
    )

    max_refs = args.max_references_per_query if args.max_references_per_query and args.max_references_per_query > 0 else None

    jobs = build_jobs(
        query=query,
        reference=reference,
        match_gh_family=args.match_gh_family,
        match_organism_type=args.match_organism_type,
        include_same_accession=args.include_same_accession,
        max_references_per_query=max_refs,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    jobs.to_csv(output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "input": str(input_path),
        "output": str(output_path),
        "input_rows": int(len(df)),
        "query_rows": int(len(query)),
        "reference_rows": int(len(reference)),
        "job_rows": int(len(jobs)),
        "query_source": args.query_source,
        "reference_source": args.reference_source,
        "match_gh_family": args.match_gh_family,
        "match_organism_type": args.match_organism_type,
        "include_same_accession": args.include_same_accession,
        "max_references_per_query": max_refs,
        "query_gh_counts": dict(Counter(query.get("gh_family", pd.Series(dtype=str)).astype(str))),
        "reference_gh_counts": dict(Counter(reference.get("gh_family", pd.Series(dtype=str)).astype(str))),
        "job_gh_counts": dict(Counter(jobs.get("query_gh_family", pd.Series(dtype=str)).astype(str))) if len(jobs) else {},
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Input rows: {len(df)}")
    logger.write(f"Query rows: {len(query)}")
    logger.write(f"Reference rows: {len(reference)}")
    logger.write(f"Job rows: {len(jobs)}")
    logger.write("Finished TM-align job preparation")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
