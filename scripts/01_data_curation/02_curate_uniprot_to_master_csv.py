#!/usr/bin/env python3
"""
Curate UniProt TSV records into a standard enzyme master table.

This script is generic and can be used for xylanase or other enzyme families.
It reads a UniProt TSV file, normalises important columns, classifies organism
type from taxonomy, extracts PDB/RefSeq cross-references, infers family labels
from configurable patterns, removes duplicate accessions, and writes a curation summary
metadata file.

Example:
    python scripts/01_curation/02_curate_uniprot_to_master_csv.py \
        --config config/xylanase_config.yaml \
        --input data/raw/uniprot_xylanase.tsv \
        --output data/curated/uniprot_xylanase_curated.csv
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys
import time
from collections import Counter

import pandas as pd


COLUMN_ALIASES = {
    "accession": ["Entry", "accession", "primary_accession", "uniprot_accession"],
    "entry_name": ["Entry Name", "entry_name", "id"],
    "protein_name": ["Protein names", "protein_name", "protein_names"],
    "gene_names": ["Gene Names", "gene_names", "genes"],
    "organism": ["Organism", "organism", "organism_name"],
    "organism_id": ["Organism (ID)", "organism_id", "taxon_id"],
    "lineage": ["Taxonomic lineage", "lineage", "taxonomic_lineage"],
    "length": ["Length", "length"],
    "sequence": ["Sequence", "sequence"],
    "pdb": ["PDB", "xref_pdb", "pdb"],
    "refseq": ["RefSeq", "xref_refseq", "refseq"],
    "cazy": ["CAZy", "xref_cazy", "cazy"],
}


DEFAULT_FAMILY_PATTERNS = {
    "GH10": [
        r"\bGH[-_ ]?10\b",
        r"glycoside hydrolase family 10",
        r"xylanase 10",
        r"xyn10",
        r"xyl10",
    ],
    "GH11": [
        r"\bGH[-_ ]?11\b",
        r"glycoside hydrolase family 11",
        r"xylanase 11",
        r"xyn11",
        r"xyl11",
    ],
}


def load_yaml_config(path: str | None) -> dict:
    if not path:
        return {}

    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        import yaml
    except ImportError as exc:
        raise ImportError(
            "PyYAML is required to use --config. Install with: pip install pyyaml"
        ) from exc

    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    if not isinstance(data, dict):
        raise ValueError(f"Config file must contain a YAML mapping: {config_path}")

    return data


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


def find_column(df: pd.DataFrame, logical_name: str) -> str | None:
    """Find a dataframe column using known aliases."""
    aliases = COLUMN_ALIASES.get(logical_name, [logical_name])

    # exact match first
    for alias in aliases:
        if alias in df.columns:
            return alias

    # case-insensitive fallback
    lower_map = {str(col).lower(): col for col in df.columns}
    for alias in aliases:
        if alias.lower() in lower_map:
            return lower_map[alias.lower()]

    return None


def get_series(df: pd.DataFrame, logical_name: str, default: str = "") -> pd.Series:
    col = find_column(df, logical_name)
    if col is None:
        return pd.Series([default] * len(df), index=df.index, dtype="object")
    return df[col].fillna("").astype(str)


def clean_text(value) -> str:
    if value is None:
        return ""
    text = str(value)
    text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def split_xrefs(value: str) -> list[str]:
    """Split UniProt xref fields such as '1ABC; 2XYZ;'."""
    value = clean_text(value)
    if not value:
        return []

    parts = re.split(r"[;,|]\s*", value)
    cleaned = []
    for part in parts:
        part = part.strip()
        if part:
            cleaned.append(part)
    return cleaned


def first_or_empty(values: list[str]) -> str:
    return values[0] if values else ""


def classify_organism_type(lineage: str, organism: str = "") -> str:
    """Classify organism type from UniProt lineage."""
    text = f"{lineage} {organism}".lower()

    if "bacteria" in text or "bacterium" in text:
        return "bacterial"

    if "fungi" in text or "fungus" in text or "ascomycota" in text or "basidiomycota" in text:
        return "fungal"

    if "archaea" in text or "archaeon" in text:
        return "archaeal"

    if "eukaryota" in text:
        return "eukaryotic_other"

    return "unknown"


def normalise_family_patterns(raw_patterns) -> dict[str, list[str]]:
    """Read family regex patterns from config or use defaults."""
    if not raw_patterns:
        return DEFAULT_FAMILY_PATTERNS

    if not isinstance(raw_patterns, dict):
        raise ValueError("curation.family_patterns must be a mapping of family -> list of patterns")

    result = {}
    for family, patterns in raw_patterns.items():
        if isinstance(patterns, str):
            result[str(family)] = [patterns]
        elif isinstance(patterns, list):
            result[str(family)] = [str(p) for p in patterns]
        else:
            raise ValueError(f"Patterns for family {family} must be string or list")
    return result


def infer_family(row_text: str, family_patterns: dict[str, list[str]]) -> str:
    """Infer family label from combined annotation text."""
    text = row_text.lower()

    for family, patterns in family_patterns.items():
        for pattern in patterns:
            if re.search(pattern, text, flags=re.IGNORECASE):
                return family

    return "unknown"


def safe_int(value) -> int | None:
    try:
        if value is None:
            return None
        text = str(value).strip()
        if not text:
            return None
        return int(float(text))
    except Exception:
        return None


def build_curated_table(
    raw: pd.DataFrame,
    family_patterns: dict[str, list[str]],
    keep_unknown_family: bool,
    accepted_organism_groups: list[str] | None,
) -> pd.DataFrame:
    """Build standard curated table from UniProt TSV dataframe."""

    accession = get_series(raw, "accession").map(clean_text)
    entry_name = get_series(raw, "entry_name").map(clean_text)
    protein_name = get_series(raw, "protein_name").map(clean_text)
    gene_names = get_series(raw, "gene_names").map(clean_text)
    organism = get_series(raw, "organism").map(clean_text)
    organism_id = get_series(raw, "organism_id").map(clean_text)
    lineage = get_series(raw, "lineage").map(clean_text)
    sequence = get_series(raw, "sequence").map(lambda x: clean_text(x).replace(" ", ""))
    length_raw = get_series(raw, "length").map(clean_text)
    pdb_raw = get_series(raw, "pdb").map(clean_text)
    refseq_raw = get_series(raw, "refseq").map(clean_text)
    cazy_raw = get_series(raw, "cazy").map(clean_text)

    curated = pd.DataFrame({
        "source": "uniprot",
        "primary_id": accession,
        "uniprot_accession": accession,
        "entry_name": entry_name,
        "protein_name": protein_name,
        "gene_names": gene_names,
        "organism": organism,
        "organism_id": organism_id,
        "lineage": lineage,
        "length": length_raw.map(safe_int),
        "sequence": sequence,
        "xref_pdb": pdb_raw,
        "xref_refseq": refseq_raw,
        "xref_cazy": cazy_raw,
    })

    curated["pdb_ids"] = curated["xref_pdb"].map(lambda x: ";".join(split_xrefs(x)))
    curated["pdb_id_first"] = curated["xref_pdb"].map(lambda x: first_or_empty(split_xrefs(x)))
    curated["refseq_ids"] = curated["xref_refseq"].map(lambda x: ";".join(split_xrefs(x)))
    curated["refseq_acc"] = curated["xref_refseq"].map(lambda x: first_or_empty(split_xrefs(x)))

    curated["organism_type"] = [
        classify_organism_type(lin, org)
        for lin, org in zip(curated["lineage"], curated["organism"])
    ]

    family_text = (
        curated["protein_name"].fillna("") + " " +
        curated["entry_name"].fillna("") + " " +
        curated["gene_names"].fillna("") + " " +
        curated["xref_cazy"].fillna("")
    )

    curated["gh_family"] = family_text.map(lambda x: infer_family(x, family_patterns))

    # Basic quality flags
    curated["has_sequence"] = curated["sequence"].str.len() > 0
    curated["has_pdb"] = curated["pdb_ids"].str.len() > 0
    curated["has_refseq"] = curated["refseq_ids"].str.len() > 0
    curated["sequence_length_observed"] = curated["sequence"].str.len()

    # Remove empty accessions and duplicate accessions
    curated = curated[curated["uniprot_accession"].astype(str).str.len() > 0].copy()
    curated = curated.drop_duplicates(subset=["uniprot_accession"], keep="first").copy()

    if accepted_organism_groups:
        accepted = {str(x).lower() for x in accepted_organism_groups}
        curated = curated[curated["organism_type"].str.lower().isin(accepted)].copy()

    if not keep_unknown_family:
        curated = curated[curated["gh_family"] != "unknown"].copy()

    # Put core columns first
    core_cols = [
        "source",
        "primary_id",
        "uniprot_accession",
        "refseq_acc",
        "organism_type",
        "organism",
        "gh_family",
        "sequence",
        "entry_name",
        "protein_name",
        "gene_names",
        "length",
        "sequence_length_observed",
        "organism_id",
        "lineage",
        "xref_pdb",
        "pdb_ids",
        "pdb_id_first",
        "xref_refseq",
        "refseq_ids",
        "xref_cazy",
        "has_sequence",
        "has_pdb",
        "has_refseq",
    ]

    other_cols = [c for c in curated.columns if c not in core_cols]
    curated = curated[core_cols + other_cols]

    return curated


def write_metadata(
    output_path: Path,
    metadata: dict,
) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def summarize_counts(df: pd.DataFrame) -> dict:
    return {
        "rows": int(len(df)),
        "organism_type_counts": dict(Counter(df.get("organism_type", pd.Series(dtype=str)).fillna("").astype(str))),
        "family_counts": dict(Counter(df.get("gh_family", pd.Series(dtype=str)).fillna("").astype(str))),
        "has_sequence_counts": dict(Counter(df.get("has_sequence", pd.Series(dtype=bool)).astype(str))),
        "has_pdb_counts": dict(Counter(df.get("has_pdb", pd.Series(dtype=bool)).astype(str))),
        "has_refseq_counts": dict(Counter(df.get("has_refseq", pd.Series(dtype=bool)).astype(str))),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Curate UniProt TSV into a standard enzyme master CSV."
    )

    parser.add_argument("--config", default=None, help="Optional YAML config file.")
    parser.add_argument("--input", required=False, help="Input UniProt TSV file.")
    parser.add_argument("--output", required=False, help="Output curated CSV file.")
    parser.add_argument("--log", required=False, help="Log file path.")
    parser.add_argument(
        "--drop-unknown-family",
        action="store_true",
        help="Drop rows where family could not be inferred.",
    )
    parser.add_argument(
        "--keep-all-organisms",
        action="store_true",
        help="Do not filter to accepted organism groups from config.",
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(
        args.input
        or deep_get(config, ["outputs", "raw_uniprot_tsv"], None),
        project_root,
    )

    output_path = resolve_path(
        args.output
        or deep_get(config, ["outputs", "curated_uniprot_csv"], None)
        or "data/curated/uniprot_curated.csv",
        project_root,
    )

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "02_curate_uniprot_to_master_csv.log"),
        project_root,
    )

    if input_path is None:
        raise ValueError("Input path is required. Use --input or config outputs.raw_uniprot_tsv.")

    if output_path is None:
        raise ValueError("Output path is required. Use --output or config outputs.curated_uniprot_csv.")

    logger = Logger(log_path)

    logger.write("Starting UniProt curation")
    logger.write(f"Input: {input_path}")
    logger.write(f"Output: {output_path}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input UniProt TSV not found: {input_path}")

    raw = pd.read_csv(input_path, sep="\t", dtype=str, keep_default_na=False)

    logger.write(f"Raw rows: {len(raw)}")
    logger.write(f"Raw columns: {list(raw.columns)}")

    family_patterns = normalise_family_patterns(
        deep_get(config, ["curation", "family_patterns"], None)
    )

    keep_unknown_family = bool(
        deep_get(config, ["curation", "keep_unknown_family"], True)
    )

    if args.drop_unknown_family:
        keep_unknown_family = False

    accepted_organism_groups = deep_get(
        config,
        ["curation", "accepted_organism_groups"],
        None,
    )

    if args.keep_all_organisms:
        accepted_organism_groups = None

    curated = build_curated_table(
        raw=raw,
        family_patterns=family_patterns,
        keep_unknown_family=keep_unknown_family,
        accepted_organism_groups=accepted_organism_groups,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    curated.to_csv(output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "input": str(input_path),
        "output": str(output_path),
        "raw_rows": int(len(raw)),
        "curated_rows": int(len(curated)),
        "raw_columns": list(raw.columns),
        "family_patterns": family_patterns,
        "keep_unknown_family": keep_unknown_family,
        "accepted_organism_groups": accepted_organism_groups,
        "summary": summarize_counts(curated),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Curated rows: {len(curated)}")
    logger.write(f"Organism counts: {metadata['summary']['organism_type_counts']}")
    logger.write(f"Family counts: {metadata['summary']['family_counts']}")
    logger.write("Finished UniProt curation successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
