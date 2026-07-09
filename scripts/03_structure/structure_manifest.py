#!/usr/bin/env python3
# Purpose: Build the structure manifest.

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import sys
import time

import pandas as pd


STRUCTURE_SOURCE_PRIORITY = {
    "pdb": 1,
    "experimental_pdb": 1,
    "mmcif": 1,
    "cif": 1,
    "swissmodel": 2,
    "alphafold": 3,
    "modeller": 4,
    "model": 5,
    "unknown": 99,
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


def normalise_accession(value) -> str:
    return clean_text(value).split("-")[0]


def normalise_source(value) -> str:
    text = clean_text(value).lower()
    if not text:
        return "unknown"
    if "pdb" in text:
        return "pdb"
    if "cif" in text or "mmcif" in text:
        return "pdb"
    if "swiss" in text:
        return "swissmodel"
    if "alpha" in text:
        return "alphafold"
    if "modeller" in text:
        return "modeller"
    return text


def path_exists(value: str) -> bool:
    text = clean_text(value)
    if not text:
        return False
    return Path(text).exists()


def choose_usable_path(row: pd.Series) -> str:
    candidates = [
        "standardized_structure_path",
        "foldx_input_pdb",
        "model_path",
        "structure_path",
        "pdb_path",
        "file_path",
        "local_path",
        "final_structure_path",
    ]

    for col in candidates:
        if col in row.index:
            val = clean_text(row.get(col, ""))
            if val and Path(val).exists():
                return val

    # If no existing file is found, still keep the first non-empty candidate for diagnostics
    for col in candidates:
        if col in row.index:
            val = clean_text(row.get(col, ""))
            if val:
                return val

    return ""


def source_priority(source: str) -> int:
    return STRUCTURE_SOURCE_PRIORITY.get(normalise_source(source), 99)


def read_csv_required(path: Path, description: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=False).fillna("")


def prepare_structure_manifest(
    master: pd.DataFrame,
    source_manifest: pd.DataFrame,
    logger: Logger,
) -> pd.DataFrame:
    if "uniprot_accession" not in master.columns:
        raise ValueError("Master table must contain uniprot_accession")

    if "uniprot_accession" not in source_manifest.columns:
        raise ValueError("Source manifest must contain uniprot_accession")

    master_accessions = set(master["uniprot_accession"].map(normalise_accession))

    source = source_manifest.copy()
    source["uniprot_accession"] = source["uniprot_accession"].map(normalise_accession)

    before_filter = len(source)
    source = source[source["uniprot_accession"].isin(master_accessions)].copy()
    after_filter = len(source)

    logger.write(f"Source manifest rows before frozen-master filter: {before_filter}")
    logger.write(f"Source manifest rows after frozen-master filter: {after_filter}")

    if "structure_source" not in source.columns:
        source["structure_source"] = "unknown"

    if "structure_id" not in source.columns:
        source["structure_id"] = source["uniprot_accession"] + "_structure"

    source["structure_source_normalized"] = source["structure_source"].map(normalise_source)
    source["structure_source_priority"] = source["structure_source_normalized"].map(source_priority)

    source["usable_structure_path"] = source.apply(choose_usable_path, axis=1)
    source["usable_structure_exists"] = source["usable_structure_path"].map(path_exists)
    source["usable_structure_status"] = source["usable_structure_exists"].map(lambda x: "available" if x else "missing_file")

    if "standardization_status" not in source.columns:
        source["standardization_status"] = ""

    if "structure_exists" not in source.columns:
        source["structure_exists"] = ""

    source["structure_manifest_row_id"] = [
        f"structure_{i:06d}" for i in range(1, len(source) + 1)
    ]

    # Put available structures first, then source priority, then accession
    source = source.sort_values(
        ["usable_structure_exists", "structure_source_priority", "uniprot_accession", "structure_id"],
        ascending=[False, True, True, True],
    ).copy()

    return source


def build_protein_availability(master: pd.DataFrame, structure_manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    master = master.copy()
    master["uniprot_accession"] = master["uniprot_accession"].map(normalise_accession)

    available = structure_manifest[structure_manifest["usable_structure_exists"] == True].copy()

    rows = []

    for _, mrow in master.iterrows():
        acc = normalise_accession(mrow.get("uniprot_accession", ""))
        group = available[available["uniprot_accession"] == acc].copy()

        row = {
            "uniprot_accession": acc,
            "organism": clean_text(mrow.get("organism", "")),
            "organism_type": clean_text(mrow.get("organism_type", "")),
            "gh_family": clean_text(mrow.get("gh_family", "")),
            "has_usable_structure": bool(len(group) > 0),
            "n_structure_entries": int(len(group)),
            "structure_sources": "",
            "best_structure_source": "",
            "best_structure_id": "",
            "best_structure_path": "",
            "best_standardization_status": "",
        }

        if len(group):
            group = group.sort_values(
                ["structure_source_priority", "structure_id"],
                ascending=[True, True],
            ).copy()

            best = group.iloc[0]

            row["structure_sources"] = ";".join(sorted(set(group["structure_source_normalized"].astype(str))))
            row["best_structure_source"] = clean_text(best.get("structure_source_normalized", ""))
            row["best_structure_id"] = clean_text(best.get("structure_id", ""))
            row["best_structure_path"] = clean_text(best.get("usable_structure_path", ""))
            row["best_standardization_status"] = clean_text(best.get("standardization_status", ""))

        rows.append(row)

    availability = pd.DataFrame(rows)
    unresolved = availability[availability["has_usable_structure"] == False].copy()

    return availability, unresolved


def summarise(structure_manifest: pd.DataFrame, availability: pd.DataFrame, unresolved: pd.DataFrame) -> dict:
    available_structures = structure_manifest[structure_manifest["usable_structure_exists"] == True].copy()

    return {
        "structure_manifest_rows": int(len(structure_manifest)),
        "available_structure_rows": int(len(available_structures)),
        "missing_file_structure_rows": int((structure_manifest["usable_structure_exists"] == False).sum()),
        "unique_accessions_with_structure_rows": int(available_structures["uniprot_accession"].nunique()) if len(available_structures) else 0,
        "protein_rows": int(len(availability)),
        "proteins_with_usable_structure": int(availability["has_usable_structure"].sum()),
        "proteins_without_usable_structure": int(len(unresolved)),
        "structure_source_counts": dict(Counter(available_structures.get("structure_source_normalized", pd.Series(dtype=str)).astype(str))),
        "organism_type_counts": dict(Counter(availability.get("organism_type", pd.Series(dtype=str)).astype(str))),
        "gh_family_counts": dict(Counter(availability.get("gh_family", pd.Series(dtype=str)).astype(str))),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build combined structure manifest for a curated master dataset.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--master", required=True, help="Curated master CSV.")
    parser.add_argument("--source-manifest", required=True, help="Existing combined/standardized structure manifest CSV.")
    parser.add_argument("--output", required=True, help="Output structure-level manifest CSV.")
    parser.add_argument("--protein-availability-output", required=True, help="Output protein-level structure availability CSV.")
    parser.add_argument("--unresolved-output", required=True, help="Output proteins without usable structure CSV.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    master_path = resolve_path(args.master, project_root)
    source_manifest_path = resolve_path(args.source_manifest, project_root)
    output_path = resolve_path(args.output, project_root)
    protein_output_path = resolve_path(args.protein_availability_output, project_root)
    unresolved_output_path = resolve_path(args.unresolved_output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "04_structure_01_build_manifest.log"),
        project_root,
    )

    logger = Logger(log_path)

    logger.write("Starting combined structure manifest build")
    logger.write(f"Master: {master_path}")
    logger.write(f"Source manifest: {source_manifest_path}")
    logger.write(f"Output: {output_path}")

    if master_path is None or source_manifest_path is None:
        raise ValueError("--master and --source-manifest are required")

    master = read_csv_required(master_path, "Master file")
    source_manifest = read_csv_required(source_manifest_path, "Source manifest")

    structure_manifest = prepare_structure_manifest(
        master=master,
        source_manifest=source_manifest,
        logger=logger,
    )

    availability, unresolved = build_protein_availability(
        master=master,
        structure_manifest=structure_manifest,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    protein_output_path.parent.mkdir(parents=True, exist_ok=True)
    unresolved_output_path.parent.mkdir(parents=True, exist_ok=True)

    structure_manifest.to_csv(output_path, index=False)
    availability.to_csv(protein_output_path, index=False)
    unresolved.to_csv(unresolved_output_path, index=False)

    metadata = {
        "script": Path(__file__).name,
        "master": str(master_path),
        "source_manifest": str(source_manifest_path),
        "output": str(output_path),
        "protein_availability_output": str(protein_output_path),
        "unresolved_output": str(unresolved_output_path),
        "summary": summarise(structure_manifest, availability, unresolved),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    summary = metadata["summary"]

    logger.write(f"Structure manifest rows: {summary['structure_manifest_rows']}")
    logger.write(f"Available structure rows: {summary['available_structure_rows']}")
    logger.write(f"Proteins with usable structure: {summary['proteins_with_usable_structure']}")
    logger.write(f"Proteins without usable structure: {summary['proteins_without_usable_structure']}")
    logger.write("Finished combined structure manifest build")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
