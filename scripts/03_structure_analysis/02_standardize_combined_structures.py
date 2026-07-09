#!/usr/bin/env python3
"""
Validate and collect standardized structures for a frozen structural dataset.

This script expects a structure manifest containing usable structure paths.
If structures are already standardized, it copies or symlinks them into a clean
output directory and writes a new standardized manifest.

It does not chemically repair structures. Repair/energy preparation belongs to
FoldX-specific stages.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import shutil
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


def safe_filename(value: str) -> str:
    value = clean_text(value)
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = value.strip("_")
    return value or "structure"


def count_structure_records(path: Path) -> tuple[int, int, int]:
    atom = 0
    hetatm = 0
    model = 0

    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("ATOM"):
                    atom += 1
                elif line.startswith("HETATM"):
                    hetatm += 1
                elif line.startswith("MODEL"):
                    model += 1
    except Exception:
        return 0, 0, 0

    return atom, hetatm, model


def infer_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        return "pdb"
    if suffix == ".cif":
        return "cif"
    if suffix == ".ent":
        return "ent"
    return suffix.replace(".", "") or "unknown"


def build_output_name(row: pd.Series, src: Path) -> str:
    accession = safe_filename(row.get("uniprot_accession", "unknown"))
    source = safe_filename(row.get("structure_source_normalized", row.get("structure_source", "unknown")))
    structure_id = safe_filename(row.get("structure_id", src.stem))
    suffix = src.suffix.lower() or ".pdb"

    return f"{accession}_{source}_{structure_id}{suffix}"


def copy_or_link(src: Path, dst: Path, mode: str, overwrite: bool) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() or dst.is_symlink():
        if overwrite:
            dst.unlink()
        else:
            return "skipped_existing"

    if mode == "copy":
        shutil.copy2(src, dst)
        return "copied"

    if mode == "symlink":
        os.symlink(src, dst)
        return "symlinked"

    if mode == "none":
        return "not_copied"

    raise ValueError(f"Unsupported copy mode: {mode}")


def process_manifest(
    manifest: pd.DataFrame,
    output_dir: Path,
    copy_mode: str,
    overwrite: bool,
    logger: Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "usable_structure_path" not in manifest.columns:
        raise ValueError("Manifest must contain usable_structure_path column")

    rows = []
    problem_rows = []

    for _, row in manifest.iterrows():
        src_text = clean_text(row.get("usable_structure_path", ""))
        out = row.to_dict()

        status = ""
        error = ""
        operation = ""
        stage3_path = ""

        if not src_text:
            status = "failed_missing_path"
            error = "usable_structure_path_empty"
        else:
            src = Path(src_text)

            if not src.exists():
                status = "failed_missing_file"
                error = "source_file_not_found"
            else:
                atom_count, hetatm_count, model_count = count_structure_records(src)
                fmt = infer_format(src)

                output_name = build_output_name(row, src)
                dst = output_dir / output_name

                try:
                    operation = copy_or_link(src, dst, mode=copy_mode, overwrite=overwrite)
                    stage3_path = str(dst) if copy_mode in {"copy", "symlink"} else str(src)

                    if fmt == "pdb" and atom_count == 0:
                        status = "warning_no_atom_records"
                    else:
                        status = "standardized_available"

                    out["structure_file_format"] = fmt
                    out["n_atom_records"] = atom_count
                    out["n_hetatm_records"] = hetatm_count
                    out["n_model_records"] = model_count

                except Exception as exc:
                    status = "failed_copy_or_link"
                    error = str(exc)

        out["stage3_standardized_structure_path"] = stage3_path
        out["stage3_standardization_status"] = status
        out["stage3_copy_mode"] = copy_mode
        out["stage3_copy_operation"] = operation
        out["stage3_standardization_error"] = error

        rows.append(out)

        if status.startswith("failed") or status.startswith("warning"):
            problem_rows.append(out)

        logger.write(f"{out.get('uniprot_accession', '')} {out.get('structure_id', '')}: {status}")

    standardized_df = pd.DataFrame(rows)
    problem_df = pd.DataFrame(problem_rows, columns=standardized_df.columns)

    return standardized_df, problem_df


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate and collect standardized structures.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--manifest", required=True, help="Input combined structure manifest CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for collected standardized structures.")
    parser.add_argument("--output-manifest", required=True, help="Output standardized manifest CSV.")
    parser.add_argument("--problem-output", required=True, help="Output failed/warning structure rows CSV.")
    parser.add_argument("--copy-mode", choices=["copy", "symlink", "none"], default="copy", help="How to collect files.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing collected files.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    manifest_path = resolve_path(args.manifest, project_root)
    output_dir = resolve_path(args.output_dir, project_root)
    output_manifest_path = resolve_path(args.output_manifest, project_root)
    problem_output_path = resolve_path(args.problem_output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "04_structure_02_standardize.log"),
        project_root,
    )

    logger = Logger(log_path)

    logger.write("Starting structure standardization/collection")
    logger.write(f"Input manifest: {manifest_path}")
    logger.write(f"Output directory: {output_dir}")
    logger.write(f"Output manifest: {output_manifest_path}")
    logger.write(f"Copy mode: {args.copy_mode}")

    if manifest_path is None or not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False).fillna("")

    standardized, problems = process_manifest(
        manifest=manifest,
        output_dir=output_dir,
        copy_mode=args.copy_mode,
        overwrite=args.overwrite,
        logger=logger,
    )

    output_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    problem_output_path.parent.mkdir(parents=True, exist_ok=True)

    standardized.to_csv(output_manifest_path, index=False)
    problems.to_csv(problem_output_path, index=False)

    status_counts = dict(Counter(standardized["stage3_standardization_status"].astype(str)))

    metadata = {
        "script": Path(__file__).name,
        "manifest": str(manifest_path),
        "output_dir": str(output_dir),
        "output_manifest": str(output_manifest_path),
        "problem_output": str(problem_output_path),
        "copy_mode": args.copy_mode,
        "input_rows": int(len(manifest)),
        "output_rows": int(len(standardized)),
        "problem_rows": int(len(problems)),
        "status_counts": status_counts,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_manifest_path.with_suffix(output_manifest_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Input rows: {len(manifest)}")
    logger.write(f"Output rows: {len(standardized)}")
    logger.write(f"Problem rows: {len(problems)}")
    logger.write(f"Status counts: {status_counts}")
    logger.write("Finished structure standardization/collection")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
