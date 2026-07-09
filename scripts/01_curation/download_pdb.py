#!/usr/bin/env python3
# Purpose: Download PDB structure files.

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import sys
import time
import urllib.request
from urllib.error import HTTPError, URLError

import pandas as pd


RCSB_DOWNLOAD_BASE = "https://files.rcsb.org/download"


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


def normalise_pdb_id(value) -> str:
    value = clean_text(value).upper()
    value = re.sub(r"[^A-Z0-9]", "", value)
    if len(value) == 4:
        return value
    return ""


def build_download_url(pdb_id: str, fmt: str) -> str:
    fmt = fmt.lower()
    if fmt == "pdb":
        return f"{RCSB_DOWNLOAD_BASE}/{pdb_id}.pdb"
    if fmt in {"cif", "mmcif"}:
        return f"{RCSB_DOWNLOAD_BASE}/{pdb_id}.cif"
    raise ValueError(f"Unsupported format: {fmt}")


def download_file(
    url: str,
    output_path: Path,
    user_agent: str,
    timeout: int,
    retries: int,
    logger: Logger,
) -> tuple[str, str, int]:
    """
    Download a file.

    Returns:
        status, error, bytes_written
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    last_error = ""

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": user_agent})

            with urllib.request.urlopen(request, timeout=timeout) as response:
                data = response.read()

            with output_path.open("wb") as handle:
                handle.write(data)

            return "downloaded", "", len(data)

        except HTTPError as exc:
            last_error = f"HTTPError {exc.code}: {exc.reason}"
            logger.write(f"{last_error} for {url}")

        except URLError as exc:
            last_error = f"URLError: {exc}"
            logger.write(f"{last_error} for {url}")

        except TimeoutError as exc:
            last_error = f"TimeoutError: {exc}"
            logger.write(f"{last_error} for {url}")

        if attempt < retries:
            time.sleep(min(5 * attempt, 30))

    return "failed", last_error, 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download PDB/mmCIF files from PDB inventory.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--inventory", required=True, help="Input PDB inventory CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for downloaded structures.")
    parser.add_argument("--manifest", required=True, help="Output download manifest CSV.")
    parser.add_argument("--log", default=None, help="Log file path.")
    parser.add_argument("--format", choices=["pdb", "cif"], default="pdb", help="Download format.")
    parser.add_argument(
        "--fallback-cif",
        action="store_true",
        help="If PDB download fails, automatically try mmCIF for the same PDB ID.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    parser.add_argument("--timeout", type=int, default=120, help="HTTP timeout in seconds.")
    parser.add_argument("--retries", type=int, default=3, help="Number of retries.")
    parser.add_argument("--sleep", type=float, default=0.05, help="Delay between requests.")
    parser.add_argument("--user-agent", default="enzyme-thermostability-pipeline/1.0", help="HTTP User-Agent.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    inventory_path = resolve_path(args.inventory, project_root)
    output_dir = resolve_path(args.output_dir, project_root)
    manifest_path = resolve_path(args.manifest, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "05_download_pdb_structures.log"),
        project_root,
    )

    if inventory_path is None or output_dir is None or manifest_path is None:
        raise ValueError("--inventory, --output-dir and --manifest are required")

    logger = Logger(log_path)

    logger.write("Starting PDB structure download")
    logger.write(f"Inventory: {inventory_path}")
    logger.write(f"Output directory: {output_dir}")
    logger.write(f"Manifest: {manifest_path}")
    logger.write(f"Format: {args.format}")

    if not inventory_path.exists():
        raise FileNotFoundError(f"PDB inventory not found: {inventory_path}")

    inv = pd.read_csv(inventory_path, dtype=str, keep_default_na=False).fillna("")

    if "pdb_id" not in inv.columns:
        raise ValueError("Inventory must contain column: pdb_id")

    inv["pdb_id"] = inv["pdb_id"].map(normalise_pdb_id)
    inv = inv[inv["pdb_id"].str.len() == 4].copy()
    inv = inv.drop_duplicates(subset=["pdb_id"], keep="first").copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    logger.write(f"Unique PDB IDs to process: {len(inv)}")

    rows = []

    extension = "pdb" if args.format == "pdb" else "cif"

    for i, row in inv.iterrows():
        pdb_id = row["pdb_id"]
        out_file = output_dir / f"{pdb_id}.{extension}"
        url = build_download_url(pdb_id, args.format)

        status = ""
        error = ""
        bytes_written = 0

        if out_file.exists() and out_file.stat().st_size > 0 and not args.overwrite:
            status = "skipped_existing"
            bytes_written = out_file.stat().st_size
        else:
            status, error, bytes_written = download_file(
                url=url,
                output_path=out_file,
                user_agent=args.user_agent,
                timeout=args.timeout,
                retries=args.retries,
                logger=logger,
            )

            # Optional fallback: if legacy PDB is unavailable, try mmCIF.
            if (
                status == "failed"
                and args.format == "pdb"
                and args.fallback_cif
            ):
                cif_url = build_download_url(pdb_id, "cif")
                cif_file = output_dir / f"{pdb_id}.cif"

                logger.write(f"{pdb_id}: PDB failed; trying mmCIF fallback")

                cif_status, cif_error, cif_bytes = download_file(
                    url=cif_url,
                    output_path=cif_file,
                    user_agent=args.user_agent,
                    timeout=args.timeout,
                    retries=args.retries,
                    logger=logger,
                )

                if cif_status == "downloaded":
                    status = "downloaded_fallback_cif"
                    error = ""
                    bytes_written = cif_bytes
                    url = cif_url
                    out_file = cif_file
                    extension = "cif"
                else:
                    error = f"PDB failed: {error}; CIF failed: {cif_error}"

        manifest_row = row.to_dict()
        manifest_row.update({
            "download_format": "cif" if status == "downloaded_fallback_cif" else args.format,
            "download_url": url,
            "local_path": str(out_file),
            "download_status": status,
            "download_error": error,
            "bytes": bytes_written,
        })
        rows.append(manifest_row)

        logger.write(f"{pdb_id}: {status}")

        if args.sleep:
            time.sleep(args.sleep)

    manifest = pd.DataFrame(rows)
    manifest.to_csv(manifest_path, index=False)

    status_counts = dict(Counter(manifest["download_status"].astype(str)))

    metadata = {
        "script": Path(__file__).name,
        "inventory": str(inventory_path),
        "output_dir": str(output_dir),
        "manifest": str(manifest_path),
        "format": args.format,
        "unique_pdb_ids": int(len(inv)),
        "status_counts": status_counts,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = manifest_path.with_suffix(manifest_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Status counts: {status_counts}")
    logger.write(f"Saved manifest: {manifest_path}")
    logger.write(f"Saved metadata: {metadata_path}")
    logger.write("Finished PDB structure download")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
