#!/usr/bin/env python3
# Purpose: Fetch RefSeq FASTA sequences.

from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError

import pandas as pd


NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


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
    if text.lower() in {"", "nan", "none", "na", "n/a", "null", "unknown"}:
        return ""
    return text


def safe_filename(value: str) -> str:
    value = clean_text(value)
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def build_efetch_url(refseq_acc: str, email: str | None = None, api_key: str | None = None) -> str:
    params = {
        "db": "protein",
        "id": refseq_acc,
        "rettype": "fasta",
        "retmode": "text",
    }

    if email:
        params["email"] = email

    if api_key:
        params["api_key"] = api_key

    return NCBI_EFETCH + "?" + urllib.parse.urlencode(params)


def fetch_text(
    url: str,
    user_agent: str,
    timeout: int,
    retries: int,
    logger: Logger,
) -> tuple[str, str]:
    """
    Fetch text from URL.

    Returns:
        text, error
    """
    last_error = ""

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": user_agent})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                data = response.read().decode("utf-8", errors="replace")
            return data, ""

        except HTTPError as exc:
            last_error = f"HTTPError {exc.code}: {exc.reason}"
            logger.write(last_error)

        except URLError as exc:
            last_error = f"URLError: {exc}"
            logger.write(last_error)

        except TimeoutError as exc:
            last_error = f"TimeoutError: {exc}"
            logger.write(last_error)

        if attempt < retries:
            time.sleep(min(5 * attempt, 30))

    return "", last_error


def validate_fasta(text: str) -> tuple[bool, int, str]:
    """
    Validate FASTA text.

    Returns:
        valid, sequence_length, error
    """
    text = text.strip()

    if not text:
        return False, 0, "empty response"

    if not text.startswith(">"):
        return False, 0, "response does not start with FASTA header"

    lines = text.splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    seq = re.sub(r"\s+", "", seq)

    if not seq:
        return False, 0, "no sequence found"

    return True, len(seq), ""


def append_to_combined_fasta(handle, fasta_text: str, metadata_prefix: str) -> None:
    lines = fasta_text.strip().splitlines()
    if not lines:
        return

    header = lines[0]
    seq_lines = lines[1:]

    if header.startswith(">"):
        header = ">" + metadata_prefix + " | " + header[1:]

    handle.write(header + "\n")
    for line in seq_lines:
        handle.write(line.rstrip() + "\n")
    handle.write("\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch RefSeq protein FASTA records.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--inventory", required=True, help="Input RefSeq inventory CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for individual FASTA files.")
    parser.add_argument("--combined-fasta", required=True, help="Combined FASTA output path.")
    parser.add_argument("--manifest", required=True, help="Fetch manifest CSV path.")
    parser.add_argument("--log", default=None, help="Log file path.")
    parser.add_argument("--email", default=None, help="Optional email for NCBI E-utilities.")
    parser.add_argument("--api-key", default=None, help="Optional NCBI API key.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing FASTA files.")
    parser.add_argument("--timeout", type=int, default=120, help="HTTP timeout.")
    parser.add_argument("--retries", type=int, default=3, help="HTTP retries.")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between requests.")
    parser.add_argument("--user-agent", default="enzyme-thermostability-pipeline/1.0", help="HTTP User-Agent.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    inventory_path = resolve_path(args.inventory, project_root)
    output_dir = resolve_path(args.output_dir, project_root)
    combined_fasta_path = resolve_path(args.combined_fasta, project_root)
    manifest_path = resolve_path(args.manifest, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "07_fetch_refseq_fasta.log"),
        project_root,
    )

    if inventory_path is None or output_dir is None or combined_fasta_path is None or manifest_path is None:
        raise ValueError("--inventory, --output-dir, --combined-fasta and --manifest are required")

    logger = Logger(log_path)

    logger.write("Starting RefSeq FASTA fetch")
    logger.write(f"Inventory: {inventory_path}")
    logger.write(f"Output directory: {output_dir}")
    logger.write(f"Combined FASTA: {combined_fasta_path}")
    logger.write(f"Manifest: {manifest_path}")

    if not inventory_path.exists():
        raise FileNotFoundError(f"RefSeq inventory not found: {inventory_path}")

    inv = pd.read_csv(inventory_path, dtype=str, keep_default_na=False).fillna("")

    if "refseq_acc" not in inv.columns:
        raise ValueError("Inventory must contain column: refseq_acc")

    inv["refseq_acc"] = inv["refseq_acc"].map(clean_text)
    inv = inv[inv["refseq_acc"].str.len() > 0].copy()
    inv = inv.drop_duplicates(subset=["refseq_acc"], keep="first").copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    combined_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    logger.write(f"Unique RefSeq accessions to fetch: {len(inv)}")

    manifest_rows = []

    with combined_fasta_path.open("w", encoding="utf-8") as combined_handle:
        for _, row in inv.iterrows():
            refseq_acc = clean_text(row.get("refseq_acc", ""))
            uniprot_accession = clean_text(row.get("uniprot_accession", ""))
            organism_type = clean_text(row.get("organism_type", ""))
            gh_family = clean_text(row.get("gh_family", ""))

            filename = safe_filename(refseq_acc) + ".fasta"
            fasta_path = output_dir / filename

            status = ""
            error = ""
            seq_len = 0
            url = build_efetch_url(refseq_acc, email=args.email, api_key=args.api_key)

            if fasta_path.exists() and fasta_path.stat().st_size > 0 and not args.overwrite:
                text = fasta_path.read_text(encoding="utf-8", errors="replace")
                valid, seq_len, error = validate_fasta(text)

                if valid:
                    status = "skipped_existing"
                    metadata_prefix = f"refseq={refseq_acc} | uniprot={uniprot_accession} | family={gh_family} | organism_type={organism_type}"
                    append_to_combined_fasta(combined_handle, text, metadata_prefix)
                else:
                    status = "existing_invalid"

            else:
                text, error = fetch_text(
                    url=url,
                    user_agent=args.user_agent,
                    timeout=args.timeout,
                    retries=args.retries,
                    logger=logger,
                )

                valid, seq_len, validation_error = validate_fasta(text)

                if valid:
                    fasta_path.write_text(text.strip() + "\n", encoding="utf-8")
                    status = "downloaded"
                    metadata_prefix = f"refseq={refseq_acc} | uniprot={uniprot_accession} | family={gh_family} | organism_type={organism_type}"
                    append_to_combined_fasta(combined_handle, text, metadata_prefix)
                    error = ""
                else:
                    status = "failed"
                    error = error or validation_error

            manifest_row = row.to_dict()
            manifest_row.update({
                "fasta_path": str(fasta_path),
                "fetch_url": url,
                "fetch_status": status,
                "fetch_error": error,
                "sequence_length": seq_len,
            })
            manifest_rows.append(manifest_row)

            logger.write(f"{refseq_acc}: {status}")

            if args.sleep:
                time.sleep(args.sleep)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(manifest_path, index=False)

    status_counts = dict(Counter(manifest["fetch_status"].astype(str)))

    metadata = {
        "script": Path(__file__).name,
        "inventory": str(inventory_path),
        "output_dir": str(output_dir),
        "combined_fasta": str(combined_fasta_path),
        "manifest": str(manifest_path),
        "unique_refseq_accessions": int(len(inv)),
        "status_counts": status_counts,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = manifest_path.with_suffix(manifest_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Status counts: {status_counts}")
    logger.write(f"Saved manifest: {manifest_path}")
    logger.write(f"Saved metadata: {metadata_path}")
    logger.write("Finished RefSeq FASTA fetch")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
