#!/usr/bin/env python3
"""
Fetch UniProtKB enzyme records using the UniProt REST streaming API.

This script is intentionally generic. It can be used for xylanase or any other
enzyme if the EC number, UniProt query, fields and output paths are supplied
through a YAML config file or command-line arguments.

Example:
    python scripts/01_curation/01_fetch_uniprot_xylanase_tsv.py \
        --config config/xylanase_config.yaml \
        --output data/raw/uniprot_xylanase.tsv \
        --log logs/01_fetch_uniprot_xylanase_tsv.log

Generic example:
    python scripts/01_curation/01_fetch_uniprot_xylanase_tsv.py \
        --ec 3.2.1.4 \
        --enzyme-name cellulase \
        --output data/raw/uniprot_cellulase.tsv
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError


DEFAULT_UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

DEFAULT_FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "lineage",
    "length",
    "sequence",
    "xref_pdb",
    "xref_refseq",
]


def load_yaml_config(path: str | None) -> dict:
    """Load YAML config if supplied."""
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
    """Safely retrieve nested dictionary values."""
    cur = data
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def resolve_path(path_value: str | None, project_root: str | None = None) -> Path | None:
    """Resolve a possibly relative output/log path."""
    if not path_value:
        return None

    path = Path(os.path.expanduser(path_value))

    if path.is_absolute():
        return path

    if project_root:
        return Path(os.path.expanduser(project_root)) / path

    return Path.cwd() / path


class Logger:
    """Simple console + file logger."""

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


def normalise_fields(fields_value) -> list[str]:
    """Normalise fields from config or CLI."""
    if fields_value is None:
        return DEFAULT_FIELDS

    if isinstance(fields_value, str):
        return [x.strip() for x in fields_value.split(",") if x.strip()]

    if isinstance(fields_value, list):
        return [str(x).strip() for x in fields_value if str(x).strip()]

    raise ValueError("Fields must be either a comma-separated string or a list.")


def build_uniprot_query(
    explicit_query: str | None,
    ec_number: str | None,
    reviewed: str | None,
) -> str:
    """Construct UniProt query from explicit query or EC number."""
    if explicit_query:
        query = explicit_query.strip()
    elif ec_number:
        query = f"ec:{ec_number}"
    else:
        raise ValueError("Either --query or --ec must be provided, or set in config.")

    reviewed_value = (reviewed or "any").lower().strip()

    if reviewed_value in {"true", "yes", "reviewed"}:
        query = f"({query}) AND reviewed:true"
    elif reviewed_value in {"false", "no", "unreviewed"}:
        query = f"({query}) AND reviewed:false"
    elif reviewed_value in {"any", "all", ""}:
        pass
    else:
        raise ValueError(
            "--reviewed must be one of: any, true, false. "
            f"Received: {reviewed}"
        )

    return query


def build_stream_url(
    stream_url: str,
    query: str,
    fields: list[str],
    output_format: str,
) -> str:
    """Build UniProt stream URL."""
    params = {
        "query": query,
        "format": output_format,
        "fields": ",".join(fields),
    }
    return stream_url + "?" + urllib.parse.urlencode(params)


def stream_download(
    url: str,
    output_path: Path,
    logger: Logger,
    user_agent: str,
    timeout: int,
    retries: int,
    chunk_size: int = 1024 * 1024,
) -> int:
    """
    Stream URL content to output_path.

    Returns:
        Number of bytes written.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    last_error = None

    for attempt in range(1, retries + 1):
        logger.write(f"Download attempt {attempt}/{retries}")
        logger.write(f"Request URL: {url}")

        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": user_agent},
            )

            bytes_written = 0

            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", "unknown")
                logger.write(f"HTTP status: {status}")

                with output_path.open("wb") as out_handle:
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        out_handle.write(chunk)
                        bytes_written += len(chunk)

            return bytes_written

        except HTTPError as exc:
            last_error = exc
            error_body = ""
            try:
                error_body = exc.read().decode("utf-8", errors="replace")
            except Exception:
                error_body = "<could not read error body>"

            logger.write(f"HTTPError: {exc.code} {exc.reason}")
            logger.write(f"Error body: {error_body[:2000]}")

        except URLError as exc:
            last_error = exc
            logger.write(f"URLError: {exc}")

        except TimeoutError as exc:
            last_error = exc
            logger.write(f"TimeoutError: {exc}")

        if attempt < retries:
            sleep_seconds = min(10 * attempt, 60)
            logger.write(f"Retrying after {sleep_seconds} seconds...")
            time.sleep(sleep_seconds)

    raise RuntimeError(f"Failed to download UniProt data after {retries} attempts: {last_error}")


def count_tsv_rows(path: Path) -> int:
    """Count data rows in a TSV file, excluding header."""
    if not path.exists():
        return 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        line_count = sum(1 for _ in handle)

    return max(0, line_count - 1)


def write_metadata(
    output_path: Path,
    metadata: dict,
) -> None:
    """Write JSON metadata beside the output file."""
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch UniProtKB enzyme records using UniProt REST streaming API."
    )

    parser.add_argument(
        "--config",
        help="Optional YAML config file.",
        default=None,
    )
    parser.add_argument(
        "--query",
        help="Explicit UniProt query. Overrides EC-derived query.",
        default=None,
    )
    parser.add_argument(
        "--ec",
        help="EC number, for example 3.2.1.8.",
        default=None,
    )
    parser.add_argument(
        "--enzyme-name",
        help="Enzyme name used only for metadata/logging.",
        default=None,
    )
    parser.add_argument(
        "--reviewed",
        choices=["any", "true", "false"],
        help="Filter reviewed UniProt records. Default from config or any.",
        default=None,
    )
    parser.add_argument(
        "--fields",
        help="Comma-separated UniProt fields. Overrides config.",
        default=None,
    )
    parser.add_argument(
        "--format",
        help="UniProt output format.",
        default=None,
    )
    parser.add_argument(
        "--output",
        help="Output TSV path.",
        default=None,
    )
    parser.add_argument(
        "--log",
        help="Log file path.",
        default=None,
    )
    parser.add_argument(
        "--timeout",
        type=int,
        help="HTTP timeout in seconds.",
        default=None,
    )
    parser.add_argument(
        "--retries",
        type=int,
        help="Number of download retries.",
        default=None,
    )
    parser.add_argument(
        "--stream-url",
        help="UniProt stream endpoint.",
        default=None,
    )
    parser.add_argument(
        "--user-agent",
        help="HTTP User-Agent string.",
        default=None,
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    enzyme_name = (
        args.enzyme_name
        or deep_get(config, ["dataset", "enzyme_name"], None)
        or "enzyme"
    )

    ec_number = (
        args.ec
        or deep_get(config, ["dataset", "ec_number"], None)
    )

    reviewed = (
        args.reviewed
        or deep_get(config, ["dataset", "reviewed"], None)
        or "any"
    )

    explicit_query = (
        args.query
        or deep_get(config, ["uniprot", "query"], None)
    )

    fields = normalise_fields(
        args.fields
        or deep_get(config, ["uniprot", "fields"], None)
    )

    output_format = (
        args.format
        or deep_get(config, ["uniprot", "format"], None)
        or "tsv"
    )

    stream_url = (
        args.stream_url
        or deep_get(config, ["uniprot", "stream_url"], None)
        or DEFAULT_UNIPROT_STREAM
    )

    timeout = int(
        args.timeout
        or deep_get(config, ["uniprot", "timeout"], None)
        or 180
    )

    retries = int(
        args.retries
        or deep_get(config, ["uniprot", "retries"], None)
        or 3
    )

    user_agent = (
        args.user_agent
        or deep_get(config, ["uniprot", "user_agent"], None)
        or "enzyme-thermostability-pipeline/1.0"
    )

    output_path = resolve_path(
        args.output
        or deep_get(config, ["outputs", "raw_uniprot_tsv"], None)
        or f"data/raw/uniprot_{enzyme_name}.tsv",
        project_root,
    )

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "01_fetch_uniprot.log"),
        project_root,
    )

    assert output_path is not None
    assert log_path is not None

    logger = Logger(log_path)

    logger.write("Starting UniProt fetch")
    logger.write(f"Enzyme name: {enzyme_name}")
    logger.write(f"EC number: {ec_number}")
    logger.write(f"Reviewed filter: {reviewed}")
    logger.write(f"Fields: {','.join(fields)}")
    logger.write(f"Output: {output_path}")

    query = build_uniprot_query(
        explicit_query=explicit_query,
        ec_number=ec_number,
        reviewed=reviewed,
    )

    url = build_stream_url(
        stream_url=stream_url,
        query=query,
        fields=fields,
        output_format=output_format,
    )

    start = time.time()

    bytes_written = stream_download(
        url=url,
        output_path=output_path,
        logger=logger,
        user_agent=user_agent,
        timeout=timeout,
        retries=retries,
    )

    elapsed = round(time.time() - start, 2)
    row_count = count_tsv_rows(output_path)

    metadata = {
        "script": Path(__file__).name,
        "enzyme_name": enzyme_name,
        "ec_number": ec_number,
        "query": query,
        "fields": fields,
        "format": output_format,
        "stream_url": stream_url,
        "output": str(output_path),
        "rows": row_count,
        "bytes_written": bytes_written,
        "elapsed_seconds": elapsed,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Rows downloaded: {row_count}")
    logger.write(f"Bytes written: {bytes_written}")
    logger.write(f"Elapsed seconds: {elapsed}")
    logger.write("Finished UniProt fetch successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
