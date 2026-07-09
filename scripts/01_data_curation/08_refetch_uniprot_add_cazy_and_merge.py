#!/usr/bin/env python3
"""
Enrich a curated enzyme master table with CAZy/GH-family annotations.

This script is generic. It can be used for xylanase or other enzyme datasets.
It can reuse an existing UniProt CAZy TSV, refetch CAZy xrefs from UniProt,
merge a local CAZy annotation CSV, and improve family labels without destroying
existing annotations.

Main inputs:
    --input             Curated UniProt CSV from script 02
    --cazy-tsv          Optional UniProt TSV containing Entry + CAZy/xref_cazy
    --annotation-csv    Optional CSV with uniprot_accession, gh_family, cazy_family
    --output            Enriched output CSV

Example:
    python scripts/01_curation/08_refetch_uniprot_add_cazy_and_merge.py \
        --config config/xylanase_config.yaml \
        --input data/curated/uniprot_xylanase_curated.csv \
        --cazy-tsv data/raw/uniprot_xylanase_with_cazy.tsv \
        --annotation-csv results/cazy_annotation.csv \
        --output data/curated/uniprot_xylanase_curated_with_cazy.csv
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
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError

import pandas as pd


DEFAULT_UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"


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


def clean_text(value) -> str:
    if value is None:
        return ""
    text = str(value)
    text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def is_unknown(value) -> bool:
    text = clean_text(value).lower()
    return text in {"", "unknown", "nan", "none", "na", "n/a", "null"}


def normalise_accession(value) -> str:
    return clean_text(value).split("-")[0]


def split_xrefs(value: str) -> list[str]:
    value = clean_text(value)
    if not value:
        return []

    parts = re.split(r"[;,|]\s*", value)
    return [p.strip() for p in parts if p.strip()]


def family_from_text(text: str, accepted_families: set[str] | None = None) -> str:
    """
    Extract a GH family from CAZy-like text.

    Examples:
        GH10 -> GH10
        GH11 -> GH11
        GH10;CBM2 -> GH10
    """
    text = clean_text(text)
    if not text:
        return "unknown"

    matches = re.findall(r"\bGH\s*[-_ ]?\s*(\d+)\b", text, flags=re.IGNORECASE)

    for m in matches:
        family = f"GH{m}"
        if accepted_families is None or family in accepted_families:
            return family

    return "unknown"


def valid_family(value: str, accepted_families: set[str] | None = None) -> bool:
    value = clean_text(value)
    if is_unknown(value):
        return False

    if accepted_families is None:
        return bool(re.match(r"^GH\d+$", value, flags=re.IGNORECASE))

    return value in accepted_families


def read_csv_or_empty(path: Path | None, logger: Logger, kind: str) -> pd.DataFrame:
    if path is None:
        logger.write(f"No {kind} path supplied")
        return pd.DataFrame()

    if not path.exists():
        logger.write(f"{kind} not found: {path}")
        return pd.DataFrame()

    logger.write(f"Reading {kind}: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=False).fillna("")


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    if df.empty:
        return None

    for c in candidates:
        if c in df.columns:
            return c

    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    return None


def fetch_uniprot_cazy_tsv(
    output_tsv: Path,
    query: str,
    stream_url: str,
    user_agent: str,
    timeout: int,
    retries: int,
    logger: Logger,
) -> None:
    """
    Fetch accession + CAZy xrefs from UniProt.

    If UniProt changes field names, this will fail clearly rather than silently.
    """
    fields = ["accession", "xref_cazy"]
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
    }
    url = stream_url + "?" + urllib.parse.urlencode(params)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    last_error = None

    for attempt in range(1, retries + 1):
        logger.write(f"Fetching UniProt CAZy TSV attempt {attempt}/{retries}")
        logger.write(f"Request URL: {url}")

        try:
            request = urllib.request.Request(url, headers={"User-Agent": user_agent})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                data = response.read()

            with output_tsv.open("wb") as handle:
                handle.write(data)

            logger.write(f"Wrote UniProt CAZy TSV: {output_tsv} bytes={len(data)}")
            return

        except HTTPError as exc:
            last_error = exc
            try:
                body = exc.read().decode("utf-8", errors="replace")
            except Exception:
                body = "<could not read body>"
            logger.write(f"HTTPError {exc.code}: {exc.reason}")
            logger.write(body[:2000])

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

    raise RuntimeError(f"Failed to fetch UniProt CAZy TSV after {retries} attempts: {last_error}")


def parse_uniprot_cazy_tsv(cazy_tsv: Path | None, logger: Logger) -> pd.DataFrame:
    """Read UniProt accession + CAZy xref table."""
    df = read_csv_or_empty(cazy_tsv, logger, "UniProt CAZy TSV")
    if df.empty:
        return pd.DataFrame(columns=["uniprot_accession", "xref_cazy_from_uniprot"])

    # TSV, not CSV; reread properly if needed
    try:
        df = pd.read_csv(cazy_tsv, sep="\t", dtype=str, keep_default_na=False).fillna("")
    except Exception:
        df = pd.read_csv(cazy_tsv, dtype=str, keep_default_na=False).fillna("")

    acc_col = find_column(df, ["Entry", "accession", "uniprot_accession"])
    cazy_col = find_column(df, ["CAZy", "xref_cazy", "Cazy", "Cross-reference (CAZy)"])

    if acc_col is None:
        raise ValueError(f"Could not detect accession column in CAZy TSV: {cazy_tsv}")

    if cazy_col is None:
        # If UniProt returns Entry and one extra column, use that extra column.
        other_cols = [c for c in df.columns if c != acc_col]
        if len(other_cols) == 1:
            cazy_col = other_cols[0]
        else:
            raise ValueError(
                f"Could not detect CAZy column in {cazy_tsv}. Columns: {list(df.columns)}"
            )

    out = pd.DataFrame({
        "uniprot_accession": df[acc_col].map(normalise_accession),
        "xref_cazy_from_uniprot": df[cazy_col].map(clean_text),
    })

    out = out.drop_duplicates(subset=["uniprot_accession"], keep="first")
    return out


def parse_annotation_csv(annotation_csv: Path | None, logger: Logger) -> pd.DataFrame:
    """Read local CAZy annotation CSV with uniprot_accession, gh_family, cazy_family."""
    df = read_csv_or_empty(annotation_csv, logger, "local CAZy annotation CSV")
    if df.empty:
        return pd.DataFrame(columns=["uniprot_accession", "gh_family_from_annotation", "cazy_family_from_annotation"])

    acc_col = find_column(df, ["uniprot_accession", "accession", "Entry"])
    gh_col = find_column(df, ["gh_family", "family", "GH_family"])
    cazy_col = find_column(df, ["cazy_family", "xref_cazy", "CAZy"])

    if acc_col is None:
        raise ValueError(f"Could not detect accession column in annotation CSV: {annotation_csv}")

    out = pd.DataFrame()
    out["uniprot_accession"] = df[acc_col].map(normalise_accession)

    if gh_col:
        out["gh_family_from_annotation"] = df[gh_col].map(clean_text)
    else:
        out["gh_family_from_annotation"] = ""

    if cazy_col:
        out["cazy_family_from_annotation"] = df[cazy_col].map(clean_text)
    else:
        out["cazy_family_from_annotation"] = ""

    out = out.drop_duplicates(subset=["uniprot_accession"], keep="first")
    return out


def enrich_master(
    master: pd.DataFrame,
    cazy_from_uniprot: pd.DataFrame,
    annotation: pd.DataFrame,
    accepted_families: set[str] | None,
) -> pd.DataFrame:
    """Merge CAZy annotations and improve gh_family labels."""
    if "uniprot_accession" not in master.columns:
        raise ValueError("Input master must contain column: uniprot_accession")

    enriched = master.copy()
    enriched["uniprot_accession"] = enriched["uniprot_accession"].map(normalise_accession)

    if "gh_family" not in enriched.columns:
        enriched["gh_family"] = "unknown"

    if "xref_cazy" not in enriched.columns:
        enriched["xref_cazy"] = ""

    enriched["gh_family_before_cazy_enrichment"] = enriched["gh_family"].fillna("").astype(str)
    enriched["xref_cazy_before_enrichment"] = enriched["xref_cazy"].fillna("").astype(str)

    if not cazy_from_uniprot.empty:
        enriched = enriched.merge(cazy_from_uniprot, on="uniprot_accession", how="left")
    else:
        enriched["xref_cazy_from_uniprot"] = ""

    if not annotation.empty:
        enriched = enriched.merge(annotation, on="uniprot_accession", how="left")
    else:
        enriched["gh_family_from_annotation"] = ""
        enriched["cazy_family_from_annotation"] = ""

    for col in ["xref_cazy_from_uniprot", "gh_family_from_annotation", "cazy_family_from_annotation"]:
        if col not in enriched.columns:
            enriched[col] = ""
        enriched[col] = enriched[col].fillna("").astype(str)

    # Build final xref_cazy
    final_cazy = []
    cazy_source = []

    for _, row in enriched.iterrows():
        existing = clean_text(row.get("xref_cazy_before_enrichment", ""))
        fetched = clean_text(row.get("xref_cazy_from_uniprot", ""))
        local = clean_text(row.get("cazy_family_from_annotation", ""))

        if existing:
            final_cazy.append(existing)
            cazy_source.append("existing")
        elif fetched:
            final_cazy.append(fetched)
            cazy_source.append("uniprot_refetch")
        elif local:
            final_cazy.append(local)
            cazy_source.append("local_annotation")
        else:
            final_cazy.append("")
            cazy_source.append("missing")

    enriched["xref_cazy"] = final_cazy
    enriched["xref_cazy_source"] = cazy_source
    enriched["has_cazy_xref"] = enriched["xref_cazy"].map(lambda x: bool(clean_text(x)))

    # Build final GH family
    final_family = []
    family_source = []

    for _, row in enriched.iterrows():
        old = clean_text(row.get("gh_family_before_cazy_enrichment", ""))
        ann = clean_text(row.get("gh_family_from_annotation", ""))
        cazy_local = clean_text(row.get("cazy_family_from_annotation", ""))
        cazy_final = clean_text(row.get("xref_cazy", ""))

        family_from_ann = ann if valid_family(ann, accepted_families) else "unknown"
        family_from_local_cazy = family_from_text(cazy_local, accepted_families)
        family_from_final_cazy = family_from_text(cazy_final, accepted_families)

        if valid_family(old, accepted_families):
            final_family.append(old)
            family_source.append("existing")
        elif valid_family(family_from_ann, accepted_families):
            final_family.append(family_from_ann)
            family_source.append("local_annotation_gh_family")
        elif valid_family(family_from_local_cazy, accepted_families):
            final_family.append(family_from_local_cazy)
            family_source.append("local_annotation_cazy_family")
        elif valid_family(family_from_final_cazy, accepted_families):
            final_family.append(family_from_final_cazy)
            family_source.append("xref_cazy")
        else:
            final_family.append("unknown")
            family_source.append("unknown")

    enriched["gh_family"] = final_family
    enriched["gh_family_source"] = family_source
    enriched["gh_family_changed_by_cazy_enrichment"] = (
        enriched["gh_family"].astype(str) != enriched["gh_family_before_cazy_enrichment"].astype(str)
    )

    enriched["cazy_status"] = enriched["has_cazy_xref"].map(lambda x: "has_cazy_xref" if x else "no_cazy_xref")

    return enriched


def summarize(df: pd.DataFrame) -> dict:
    return {
        "rows": int(len(df)),
        "family_counts": dict(Counter(df.get("gh_family", pd.Series(dtype=str)).fillna("").astype(str))),
        "family_source_counts": dict(Counter(df.get("gh_family_source", pd.Series(dtype=str)).fillna("").astype(str))),
        "cazy_status_counts": dict(Counter(df.get("cazy_status", pd.Series(dtype=str)).fillna("").astype(str))),
        "changed_family_count": int(df.get("gh_family_changed_by_cazy_enrichment", pd.Series(dtype=bool)).sum()),
    }


def write_metadata(output_path: Path, metadata: dict) -> None:
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Enrich curated enzyme table with CAZy/GH-family annotations."
    )

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--input", required=True, help="Input curated CSV.")
    parser.add_argument("--output", required=True, help="Output enriched CSV.")
    parser.add_argument("--cazy-tsv", default=None, help="UniProt CAZy TSV path.")
    parser.add_argument("--annotation-csv", default=None, help="Local CAZy annotation CSV.")
    parser.add_argument("--log", default=None, help="Log file path.")
    parser.add_argument("--force-fetch", action="store_true", help="Refetch UniProt CAZy TSV even if --cazy-tsv exists.")
    parser.add_argument("--no-fetch", action="store_true", help="Do not fetch UniProt CAZy TSV if missing.")
    parser.add_argument("--query", default=None, help="UniProt query for CAZy refetch.")
    parser.add_argument("--timeout", type=int, default=None, help="HTTP timeout.")
    parser.add_argument("--retries", type=int, default=None, help="HTTP retries.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    input_path = resolve_path(args.input, project_root)
    output_path = resolve_path(args.output, project_root)

    cazy_tsv = resolve_path(
        args.cazy_tsv
        or deep_get(config, ["outputs", "uniprot_cazy_tsv"], None)
        or "data/raw/uniprot_xylanase_with_cazy.tsv",
        project_root,
    )

    annotation_csv = resolve_path(
        args.annotation_csv
        or deep_get(config, ["outputs", "cazy_annotation_csv"], None)
        or "results/cazy_annotation.csv",
        project_root,
    )

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "08_refetch_uniprot_add_cazy_and_merge.log"),
        project_root,
    )

    logger = Logger(log_path)

    if input_path is None or output_path is None:
        raise ValueError("--input and --output are required")

    logger.write("Starting CAZy/GH-family enrichment")
    logger.write(f"Input: {input_path}")
    logger.write(f"Output: {output_path}")
    logger.write(f"CAZy TSV: {cazy_tsv}")
    logger.write(f"Annotation CSV: {annotation_csv}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    stream_url = deep_get(config, ["uniprot", "stream_url"], DEFAULT_UNIPROT_STREAM)
    query = args.query or deep_get(config, ["uniprot", "query"], None) or deep_get(config, ["dataset", "ec_number"], None)
    if query and not str(query).startswith("ec:") and re.match(r"^\d+\.\d+\.\d+\.\d+$", str(query)):
        query = f"ec:{query}"
    if not query:
        query = "ec:3.2.1.8"

    user_agent = deep_get(config, ["uniprot", "user_agent"], "enzyme-thermostability-pipeline/1.0")
    timeout = int(args.timeout or deep_get(config, ["uniprot", "timeout"], 180))
    retries = int(args.retries or deep_get(config, ["uniprot", "retries"], 3))

    accepted_families_raw = deep_get(config, ["dataset", "families"], None)
    accepted_families = set(accepted_families_raw) if accepted_families_raw else None

    # Fetch CAZy TSV if requested or missing
    if cazy_tsv is not None:
        if args.force_fetch or (not cazy_tsv.exists() and not args.no_fetch):
            fetch_uniprot_cazy_tsv(
                output_tsv=cazy_tsv,
                query=query,
                stream_url=stream_url,
                user_agent=user_agent,
                timeout=timeout,
                retries=retries,
                logger=logger,
            )
        elif cazy_tsv.exists():
            logger.write(f"Using existing CAZy TSV: {cazy_tsv}")
        elif args.no_fetch:
            logger.write("CAZy TSV missing and --no-fetch was used; continuing without UniProt CAZy TSV")

    master = pd.read_csv(input_path, dtype=str, keep_default_na=False).fillna("")
    logger.write(f"Input rows: {len(master)}")
    before_summary = summarize(master)

    cazy_from_uniprot = parse_uniprot_cazy_tsv(cazy_tsv, logger)
    logger.write(f"UniProt CAZy rows: {len(cazy_from_uniprot)}")

    annotation = parse_annotation_csv(annotation_csv, logger)
    logger.write(f"Local annotation rows: {len(annotation)}")

    enriched = enrich_master(
        master=master,
        cazy_from_uniprot=cazy_from_uniprot,
        annotation=annotation,
        accepted_families=accepted_families,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(output_path, index=False)

    after_summary = summarize(enriched)

    metadata = {
        "script": Path(__file__).name,
        "input": str(input_path),
        "output": str(output_path),
        "cazy_tsv": str(cazy_tsv) if cazy_tsv else None,
        "annotation_csv": str(annotation_csv) if annotation_csv else None,
        "query": query,
        "accepted_families": sorted(accepted_families) if accepted_families else None,
        "before_summary": before_summary,
        "after_summary": after_summary,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    write_metadata(output_path, metadata)

    logger.write(f"Output rows: {len(enriched)}")
    logger.write(f"Before family counts: {before_summary.get('family_counts')}")
    logger.write(f"After family counts: {after_summary.get('family_counts')}")
    logger.write(f"Family source counts: {after_summary.get('family_source_counts')}")
    logger.write(f"Changed family count: {after_summary.get('changed_family_count')}")
    logger.write("Finished CAZy/GH-family enrichment successfully")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
