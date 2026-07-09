#!/usr/bin/env python3
# Purpose: Run TM-align comparisons.

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import pandas as pd


RESULT_COLUMNS = [
    "job_id",
    "query_uniprot_accession",
    "query_structure_id",
    "query_structure_source",
    "query_organism_type",
    "query_gh_family",
    "query_structure_path",
    "reference_uniprot_accession",
    "reference_structure_id",
    "reference_structure_source",
    "reference_organism_type",
    "reference_gh_family",
    "reference_structure_path",
    "same_gh_family",
    "same_organism_type",
    "same_uniprot_accession",
    "tmalign_status",
    "return_code",
    "aligned_length",
    "rmsd",
    "seq_id_aligned",
    "tm_score_query",
    "tm_score_reference",
    "tm_score_max",
    "tm_score_min",
    "runtime_seconds",
    "error_message",
]


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


def find_tmalign_binary(user_value: str | None) -> str:
    candidates = []

    if user_value:
        candidates.append(user_value)

    candidates.extend(["TMalign", "TM-align", "tmalign", "TMalign.exe"])

    for candidate in candidates:
        if not candidate:
            continue

        expanded = os.path.expanduser(candidate)

        if Path(expanded).exists():
            return str(Path(expanded).resolve())

        resolved = shutil.which(candidate)
        if resolved:
            return resolved

    raise FileNotFoundError(
        "TM-align binary not found. Provide --tmalign-bin /path/to/TMalign"
    )


def parse_tmalign_output(stdout: str) -> dict:
    result = {
        "aligned_length": "",
        "rmsd": "",
        "seq_id_aligned": "",
        "tm_score_query": "",
        "tm_score_reference": "",
        "tm_score_max": "",
        "tm_score_min": "",
    }

    aligned_match = re.search(
        r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+),\s*Seq_ID=n_identical/n_aligned=\s*([0-9.]+)",
        stdout,
    )

    if aligned_match:
        result["aligned_length"] = aligned_match.group(1)
        result["rmsd"] = aligned_match.group(2)
        result["seq_id_aligned"] = aligned_match.group(3)

    tm_scores = re.findall(r"TM-score=\s*([0-9.]+)", stdout)

    if len(tm_scores) >= 1:
        result["tm_score_query"] = tm_scores[0]

    if len(tm_scores) >= 2:
        result["tm_score_reference"] = tm_scores[1]

    numeric_scores = []

    for value in tm_scores[:2]:
        try:
            numeric_scores.append(float(value))
        except Exception:
            pass

    if numeric_scores:
        result["tm_score_max"] = max(numeric_scores)
        result["tm_score_min"] = min(numeric_scores)

    return result


def row_to_result_base(row: dict) -> dict:
    base = {col: "" for col in RESULT_COLUMNS}

    for col in RESULT_COLUMNS:
        if col in row:
            base[col] = clean_text(row.get(col, ""))

    return base


def run_one_job(row: dict, tmalign_bin: str, timeout: int) -> dict:
    start = time.time()

    result = row_to_result_base(row)

    query_path = clean_text(row.get("query_structure_path", ""))
    reference_path = clean_text(row.get("reference_structure_path", ""))

    if not query_path or not Path(query_path).exists():
        result["tmalign_status"] = "failed_missing_query"
        result["error_message"] = f"Query file missing: {query_path}"
        result["runtime_seconds"] = round(time.time() - start, 4)
        return result

    if not reference_path or not Path(reference_path).exists():
        result["tmalign_status"] = "failed_missing_reference"
        result["error_message"] = f"Reference file missing: {reference_path}"
        result["runtime_seconds"] = round(time.time() - start, 4)
        return result

    cmd = [tmalign_bin, query_path, reference_path]

    try:
        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            timeout=timeout,
        )

        result["return_code"] = proc.returncode

        parsed = parse_tmalign_output(proc.stdout)

        for key, value in parsed.items():
            result[key] = value

        if proc.returncode != 0:
            result["tmalign_status"] = "failed_return_code"
            result["error_message"] = clean_text(proc.stderr)[:500]
        elif not result["tm_score_max"]:
            result["tmalign_status"] = "failed_parse"
            result["error_message"] = "TM-align completed but TM-score was not parsed"
        else:
            result["tmalign_status"] = "success"

    except subprocess.TimeoutExpired:
        result["tmalign_status"] = "failed_timeout"
        result["error_message"] = f"Timed out after {timeout} seconds"

    except Exception as exc:
        result["tmalign_status"] = "failed_exception"
        result["error_message"] = str(exc)

    result["runtime_seconds"] = round(time.time() - start, 4)

    return result


def load_completed_job_ids(output_path: Path) -> set[str]:
    if not output_path.exists() or output_path.stat().st_size == 0:
        return set()

    try:
        df = pd.read_csv(output_path, dtype=str, usecols=["job_id", "tmalign_status"]).fillna("")
        completed = set(df[df["tmalign_status"] == "success"]["job_id"].astype(str))
        return completed
    except Exception:
        return set()


def append_results(output_path: Path, rows: list[dict]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    file_exists = output_path.exists() and output_path.stat().st_size > 0

    with output_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS)

        if not file_exists:
            writer.writeheader()

        for row in rows:
            writer.writerow({col: row.get(col, "") for col in RESULT_COLUMNS})


def write_failed_output(output_path: Path, failed_path: Path) -> None:
    if not output_path.exists() or output_path.stat().st_size == 0:
        pd.DataFrame(columns=RESULT_COLUMNS).to_csv(failed_path, index=False)
        return

    df = pd.read_csv(output_path, dtype=str).fillna("")
    failed = df[df["tmalign_status"] != "success"].copy()
    failed_path.parent.mkdir(parents=True, exist_ok=True)
    failed.to_csv(failed_path, index=False)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run TM-align jobs.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--jobs", required=True, help="TM-align job manifest CSV.")
    parser.add_argument("--output", required=True, help="TM-align parsed result CSV.")
    parser.add_argument("--failed-output", required=True, help="Failed TM-align jobs CSV.")
    parser.add_argument("--tmalign-bin", default=None, help="Path/name of TM-align binary.")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers.")
    parser.add_argument("--timeout", type=int, default=120, help="Timeout per job in seconds.")
    parser.add_argument("--limit", type=int, default=0, help="Optional pilot limit; 0 means all jobs.")
    parser.add_argument("--resume", action="store_true", help="Skip already successful job IDs in output.")
    parser.add_argument("--flush-every", type=int, default=100, help="Write results every N completed jobs.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    jobs_path = resolve_path(args.jobs, project_root)
    output_path = resolve_path(args.output, project_root)
    failed_output_path = resolve_path(args.failed_output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "04_structure_05_run_tmalign.log"),
        project_root,
    )

    logger = Logger(log_path)

    tmalign_bin = find_tmalign_binary(args.tmalign_bin)

    logger.write("Starting TM-align run")
    logger.write(f"Jobs: {jobs_path}")
    logger.write(f"Output: {output_path}")
    logger.write(f"Failed output: {failed_output_path}")
    logger.write(f"TM-align binary: {tmalign_bin}")
    logger.write(f"Workers: {args.workers}")
    logger.write(f"Timeout: {args.timeout}")
    logger.write(f"Limit: {args.limit if args.limit else 'all'}")
    logger.write(f"Resume: {args.resume}")

    if jobs_path is None or not jobs_path.exists():
        raise FileNotFoundError(f"Jobs file not found: {jobs_path}")

    jobs = pd.read_csv(jobs_path, dtype=str, keep_default_na=False).fillna("")

    if args.limit and args.limit > 0:
        jobs = jobs.head(args.limit).copy()

    completed = set()

    if args.resume:
        completed = load_completed_job_ids(output_path)
        if completed:
            jobs = jobs[~jobs["job_id"].astype(str).isin(completed)].copy()

    logger.write(f"Jobs to run: {len(jobs)}")
    logger.write(f"Already successful jobs skipped: {len(completed)}")

    if len(jobs) == 0:
        write_failed_output(output_path, failed_output_path)
        logger.write("No jobs to run.")
        return 0

    buffer = []
    total = len(jobs)
    success = 0
    failed = 0
    processed = 0

    job_dicts = jobs.to_dict(orient="records")

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [
            executor.submit(run_one_job, row, tmalign_bin, args.timeout)
            for row in job_dicts
        ]

        for future in as_completed(futures):
            result = future.result()
            buffer.append(result)
            processed += 1

            if result.get("tmalign_status") == "success":
                success += 1
            else:
                failed += 1

            if len(buffer) >= args.flush_every:
                append_results(output_path, buffer)
                buffer = []

            if processed % 500 == 0 or processed == total:
                logger.write(f"Processed {processed}/{total}; success={success}; failed={failed}")

    if buffer:
        append_results(output_path, buffer)

    write_failed_output(output_path, failed_output_path)

    metadata = {
        "script": Path(__file__).name,
        "jobs": str(jobs_path),
        "output": str(output_path),
        "failed_output": str(failed_output_path),
        "tmalign_bin": tmalign_bin,
        "workers": args.workers,
        "timeout": args.timeout,
        "limit": args.limit if args.limit else None,
        "resume": args.resume,
        "jobs_run_this_call": int(total),
        "success_this_call": int(success),
        "failed_this_call": int(failed),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Finished TM-align run; success={success}; failed={failed}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
