#!/usr/bin/env python3
# Purpose: Extract generated MD metrics.
from __future__ import annotations

import argparse
from pathlib import Path
import re
import math
import statistics
import pandas as pd


TEMP_RE = re.compile(r"^\d{3}K$")


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None


def infer_system_and_temp(path: Path):
    parts = list(path.parts)

    temp = ""
    system_id = ""

    for i, p in enumerate(parts):
        if TEMP_RE.match(p):
            temp = p
            if i > 0:
                system_id = parts[i - 1]
            break

    if not system_id:
        system_id = path.parent.name

    return system_id, temp


def classify_metric(path: Path):
    name = path.name.lower()
    stem = path.stem.lower()

    if path.suffix.lower() == ".xpm":
        if "dssp" in name or stem in {"ss", "secondary", "secondary_structure"} or "ss_" in name:
            return "dssp_xpm"

    if path.suffix.lower() != ".xvg":
        return None

    if "rmsf" in name:
        return "rmsf"

    if "rmsd" in name:
        return "rmsd"

    if "gyrate" in name or stem in {"rg", "md_rg"} or "_rg" in stem or "rg_" in stem:
        return "rg"

    if "hbond" in name or "hbnum" in name or "hbonds" in name:
        return "hbond"

    if "sasa" in name or "area" in name:
        return "sasa"

    if "dssp" in name or "scount" in name or "ss_count" in name or "secondary" in name:
        return "dssp"

    return None


def parse_xvg(path: Path):
    rows = []

    with path.open("r", errors="ignore") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#") or line.startswith("@"):
                continue

            parts = line.split()
            values = [safe_float(x) for x in parts]

            if any(v is None for v in values):
                continue

            if len(values) >= 2:
                rows.append(values)

    return rows


def summarise_values(values):
    clean = [v for v in values if v is not None and not math.isnan(v)]

    if not clean:
        return {
            "n_points": 0,
            "mean_all": "",
            "sd_all": "",
            "min_value": "",
            "max_value": "",
            "final_value": "",
            "late50_mean": "",
            "late50_sd": "",
        }

    late = clean[len(clean)//2:]

    return {
        "n_points": len(clean),
        "mean_all": statistics.mean(clean),
        "sd_all": statistics.stdev(clean) if len(clean) > 1 else 0.0,
        "min_value": min(clean),
        "max_value": max(clean),
        "final_value": clean[-1],
        "late50_mean": statistics.mean(late) if late else "",
        "late50_sd": statistics.stdev(late) if len(late) > 1 else 0.0,
    }


def summarise_xvg(path: Path, metric: str, root: Path):
    rows = parse_xvg(path)

    system_id, temperature = infer_system_and_temp(path)

    rel = str(path)
    try:
        rel = str(path.relative_to(Path.cwd()))
    except Exception:
        pass

    if not rows:
        return [{
            "source_root": str(root),
            "system_id": system_id,
            "temperature": temperature,
            "metric": metric,
            "metric_component": "unparsed",
            "source_file": rel,
            "n_numeric_columns": 0,
            "time_start": "",
            "time_end": "",
            "n_points": 0,
            "mean_all": "",
            "sd_all": "",
            "min_value": "",
            "max_value": "",
            "final_value": "",
            "late50_mean": "",
            "late50_sd": "",
        }]

    ncols = max(len(r) for r in rows)
    time_values = [r[0] for r in rows if len(r) >= 1]

    output = []

    if metric in {"rmsd", "rmsf", "rg", "hbond", "sasa"}:
        value_columns = [1]
    else:
        value_columns = list(range(1, ncols))

    for col in value_columns:
        values = [r[col] for r in rows if len(r) > col]
        summary = summarise_values(values)

        component = metric
        if metric == "dssp":
            component = f"dssp_col{col}"

        output.append({
            "source_root": str(root),
            "system_id": system_id,
            "temperature": temperature,
            "metric": metric,
            "metric_component": component,
            "source_file": rel,
            "n_numeric_columns": ncols,
            "time_start": min(time_values) if time_values else "",
            "time_end": max(time_values) if time_values else "",
            **summary,
        })

    return output


def parse_xpm_dssp(path: Path, root: Path):
    system_id, temperature = infer_system_and_temp(path)

    rel = str(path)
    try:
        rel = str(path.relative_to(Path.cwd()))
    except Exception:
        pass

    text = path.read_text(errors="ignore").splitlines()

    symbol_map = {}
    matrix_lines = []

    for line in text:
        line = line.strip()

        legend = re.search(r'"(.{1})\s+c\s+[^"]*"\s*/\*\s*"([^"]+)"\s*\*/', line)
        if legend:
            symbol_map[legend.group(1)] = legend.group(2)

        if line.startswith('"') and not line.startswith('"/*'):
            content = line.strip().strip(",").strip('"')
            if content and not content.startswith("/*"):
                matrix_lines.append(content)

    counts = {}
    total = 0

    for row in matrix_lines:
        for ch in row:
            label = symbol_map.get(ch, ch)
            counts[label] = counts.get(label, 0) + 1
            total += 1

    output = []

    if total == 0:
        output.append({
            "source_root": str(root),
            "system_id": system_id,
            "temperature": temperature,
            "metric": "dssp_xpm",
            "metric_component": "unparsed",
            "source_file": rel,
            "n_numeric_columns": "",
            "time_start": "",
            "time_end": "",
            "n_points": 0,
            "mean_all": "",
            "sd_all": "",
            "min_value": "",
            "max_value": "",
            "final_value": "",
            "late50_mean": "",
            "late50_sd": "",
        })
        return output

    for label, count in sorted(counts.items()):
        output.append({
            "source_root": str(root),
            "system_id": system_id,
            "temperature": temperature,
            "metric": "dssp_xpm",
            "metric_component": f"dssp_fraction_{label}",
            "source_file": rel,
            "n_numeric_columns": "",
            "time_start": "",
            "time_end": "",
            "n_points": total,
            "mean_all": count / total,
            "sd_all": "",
            "min_value": "",
            "max_value": "",
            "final_value": "",
            "late50_mean": "",
            "late50_sd": "",
        })

    return output


def file_priority(row):
    path = str(row["source_file"]).lower()
    score = 0

    if "/analysis/" in path:
        score += 100000

    if "protein" in path:
        score += 10000

    if "backbone" in path:
        score += 5000

    if "rmsd_protein" in path or "rmsd_backbone" in path:
        score += 1000

    if "gyrate_protein" in path or "rg_protein" in path:
        score += 1000

    if "hbnum_protein" in path or "hbond_protein" in path:
        score += 1000

    try:
        score += int(row.get("n_points", 0))
    except Exception:
        pass

    return score


def make_primary_and_wide(long_df: pd.DataFrame):
    if long_df.empty:
        return long_df, pd.DataFrame()

    df = long_df.copy()
    df["selection_score"] = df.apply(file_priority, axis=1)

    sort_cols = [
        "system_id",
        "temperature",
        "metric_component",
        "selection_score",
        "n_points",
    ]

    df = df.sort_values(sort_cols, ascending=[True, True, True, False, False])

    primary = df.drop_duplicates(
        subset=["system_id", "temperature", "metric_component"],
        keep="first"
    ).copy()

    compact_metrics = primary[
        primary["metric_component"].isin(["rmsd", "rmsf", "rg", "hbond", "sasa"])
    ].copy()

    if compact_metrics.empty:
        return primary, pd.DataFrame()

    wide = compact_metrics.pivot_table(
        index=["system_id", "temperature"],
        columns="metric_component",
        values=["mean_all", "sd_all", "late50_mean", "final_value", "source_file"],
        aggfunc="first"
    )

    wide.columns = [f"{metric}_{value}" for value, metric in wide.columns]
    wide = wide.reset_index()

    return primary, wide


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--roots",
        nargs="+",
        default=[
            "md_10ns/systems",
            "md_tier2_wt_mutant_compact/systems",
            "md_tier2_wt_mutant_compact/analysis",
        ],
        help="Folders to scan recursively."
    )
    parser.add_argument("--outdir", default="results/md_generated_metric_extraction")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    inventory_rows = []

    for root_str in args.roots:
        root = Path(root_str)

        if not root.exists():
            inventory_rows.append({
                "root": str(root),
                "status": "missing",
                "n_candidate_files": 0,
            })
            continue

        files = []
        for suffix in ("*.xvg", "*.xpm"):
            files.extend(root.rglob(suffix))

        n_used = 0

        for path in sorted(files):
            metric = classify_metric(path)

            if metric is None:
                continue

            n_used += 1

            if metric == "dssp_xpm":
                all_rows.extend(parse_xpm_dssp(path, root))
            else:
                all_rows.extend(summarise_xvg(path, metric, root))

        inventory_rows.append({
            "root": str(root),
            "status": "scanned",
            "n_candidate_files": len(files),
            "n_metric_files_used": n_used,
        })

    long_df = pd.DataFrame(all_rows)
    inventory_df = pd.DataFrame(inventory_rows)

    if not long_df.empty:
        long_df = long_df.sort_values(
            ["system_id", "temperature", "metric", "metric_component", "source_file"]
        )

    primary_df, wide_df = make_primary_and_wide(long_df)

    long_csv = outdir / "all_md_generated_metrics_long.csv"
    primary_csv = outdir / "all_md_generated_metrics_primary_selected.csv"
    wide_csv = outdir / "all_md_generated_metrics_compact_wide.csv"
    inventory_csv = outdir / "md_generated_metric_file_inventory.csv"

    long_df.to_csv(long_csv, index=False)
    primary_df.to_csv(primary_csv, index=False)
    wide_df.to_csv(wide_csv, index=False)
    inventory_df.to_csv(inventory_csv, index=False)

    long_md = outdir / "all_md_generated_metrics_primary_selected_preview.md"
    wide_md = outdir / "all_md_generated_metrics_compact_wide_preview.md"

    primary_df.head(80).to_markdown(long_md, index=False)
    wide_df.head(80).to_markdown(wide_md, index=False)

    print("Saved:")
    print(long_csv)
    print(primary_csv)
    print(wide_csv)
    print(inventory_csv)
    print(long_md)
    print(wide_md)

    print()
    print("Inventory:")
    print(inventory_df.to_string(index=False))

    print()
    print("Metric counts:")
    if not long_df.empty:
        print(long_df["metric"].value_counts(dropna=False).to_string())
    else:
        print("No metric files found.")


if __name__ == "__main__":
    main()
