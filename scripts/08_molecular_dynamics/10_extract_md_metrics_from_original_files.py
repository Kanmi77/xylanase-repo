#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import math
import re
import statistics as stats

BASE = Path.home() / "xylanase-thesis"

OUT = BASE / "results" / "md_metrics_extracted_from_original_files"
PDB_OUT = OUT / "phase1_pdb_referenced_baseline"
MOD_OUT = OUT / "phase2_modeller_wt_mutant_validation"
MAN_OUT = OUT / "manifests"

for d in [PDB_OUT, MOD_OUT, MAN_OUT]:
    d.mkdir(parents=True, exist_ok=True)

PDB_IDS = ["1O8S", "1T6G", "1VBU", "2C79", "3NIY", "7K4X", "1VBR", "7K4P", "8RD5"]
PDB_TEMPS = ["333K", "353K", "373K"]

MOD_PROTEINS = ["D5UGW9", "Q60041", "B0ZSE5", "Q9HGX1", "Q9HEN6", "E3WF08"]
MOD_TEMPS = ["333K", "373K"]

TIME_SERIES_PATTERNS = {
    "rmsd": ["*rmsd*.xvg"],
    "rg": ["*rg*.xvg", "*gyrate*.xvg"],
    "hbond": ["*hbond*.xvg", "*hbnum*.xvg"],
    "sasa": ["*sasa*.xvg"],
}

RMSF_PATTERNS = ["*rmsf*.xvg"]
SS_PATTERNS = ["ss_counts.xvg", "*ss_counts*.xvg"]


def is_gromacs_backup_file(path):
    name = path.name
    return name.startswith("#") and name.endswith("#")


def read_xvg(path):
    rows = []
    comments = []
    with open(path, "r", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("@"):
                comments.append(line)
                continue
            parts = line.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) >= 2:
                rows.append(vals)
    return rows, comments


def mean_sd(values):
    values = [v for v in values if v is not None and not math.isnan(v)]
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return stats.mean(values), stats.stdev(values)


def summarize_timeseries(path, phase, system_id, protein, state, temperature, metric, source_root):
    rows, comments = read_xvg(path)
    if not rows:
        return None, []

    values = [r[1] for r in rows if len(r) >= 2]
    times = [r[0] for r in rows if len(r) >= 2]

    if not values:
        return None, []

    mean_all, sd_all = mean_sd(values)
    late_start_index = len(values) // 2
    late_values = values[late_start_index:]
    late_mean, late_sd = mean_sd(late_values)

    summary = {
        "phase": phase,
        "system_id": system_id,
        "protein": protein,
        "state": state,
        "temperature": temperature,
        "metric": metric,
        "source_root": source_root,
        "source_file": str(path.relative_to(BASE)),
        "n_points": len(values),
        "time_start_ps_or_x": times[0] if times else None,
        "time_end_ps_or_x": times[-1] if times else None,
        "mean_all": mean_all,
        "sd_all": sd_all,
        "min_value": min(values),
        "max_value": max(values),
        "final_value": values[-1],
        "late50_mean": late_mean,
        "late50_sd": late_sd,
        "late50_start_index": late_start_index,
        "n_numeric_columns": max(len(r) for r in rows),
    }

    long_rows = []
    for r in rows:
        if len(r) >= 2:
            long_rows.append({
                "phase": phase,
                "system_id": system_id,
                "protein": protein,
                "state": state,
                "temperature": temperature,
                "metric": metric,
                "source_file": str(path.relative_to(BASE)),
                "x_or_time": r[0],
                "value": r[1],
            })

    return summary, long_rows


def summarize_rmsf(path, phase, system_id, protein, state, temperature, source_root):
    rows, comments = read_xvg(path)
    if not rows:
        return None, []

    residue = [r[0] for r in rows if len(r) >= 2]
    values = [r[1] for r in rows if len(r) >= 2]

    if not values:
        return None, []

    mean_all, sd_all = mean_sd(values)
    max_idx = max(range(len(values)), key=lambda i: values[i])

    summary = {
        "phase": phase,
        "system_id": system_id,
        "protein": protein,
        "state": state,
        "temperature": temperature,
        "metric": "rmsf",
        "source_root": source_root,
        "source_file": str(path.relative_to(BASE)),
        "n_points": len(values),
        "time_start_ps_or_x": residue[0] if residue else None,
        "time_end_ps_or_x": residue[-1] if residue else None,
        "mean_all": mean_all,
        "sd_all": sd_all,
        "min_value": min(values),
        "max_value": max(values),
        "final_value": values[-1],
        "late50_mean": "",
        "late50_sd": "",
        "late50_start_index": "",
        "n_numeric_columns": max(len(r) for r in rows),
        "max_rmsf_residue_or_index": residue[max_idx],
        "max_rmsf_value": values[max_idx],
    }

    long_rows = []
    for r in rows:
        if len(r) >= 2:
            long_rows.append({
                "phase": phase,
                "system_id": system_id,
                "protein": protein,
                "state": state,
                "temperature": temperature,
                "metric": "rmsf",
                "source_file": str(path.relative_to(BASE)),
                "residue_or_index": r[0],
                "value": r[1],
            })

    return summary, long_rows


def summarize_ss_counts(path, phase, system_id, protein, state, temperature, source_root):
    rows, comments = read_xvg(path)
    if not rows:
        return None

    ncols = max(len(r) for r in rows)
    out = {
        "phase": phase,
        "system_id": system_id,
        "protein": protein,
        "state": state,
        "temperature": temperature,
        "metric": "secondary_structure_counts",
        "source_root": source_root,
        "source_file": str(path.relative_to(BASE)),
        "n_points": len(rows),
        "n_numeric_columns": ncols,
    }

    for col in range(1, ncols):
        vals = [r[col] for r in rows if len(r) > col]
        if vals:
            m, s = mean_sd(vals)
            out[f"col{col}_mean"] = m
            out[f"col{col}_sd"] = s
            out[f"col{col}_min"] = min(vals)
            out[f"col{col}_max"] = max(vals)
            out[f"col{col}_final"] = vals[-1]

    return out


def infer_temperature(path):
    parts = path.parts
    for p in parts:
        if re.fullmatch(r"\d{3}K", p):
            return p
    m = re.search(r"(\d{3}K)", str(path))
    return m.group(1) if m else ""


def infer_pdb_system(path):
    s = str(path)
    for pdb in PDB_IDS:
        if pdb in s:
            return pdb
    return ""


def infer_modeller_system(path):
    s = str(path)
    protein = ""
    for p in MOD_PROTEINS:
        if p in s:
            protein = p
            break

    state = ""
    if re.search(r"(^|[_/\-])WT([_/\-]|$)", s):
        state = "WT"
    elif re.search(r"(^|[_/\-])MUT([_/\-]|$)", s) or "mutant" in s.lower():
        state = "MUT"

    system_id = ""
    parts = list(path.parts)
    temp = infer_temperature(path)
    if temp and temp in parts:
        idx = parts.index(temp)
        if idx > 0:
            system_id = parts[idx - 1]
    if not system_id and protein:
        system_id = protein

    return system_id, protein, state


def classify_metric_file(path):
    name = path.name.lower()
    if "rmsf" in name:
        return "rmsf"
    if "rmsd" in name:
        return "rmsd"
    if "gyrate" in name or re.search(r"(^|_)rg(_|\.)", name) or name == "md_rg.xvg":
        return "rg"
    if "hbond" in name or "hbnum" in name:
        return "hbond"
    if "sasa" in name:
        return "sasa"
    if "ss_counts" in name:
        return "secondary_structure_counts"
    return ""


def collect_phase1_files():
    roots = [
        BASE / "md_10ns" / "systems",
        BASE / "md" / "systems",
        BASE / "md" / "_archive",
    ]

    files = []
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if is_gromacs_backup_file(path):
                continue
            if infer_pdb_system(path) == "":
                continue
            metric = classify_metric_file(path)
            if metric:
                files.append((path, root))
    return files


def collect_phase2_files():
    roots = [
        BASE / "md_tier2_wt_mutant_compact" / "systems",
        BASE / "md_tier2_wt_mutant_compact" / "analysis",
    ]

    files = []
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if is_gromacs_backup_file(path):
                continue
            metric = classify_metric_file(path)
            if metric:
                files.append((path, root))
    return files


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        path.write_text("")
        return

    fieldnames = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_csv_gz(path, rows):
    rows = list(rows)
    if not rows:
        with gzip.open(path, "wt", newline="") as handle:
            handle.write("")
        return

    fieldnames = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)

    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_wide(summary_rows):
    grouped = {}
    for r in summary_rows:
        key = (
            r.get("phase", ""),
            r.get("system_id", ""),
            r.get("protein", ""),
            r.get("state", ""),
            r.get("temperature", ""),
            r.get("source_root", ""),
        )
        if key not in grouped:
            grouped[key] = {
                "phase": key[0],
                "system_id": key[1],
                "protein": key[2],
                "state": key[3],
                "temperature": key[4],
                "source_root": key[5],
            }

        metric = r.get("metric", "")
        if metric in ["rmsd", "rg", "hbond", "sasa", "rmsf"]:
            for col in ["mean_all", "sd_all", "min_value", "max_value", "final_value", "late50_mean", "late50_sd", "n_points", "source_file"]:
                grouped[key][f"{metric}_{col}"] = r.get(col, "")

        if metric == "rmsf":
            grouped[key]["rmsf_max_residue_or_index"] = r.get("max_rmsf_residue_or_index", "")
            grouped[key]["rmsf_max_rmsf_value"] = r.get("max_rmsf_value", "")

    return list(grouped.values())


def phase1_extract():
    files = collect_phase1_files()
    manifest = []
    summaries = []
    timeseries_long = []
    rmsf_long = []
    ss_summaries = []

    for path, root in files:
        pdb = infer_pdb_system(path)
        temp = infer_temperature(path)
        metric = classify_metric_file(path)
        system_id = pdb

        manifest.append({
            "phase": "phase1_pdb_referenced_baseline",
            "system_id": system_id,
            "protein": pdb,
            "state": "PDB_WT",
            "temperature": temp,
            "metric": metric,
            "source_root": str(root.relative_to(BASE)),
            "source_file": str(path.relative_to(BASE)),
            "size_bytes": path.stat().st_size,
        })

        if metric == "rmsf":
            summary, long_rows = summarize_rmsf(
                path, "phase1_pdb_referenced_baseline", system_id, pdb, "PDB_WT", temp, str(root.relative_to(BASE))
            )
            if summary:
                summaries.append(summary)
                rmsf_long.extend(long_rows)

        elif metric == "secondary_structure_counts":
            ss = summarize_ss_counts(
                path, "phase1_pdb_referenced_baseline", system_id, pdb, "PDB_WT", temp, str(root.relative_to(BASE))
            )
            if ss:
                ss_summaries.append(ss)

        else:
            summary, long_rows = summarize_timeseries(
                path, "phase1_pdb_referenced_baseline", system_id, pdb, "PDB_WT", temp, metric, str(root.relative_to(BASE))
            )
            if summary:
                summaries.append(summary)
                timeseries_long.extend(long_rows)

    write_csv(PDB_OUT / "phase1_pdb_md_metric_summary_long.csv", summaries)
    write_csv(PDB_OUT / "phase1_pdb_md_metric_summary_wide.csv", make_wide(summaries))
    write_csv(PDB_OUT / "phase1_pdb_secondary_structure_counts_summary.csv", ss_summaries)
    write_csv(PDB_OUT / "phase1_pdb_original_metric_file_manifest.csv", manifest)
    write_csv_gz(PDB_OUT / "phase1_pdb_md_timeseries_long.csv.gz", timeseries_long)
    write_csv_gz(PDB_OUT / "phase1_pdb_rmsf_residue_long.csv.gz", rmsf_long)

    coverage_rows = []
    for pdb in PDB_IDS:
        for temp in PDB_TEMPS:
            present_metrics = sorted(set(
                r["metric"] for r in manifest
                if r["system_id"] == pdb and r["temperature"] == temp
            ))
            coverage_rows.append({
                "phase": "phase1_pdb_referenced_baseline",
                "system_id": pdb,
                "temperature": temp,
                "metrics_found": ";".join(present_metrics),
                "has_rmsd": "rmsd" in present_metrics,
                "has_rmsf": "rmsf" in present_metrics,
                "has_rg": "rg" in present_metrics,
                "has_hbond": "hbond" in present_metrics,
                "has_sasa": "sasa" in present_metrics,
                "has_secondary_structure_counts": "secondary_structure_counts" in present_metrics,
            })
    write_csv(PDB_OUT / "phase1_pdb_expected_coverage.csv", coverage_rows)

    return files, summaries, manifest


def phase2_extract():
    files = collect_phase2_files()
    manifest = []
    summaries = []
    timeseries_long = []
    rmsf_long = []
    ss_summaries = []

    for path, root in files:
        system_id, protein, state = infer_modeller_system(path)
        temp = infer_temperature(path)
        metric = classify_metric_file(path)

        if not protein:
            protein = system_id

        manifest.append({
            "phase": "phase2_modeller_wt_mutant_validation",
            "system_id": system_id,
            "protein": protein,
            "state": state,
            "temperature": temp,
            "metric": metric,
            "source_root": str(root.relative_to(BASE)),
            "source_file": str(path.relative_to(BASE)),
            "size_bytes": path.stat().st_size,
        })

        if metric == "rmsf":
            summary, long_rows = summarize_rmsf(
                path, "phase2_modeller_wt_mutant_validation", system_id, protein, state, temp, str(root.relative_to(BASE))
            )
            if summary:
                summaries.append(summary)
                rmsf_long.extend(long_rows)

        elif metric == "secondary_structure_counts":
            ss = summarize_ss_counts(
                path, "phase2_modeller_wt_mutant_validation", system_id, protein, state, temp, str(root.relative_to(BASE))
            )
            if ss:
                ss_summaries.append(ss)

        else:
            summary, long_rows = summarize_timeseries(
                path, "phase2_modeller_wt_mutant_validation", system_id, protein, state, temp, metric, str(root.relative_to(BASE))
            )
            if summary:
                summaries.append(summary)
                timeseries_long.extend(long_rows)

    write_csv(MOD_OUT / "phase2_modeller_md_metric_summary_long.csv", summaries)
    write_csv(MOD_OUT / "phase2_modeller_md_metric_summary_wide.csv", make_wide(summaries))
    write_csv(MOD_OUT / "phase2_modeller_secondary_structure_counts_summary.csv", ss_summaries)
    write_csv(MOD_OUT / "phase2_modeller_original_metric_file_manifest.csv", manifest)
    write_csv_gz(MOD_OUT / "phase2_modeller_md_timeseries_long.csv.gz", timeseries_long)
    write_csv_gz(MOD_OUT / "phase2_modeller_rmsf_residue_long.csv.gz", rmsf_long)

    coverage_rows = []
    for protein in MOD_PROTEINS:
        for state in ["WT", "MUT"]:
            for temp in MOD_TEMPS:
                present_metrics = sorted(set(
                    r["metric"] for r in manifest
                    if r["protein"] == protein and r["state"] == state and r["temperature"] == temp
                ))
                coverage_rows.append({
                    "phase": "phase2_modeller_wt_mutant_validation",
                    "protein": protein,
                    "state": state,
                    "temperature": temp,
                    "metrics_found": ";".join(present_metrics),
                    "has_rmsd": "rmsd" in present_metrics,
                    "has_rmsf": "rmsf" in present_metrics,
                    "has_rg": "rg" in present_metrics,
                    "has_hbond": "hbond" in present_metrics,
                    "has_sasa": "sasa" in present_metrics,
                    "has_secondary_structure_counts": "secondary_structure_counts" in present_metrics,
                })
    write_csv(MOD_OUT / "phase2_modeller_expected_coverage.csv", coverage_rows)

    return files, summaries, manifest


def main():
    pdb_files, pdb_summaries, pdb_manifest = phase1_extract()
    mod_files, mod_summaries, mod_manifest = phase2_extract()

    readme = f"""# MD metrics extracted from original analysis files

This directory contains molecular dynamics result tables extracted directly from original MD analysis files on the VM.

It does not use previously generated comparison CSV files as input.

## Phase 1: PDB-referenced baseline MD

Output folder:
phase1_pdb_referenced_baseline

Expected PDB systems:
{", ".join(PDB_IDS)}

Expected temperatures:
{", ".join(PDB_TEMPS)}

Main extracted files:
- phase1_pdb_md_metric_summary_long.csv
- phase1_pdb_md_metric_summary_wide.csv
- phase1_pdb_md_timeseries_long.csv.gz
- phase1_pdb_rmsf_residue_long.csv.gz
- phase1_pdb_secondary_structure_counts_summary.csv
- phase1_pdb_original_metric_file_manifest.csv
- phase1_pdb_expected_coverage.csv

## Phase 2: MODELLER / homology WT-mutant validation

Output folder:
phase2_modeller_wt_mutant_validation

Expected proteins:
{", ".join(MOD_PROTEINS)}

Expected temperatures:
{", ".join(MOD_TEMPS)}

Main extracted files:
- phase2_modeller_md_metric_summary_long.csv
- phase2_modeller_md_metric_summary_wide.csv
- phase2_modeller_md_timeseries_long.csv.gz
- phase2_modeller_rmsf_residue_long.csv.gz
- phase2_modeller_secondary_structure_counts_summary.csv
- phase2_modeller_original_metric_file_manifest.csv
- phase2_modeller_expected_coverage.csv

## Metrics

For time-series metrics, the script extracts:
- number of points
- start and end time/x value
- mean
- standard deviation
- minimum
- maximum
- final value
- late 50% mean
- late 50% standard deviation

For RMSF, the script extracts:
- mean RMSF
- standard deviation
- minimum RMSF
- maximum RMSF
- residue/index with maximum RMSF

## Extraction counts

PDB original metric files found: {len(pdb_files)}
PDB summary rows generated: {len(pdb_summaries)}
MODELLER original metric files found: {len(mod_files)}
MODELLER summary rows generated: {len(mod_summaries)}
"""
    (OUT / "README.md").write_text(readme)

    write_csv(MAN_OUT / "all_original_md_metric_files_manifest.csv", pdb_manifest + mod_manifest)

    print("Done.")
    print(f"Output: {OUT}")
    print(f"PDB original metric files found: {len(pdb_files)}")
    print(f"PDB summary rows generated: {len(pdb_summaries)}")
    print(f"MODELLER original metric files found: {len(mod_files)}")
    print(f"MODELLER summary rows generated: {len(mod_summaries)}")


if __name__ == "__main__":
    main()
