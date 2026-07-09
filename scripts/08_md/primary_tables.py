#!/usr/bin/env python3
# Purpose: Build primary MD tables.

from pathlib import Path
import csv

BASE = Path.home() / "xylanase-thesis"
ROOT = BASE / "results" / "md_metrics_extracted_from_original_files"

PHASES = {
    "phase1_pdb_referenced_baseline": ROOT / "phase1_pdb_referenced_baseline",
    "phase2_modeller_wt_mutant_validation": ROOT / "phase2_modeller_wt_mutant_validation",
}


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        path.write_text("")
        return

    fields = []
    seen = set()
    for row in rows:
        for k in row:
            if k not in seen:
                seen.add(k)
                fields.append(k)

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def float_or_zero(x):
    try:
        return float(x)
    except Exception:
        return 0.0


def primary_score(row):
    """
    Higher score = better primary result file.

    Rules:
    1. Prefer files inside /analysis/ when available.
    2. Prefer longer files with more points.
    3. Prefer clean final names over alternate copies such as xvg1.xvg.
    """
    sf = row.get("source_file", "")
    metric = row.get("metric", "")
    n_points = float_or_zero(row.get("n_points", 0))

    score = 0.0

    if "/analysis/" in sf:
        score += 100000

    score += n_points

    # Prefer clean standard analysis names.
    preferred_names = {
        "rmsd": ["rmsd_backbone.xvg", "rmsd_protein.xvg", "md_rmsd_backbone.xvg"],
        "rmsf": ["rmsf_backbone.xvg", "rmsf_protein.xvg", "md_rmsf_backbone.xvg"],
        "rg": ["rg_protein.xvg", "gyrate_protein.xvg", "md_rg.xvg"],
        "hbond": ["hbond_protein.xvg", "hbnum_protein.xvg"],
        "sasa": ["sasa_protein.xvg"],
    }

    for i, name in enumerate(preferred_names.get(metric, [])):
        if sf.endswith(name):
            score += 1000 - i

    # Penalise duplicated alternative SASA filename.
    if sf.endswith("sasa_protein.xvg1.xvg"):
        score -= 500

    return score


def make_primary(summary_rows):
    groups = {}

    for row in summary_rows:
        key = (
            row.get("phase", ""),
            row.get("system_id", ""),
            row.get("protein", ""),
            row.get("state", ""),
            row.get("temperature", ""),
            row.get("metric", ""),
        )
        groups.setdefault(key, []).append(row)

    primary = []
    duplicates = []

    for key, rows in groups.items():
        ranked = sorted(rows, key=primary_score, reverse=True)
        selected = dict(ranked[0])
        selected["primary_selection_score"] = primary_score(selected)
        selected["n_duplicate_candidates_for_same_metric"] = len(rows)
        primary.append(selected)

        for duplicate in ranked[1:]:
            dup = dict(duplicate)
            dup["primary_selected_source_file"] = selected.get("source_file", "")
            dup["primary_selection_score"] = primary_score(selected)
            dup["duplicate_selection_score"] = primary_score(dup)
            duplicates.append(dup)

    primary = sorted(
        primary,
        key=lambda r: (
            r.get("phase", ""),
            r.get("system_id", ""),
            r.get("protein", ""),
            r.get("state", ""),
            r.get("temperature", ""),
            r.get("metric", ""),
        ),
    )

    duplicates = sorted(
        duplicates,
        key=lambda r: (
            r.get("phase", ""),
            r.get("system_id", ""),
            r.get("protein", ""),
            r.get("state", ""),
            r.get("temperature", ""),
            r.get("metric", ""),
            r.get("source_file", ""),
        ),
    )

    return primary, duplicates


def make_wide(primary_rows):
    grouped = {}

    for row in primary_rows:
        key = (
            row.get("phase", ""),
            row.get("system_id", ""),
            row.get("protein", ""),
            row.get("state", ""),
            row.get("temperature", ""),
        )

        grouped.setdefault(key, {
            "phase": key[0],
            "system_id": key[1],
            "protein": key[2],
            "state": key[3],
            "temperature": key[4],
        })

        metric = row.get("metric", "")

        for col in [
            "n_points",
            "time_start_ps_or_x",
            "time_end_ps_or_x",
            "mean_all",
            "sd_all",
            "min_value",
            "max_value",
            "final_value",
            "late50_mean",
            "late50_sd",
            "source_file",
        ]:
            grouped[key][f"{metric}_{col}"] = row.get(col, "")

        if metric == "rmsf":
            grouped[key]["rmsf_max_residue_or_index"] = row.get("max_rmsf_residue_or_index", "")
            grouped[key]["rmsf_max_rmsf_value"] = row.get("max_rmsf_value", "")

    return list(grouped.values())


def main():
    for phase, folder in PHASES.items():
        if phase == "phase1_pdb_referenced_baseline":
            input_file = folder / "phase1_pdb_md_metric_summary_long.csv"
            primary_long = folder / "phase1_pdb_primary_md_metric_summary_long.csv"
            primary_wide = folder / "phase1_pdb_primary_md_metric_summary_wide.csv"
            duplicates_file = folder / "phase1_pdb_duplicate_metric_files_not_used_in_primary_table.csv"
        else:
            input_file = folder / "phase2_modeller_md_metric_summary_long.csv"
            primary_long = folder / "phase2_modeller_primary_md_metric_summary_long.csv"
            primary_wide = folder / "phase2_modeller_primary_md_metric_summary_wide.csv"
            duplicates_file = folder / "phase2_modeller_duplicate_metric_files_not_used_in_primary_table.csv"

        rows = read_csv(input_file)
        primary, duplicates = make_primary(rows)

        write_csv(primary_long, primary)
        write_csv(primary_wide, make_wide(primary))
        write_csv(duplicates_file, duplicates)

        print(f"{phase}")
        print(f"  Input summary rows: {len(rows)}")
        print(f"  Primary rows: {len(primary)}")
        print(f"  Duplicate/non-primary rows: {len(duplicates)}")
        print(f"  Wrote: {primary_long.relative_to(BASE)}")
        print(f"  Wrote: {primary_wide.relative_to(BASE)}")
        print(f"  Wrote: {duplicates_file.relative_to(BASE)}")


if __name__ == "__main__":
    main()
