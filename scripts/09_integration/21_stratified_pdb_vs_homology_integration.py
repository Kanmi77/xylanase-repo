#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import re
from datetime import datetime

import numpy as np
import pandas as pd


ROOT = Path(".").resolve()


def read_csv(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def find_col(df, regexes):
    for rgx in regexes:
        for c in df.columns:
            if re.search(rgx, c, flags=re.I):
                return c
    return None


def find_numeric_col(df, regexes):
    hits = []
    for c in df.columns:
        for rgx in regexes:
            if re.search(rgx, c, flags=re.I):
                vals = pd.to_numeric(df[c], errors="coerce")
                n = vals.notna().sum()
                if n > 0:
                    hits.append((c, n))
                break
    if not hits:
        return None
    hits = sorted(hits, key=lambda x: x[1], reverse=True)
    return hits[0][0]


def standardize_ids(df):
    df = df.copy()

    uniprot_col = find_col(df, [
        r"^uniprot$",
        r"uniprot.*accession",
        r"accession$",
        r"^protein$",
        r"protein_id",
    ])

    pdb_col = find_col(df, [
        r"^pdb_id$",
        r"pdb.*id",
        r"structure_id",
        r"^pdb$",
    ])

    mutation_col = find_col(df, [
        r"^mutation$",
        r"mutation_code",
        r"foldx_mutation_code",
        r"mutant",
        r"variant",
        r"substitution",
    ])

    candidate_col = find_col(df, [
        r"candidate_id",
        r"system_id",
        r"receptor",
        r"structure_name",
        r"mutant_model",
    ])

    if uniprot_col:
        df["std_uniprot"] = df[uniprot_col].astype(str).str.strip()
    else:
        df["std_uniprot"] = np.nan

    if pdb_col:
        df["std_pdb_id"] = df[pdb_col].astype(str).str.upper().str.strip()
    else:
        df["std_pdb_id"] = np.nan

    if mutation_col:
        df["std_mutation"] = df[mutation_col].astype(str).str.strip()
    else:
        df["std_mutation"] = np.nan

    if candidate_col:
        df["std_candidate_id"] = df[candidate_col].astype(str).str.strip()
    else:
        df["std_candidate_id"] = np.nan

    for c in ["std_uniprot", "std_pdb_id", "std_mutation", "std_candidate_id"]:
        df[c] = df[c].replace({"nan": np.nan, "None": np.nan, "": np.nan})

    return df


def robust_minmax(series, invert=False):
    s = pd.to_numeric(series, errors="coerce")

    if invert:
        s = -s

    valid = s.dropna()

    if valid.empty:
        return pd.Series(np.nan, index=s.index)

    lo = valid.quantile(0.05)
    hi = valid.quantile(0.95)

    if hi == lo:
        return pd.Series(0.5, index=s.index)

    clipped = s.clip(lo, hi)

    return (clipped - lo) / (hi - lo)


def discover_ml_score(df):
    ml_prob_col = find_numeric_col(df, [
        r"therm.*prob",
        r"prob.*therm",
        r"stable.*prob",
        r"prob.*stable",
        r"class.*prob",
        r"ml.*score",
        r"prediction.*score",
    ])

    ml_label_col = find_col(df, [
        r"pred.*thermal",
        r"thermal.*pred",
        r"predicted.*class",
        r"ml.*label",
        r"thermal_class",
    ])

    return ml_prob_col, ml_label_col


def merge_ml(branch_df, ml_paths, report):
    out = branch_df.copy()

    valid_ml_paths = []
    for p in ml_paths:
        path = Path(p)
        if not path.exists():
            continue
        if path.suffix.lower() != ".csv":
            continue
        valid_ml_paths.append(path)

    if not valid_ml_paths:
        report.append("- No valid ML CSV files were provided or discovered.")
        return out

    for i, path in enumerate(valid_ml_paths, start=1):
        try:
            ml = read_csv(path)
            ml = standardize_ids(ml)

            ml_prob_col, ml_label_col = discover_ml_score(ml)

            keep = ["std_uniprot", "std_pdb_id", "std_mutation"]
            rename = {}

            if ml_prob_col:
                keep.append(ml_prob_col)
                rename[ml_prob_col] = f"ml{i}_probability_score"

            if ml_label_col:
                keep.append(ml_label_col)
                rename[ml_label_col] = f"ml{i}_predicted_label"

            keep = [c for c in keep if c in ml.columns]
            ml_small = ml[keep].copy().rename(columns=rename)

            merge_key = None

            for key in ["std_uniprot", "std_pdb_id"]:
                if key in out.columns and key in ml_small.columns:
                    if out[key].notna().sum() > 0 and ml_small[key].notna().sum() > 0:
                        merge_key = key
                        break

            if merge_key is None:
                report.append(f"- ML file `{path}` could not be merged: no shared key.")
                continue

            # Collapse ML to one row per key.
            numeric_cols = ml_small.select_dtypes(include=[np.number]).columns.tolist()
            non_numeric_cols = [c for c in ml_small.columns if c not in numeric_cols and c != merge_key]

            agg = {}
            for c in numeric_cols:
                agg[c] = "mean"
            for c in non_numeric_cols:
                agg[c] = "first"

            ml_small = ml_small[ml_small[merge_key].notna()]
            ml_small = ml_small.groupby(merge_key, dropna=False).agg(agg).reset_index()

            before = len(out)
            out = out.merge(ml_small, on=merge_key, how="left")
            after = len(out)

            report.append(
                f"- ML file `{path}` merged by `{merge_key}`: rows {before} -> {after}."
            )

        except Exception as e:
            report.append(f"- Failed to merge ML file `{path}`: {e}")

    return out


def score_branch(df, branch_name):
    df = df.copy()

    foldx_col = find_numeric_col(df, [
        r"^foldx_ddg$",
        r"^ddg$",
        r"foldx.*ddg",
        r"delta.*g",
        r"mutant_ddg",
        r"foldx_energy_per_residue",
    ])

    docking_col = find_numeric_col(df, [
        r"delta_binding_mut_minus_wt_mean",
        r"delta_affinity_mut_minus_wt",
        r"delta.*binding.*mean",
        r"delta.*mut.*wt",
        r"mut.*minus.*wt",
        r"docking.*delta",
        r"vina.*delta",
    ])

    # Fallback if no delta exists: direct docking/binding score, lower/more negative is better.
    direct_docking_col = find_numeric_col(df, [
        r"mut_binding_energy_mean",
        r"mut_mean_top3_binding_mean",
        r"mutant_affinity",
        r"vina_affinity",
        r"docking_score",
        r"binding_energy",
        r"xylobiose",
        r"xylotriose",
    ])

    if foldx_col:
        df["score_foldx"] = robust_minmax(df[foldx_col], invert=True)
    else:
        df["score_foldx"] = np.nan

    if docking_col:
        df["score_docking"] = robust_minmax(df[docking_col], invert=True)
        docking_used = docking_col
    elif direct_docking_col:
        df["score_docking"] = robust_minmax(df[direct_docking_col], invert=True)
        docking_used = direct_docking_col
    else:
        df["score_docking"] = np.nan
        docking_used = None

    ml_prob_col = find_numeric_col(df, [
        r"ml\d+_probability_score",
        r"therm.*prob",
        r"prob.*therm",
        r"stable.*prob",
        r"prob.*stable",
        r"ml.*score",
    ])

    ml_label_col = find_col(df, [
        r"ml\d+_predicted_label",
        r"pred.*thermal",
        r"predicted.*class",
        r"thermal_class",
    ])

    if ml_prob_col:
        df["score_ml"] = robust_minmax(df[ml_prob_col], invert=False)
        ml_used = ml_prob_col
    elif ml_label_col:
        labels = df[ml_label_col].astype(str).str.lower()
        df["score_ml"] = np.where(
            labels.str.contains("thermo|stable|high"),
            1.0,
            np.where(labels.str.contains("meso|lower|low"), 0.0, np.nan)
        )
        ml_used = ml_label_col
    else:
        df["score_ml"] = np.nan
        ml_used = None

    weights = {
        "score_foldx": 0.40,
        "score_docking": 0.35,
        "score_ml": 0.25,
    }

    weighted = pd.Series(0.0, index=df.index)
    used_weight = pd.Series(0.0, index=df.index)

    for col, w in weights.items():
        valid = df[col].notna()
        weighted.loc[valid] += df.loc[valid, col] * w
        used_weight.loc[valid] += w

    df["branch_integrated_score"] = np.where(
        used_weight > 0,
        weighted / used_weight,
        np.nan
    )

    df["evidence_count"] = df[["score_foldx", "score_docking", "score_ml"]].notna().sum(axis=1)

    df["source_branch"] = branch_name

    df["branch_rank"] = df["branch_integrated_score"].rank(
        method="dense",
        ascending=False
    )

    df = df.sort_values(
        ["branch_integrated_score", "evidence_count"],
        ascending=[False, False]
    )

    used = {
        "foldx_column_used": foldx_col,
        "docking_column_used": docking_used,
        "ml_column_used": ml_used,
        "weights": weights,
    }

    return df, used


def assign_evidence_tier(row):
    branch = row.get("source_branch", "")

    has_foldx = pd.notna(row.get("score_foldx"))
    has_docking = pd.notna(row.get("score_docking"))
    has_ml = pd.notna(row.get("score_ml"))

    if branch == "PDB":
        if has_foldx and has_docking:
            return "Tier A - PDB-supported primary"
        return "Tier A/B - PDB-supported incomplete"

    if branch == "Homology":
        if has_foldx and has_docking and has_ml:
            return "Tier B - model-supported multi-method"
        if has_foldx and has_docking:
            return "Tier B - model-supported structural"
        return "Tier C - model-supported exploratory"

    return "Unclassified"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-file", required=True)
    parser.add_argument("--homology-file", required=True)
    parser.add_argument("--ml", nargs="*", default=[])
    parser.add_argument(
        "--outdir",
        default="results/integration/stratified_pdb_vs_homology_new_ml_pdb_foldx_docking"
    )
    args = parser.parse_args()

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    report = []
    sources = {}

    report.append("# Stratified PDB vs homology integration report\n")
    report.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n")

    pdb_path = ROOT / args.pdb_file
    homology_path = ROOT / args.homology_file

    pdb_df = read_csv(pdb_path)
    homology_df = read_csv(homology_path)

    pdb_df = standardize_ids(pdb_df)
    homology_df = standardize_ids(homology_df)

    report.append("## Source files\n")
    report.append(f"- PDB branch file: `{args.pdb_file}`")
    report.append(f"- Homology branch file: `{args.homology_file}`")
    report.append(f"- PDB input rows: {len(pdb_df)}")
    report.append(f"- Homology input rows: {len(homology_df)}")

    ml_paths = []

    for item in args.ml:
        p = Path(item)
        if p.exists() and p.suffix.lower() == ".csv":
            ml_paths.append(str(p))

    report.append(f"- ML CSV files supplied/found: {len(ml_paths)}")
    for p in ml_paths:
        report.append(f"  - `{p}`")

    report.append("\n## ML merge log\n")
    pdb_df = merge_ml(pdb_df, ml_paths, report)
    homology_df = merge_ml(homology_df, ml_paths, report)

    pdb_ranked, pdb_used = score_branch(pdb_df, "PDB")
    homology_ranked, homology_used = score_branch(homology_df, "Homology")

    combined = pd.concat([pdb_ranked, homology_ranked], ignore_index=True, sort=False)
    combined["evidence_tier"] = combined.apply(assign_evidence_tier, axis=1)

    tier_order = {
        "Tier A - PDB-supported primary": 1,
        "Tier A/B - PDB-supported incomplete": 2,
        "Tier B - model-supported multi-method": 3,
        "Tier B - model-supported structural": 4,
        "Tier C - model-supported exploratory": 5,
        "Unclassified": 9,
    }

    combined["evidence_tier_order"] = combined["evidence_tier"].map(tier_order).fillna(9)

    combined = combined.sort_values(
        ["evidence_tier_order", "branch_rank", "branch_integrated_score"],
        ascending=[True, True, False]
    )

    pdb_out = outdir / "pdb_branch_integrated.csv"
    homology_out = outdir / "homology_branch_integrated.csv"
    combined_out = outdir / "final_evidence_tiered_summary.csv"
    top_pdb_out = outdir / "top15_pdb_branch.csv"
    top_homology_out = outdir / "top15_homology_branch.csv"
    report_out = outdir / "STRATIFIED_INTEGRATION_REPORT.md"
    sources_out = outdir / "integration_source_files.json"

    pdb_ranked.to_csv(pdb_out, index=False)
    homology_ranked.to_csv(homology_out, index=False)
    combined.to_csv(combined_out, index=False)
    pdb_ranked.head(15).to_csv(top_pdb_out, index=False)
    homology_ranked.head(15).to_csv(top_homology_out, index=False)

    sources = {
        "pdb_file": args.pdb_file,
        "homology_file": args.homology_file,
        "ml_files": ml_paths,
        "pdb_score_columns_used": pdb_used,
        "homology_score_columns_used": homology_used,
    }

    with open(sources_out, "w") as f:
        json.dump(sources, f, indent=2)

    report.append("\n## Score columns used\n")
    report.append("\n### PDB branch\n")
    for k, v in pdb_used.items():
        report.append(f"- `{k}`: `{v}`")

    report.append("\n### Homology branch\n")
    for k, v in homology_used.items():
        report.append(f"- `{k}`: `{v}`")

    report.append("\n## Output files\n")
    for p in [pdb_out, homology_out, combined_out, top_pdb_out, top_homology_out, report_out, sources_out]:
        report.append(f"- `{p.relative_to(ROOT)}`")

    preview_cols = [
        c for c in [
            "source_branch",
            "evidence_tier",
            "protein",
            "uniprot_accession",
            "mutation",
            "candidate_id",
            "std_uniprot",
            "std_pdb_id",
            "std_mutation",
            "branch_rank",
            "branch_integrated_score",
            "score_foldx",
            "score_docking",
            "score_ml",
            "evidence_count",
        ]
        if c in combined.columns
    ]

    report.append("\n## Final evidence-tiered preview\n")
    report.append(combined[preview_cols].head(40).to_markdown(index=False))

    with open(report_out, "w") as f:
        f.write("\n".join(report))
        f.write("\n")

    print(f"Saved: {pdb_out}")
    print(f"Saved: {homology_out}")
    print(f"Saved: {combined_out}")
    print(f"Saved: {top_pdb_out}")
    print(f"Saved: {top_homology_out}")
    print(f"Saved: {report_out}")
    print(f"Saved: {sources_out}")

    print("\nFinal evidence-tiered preview:")
    print(combined[preview_cols].head(40).to_string(index=False))


if __name__ == "__main__":
    main()
