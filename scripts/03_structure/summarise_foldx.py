#!/usr/bin/env python3
# Purpose: Summarise FoldX structural features.

from pathlib import Path
import pandas as pd
import numpy as np
import json

OUTDIR = Path("results/structural_features_foldx_summary")
OUTDIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_FILES = [
    "results/structures/structural_features.csv",
    "results/structural_features.csv",
    "data/curated/structural_features.csv",
    "results/foldx/foldx_normalized.csv",
    "results/foldx/foldx_with_structure_full.csv",
    "results/foldx/foldx_structural_ranked.csv",
    "results/foldx/modeller_foldx_stability.csv",
    "results/foldx4_wt_stability.csv",
]

existing = [Path(p) for p in CANDIDATE_FILES if Path(p).exists()]

if not existing:
    raise SystemExit("No candidate files found. Run find command and adjust CANDIDATE_FILES in the script.")

def read_any(path):
    if path.suffix.lower() in [".tsv", ".txt"]:
        return pd.read_csv(path, sep=None, engine="python")
    return pd.read_csv(path)

def find_col(df, options):
    lower = {c.lower(): c for c in df.columns}
    for opt in options:
        if opt.lower() in lower:
            return lower[opt.lower()]
    for c in df.columns:
        cl = c.lower()
        for opt in options:
            if opt.lower() in cl:
                return c
    return None

loaded = {}
for p in existing:
    try:
        df = read_any(p)
        loaded[str(p)] = df
        print(f"Loaded {p}: {df.shape}")
    except Exception as e:
        print(f"Skipped {p}: {e}")

# Choose likely structural feature file: one containing hbond/salt/sasa columns
feature_df = None
feature_path = None
for p, df in loaded.items():
    cols = " ".join(df.columns).lower()
    if ("sasa" in cols and ("hbond" in cols or "h_bond" in cols or "hydrogen" in cols or "salt" in cols)):
        feature_df = df.copy()
        feature_path = p
        break

# Choose FoldX file: one containing foldx/energy columns
foldx_df = None
foldx_path = None
for p, df in loaded.items():
    cols = " ".join(df.columns).lower()
    if "foldx" in cols or "energy" in cols or "stability" in cols:
        foldx_df = df.copy()
        foldx_path = p
        break

report_lines = []
report_lines.append("# Structural features and FoldX summary report\n")

if feature_df is not None:
    df = feature_df
    report_lines.append(f"## Structural feature file\n")
    report_lines.append(f"- Source file: `{feature_path}`")
    report_lines.append(f"- Rows: {len(df)}")
    report_lines.append(f"- Columns: {list(df.columns)}\n")

    accession_col = find_col(df, ["uniprot_accession", "accession", "protein", "query_accession"])
    organism_col = find_col(df, ["organism_type", "query_organism_type", "organism"])
    gh_col = find_col(df, ["gh_family", "query_gh_family", "family"])

    length_col = find_col(df, ["chain_length", "length", "sequence_length", "n_residues", "residue_count"])
    hbond_col = find_col(df, ["hbond_count", "h_bond_count", "hydrogen_bond_count", "hbond_proxy_count", "hydrogen_bonds"])
    salt_col = find_col(df, ["salt_bridge_count", "salt_bridges", "saltbridge_count"])
    disulfide_col = find_col(df, ["disulfide_count", "disulfide_bond_count", "disulfides"])
    sasa_col = find_col(df, ["total_sasa", "sasa", "sasa_total"])

    report_lines.append("### Detected structural columns\n")
    for name, col in [
        ("accession", accession_col),
        ("organism_type", organism_col),
        ("gh_family", gh_col),
        ("chain_length", length_col),
        ("hbond_count", hbond_col),
        ("salt_bridge_count", salt_col),
        ("disulfide_count", disulfide_col),
        ("total_sasa", sasa_col),
    ]:
        report_lines.append(f"- {name}: `{col}`")
    report_lines.append("")

    summary_rows = []
    for label, col in [
        ("Chain length", length_col),
        ("Hydrogen-bond proxy count", hbond_col),
        ("Salt-bridge count", salt_col),
        ("Disulfide-bond count", disulfide_col),
        ("Total SASA", sasa_col),
    ]:
        if col and col in df.columns:
            x = pd.to_numeric(df[col], errors="coerce").dropna()
            summary_rows.append({
                "feature": label,
                "valid_records": int(x.shape[0]),
                "minimum": x.min(),
                "mean": x.mean(),
                "median": x.median(),
                "maximum": x.max(),
            })

    raw_summary = pd.DataFrame(summary_rows)
    raw_summary.to_csv(OUTDIR / "structural_feature_raw_summary.csv", index=False)

    report_lines.append("### Raw structural feature summary\n")
    report_lines.append(raw_summary.to_markdown(index=False))
    report_lines.append("")

    # Normalised columns
    if length_col:
        length = pd.to_numeric(df[length_col], errors="coerce")
        if hbond_col:
            df["hbond_per_residue_calc"] = pd.to_numeric(df[hbond_col], errors="coerce") / length
        if salt_col:
            df["salt_bridge_per_residue_calc"] = pd.to_numeric(df[salt_col], errors="coerce") / length
        if disulfide_col:
            df["disulfide_per_residue_calc"] = pd.to_numeric(df[disulfide_col], errors="coerce") / length
        if sasa_col:
            df["sasa_per_residue_calc"] = pd.to_numeric(df[sasa_col], errors="coerce") / length

    norm_cols = [
        c for c in ["hbond_per_residue_calc", "salt_bridge_per_residue_calc", "disulfide_per_residue_calc", "sasa_per_residue_calc"]
        if c in df.columns
    ]

    norm_summary_rows = []
    for col in norm_cols:
        x = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        norm_summary_rows.append({
            "feature": col,
            "valid_records": int(x.shape[0]),
            "minimum": x.min(),
            "mean": x.mean(),
            "median": x.median(),
            "maximum": x.max(),
        })

    norm_summary = pd.DataFrame(norm_summary_rows)
    norm_summary.to_csv(OUTDIR / "structural_feature_normalised_summary.csv", index=False)

    report_lines.append("### Normalised structural feature summary\n")
    report_lines.append(norm_summary.to_markdown(index=False))
    report_lines.append("")

    if organism_col and gh_col and norm_cols:
        group = df.copy()
        group["group"] = group[organism_col].astype(str) + " " + group[gh_col].astype(str)
        group_summary = (
            group.groupby("group")
            .agg(
                records=("group", "size"),
                mean_chain_length=(length_col, "mean") if length_col else ("group", "size"),
                mean_hbond_per_residue=("hbond_per_residue_calc", "mean") if "hbond_per_residue_calc" in group.columns else ("group", "size"),
                mean_salt_bridge_per_residue=("salt_bridge_per_residue_calc", "mean") if "salt_bridge_per_residue_calc" in group.columns else ("group", "size"),
                mean_disulfide_per_residue=("disulfide_per_residue_calc", "mean") if "disulfide_per_residue_calc" in group.columns else ("group", "size"),
                mean_sasa_per_residue=("sasa_per_residue_calc", "mean") if "sasa_per_residue_calc" in group.columns else ("group", "size"),
            )
            .reset_index()
        )
        group_summary.to_csv(OUTDIR / "normalised_feature_profiles_by_group.csv", index=False)

        report_lines.append("### Normalised feature profiles by group\n")
        report_lines.append(group_summary.to_markdown(index=False))
        report_lines.append("")

    df.to_csv(OUTDIR / "structural_features_with_calculated_normalised_values.csv", index=False)

else:
    report_lines.append("## Structural feature file\n")
    report_lines.append("- No structural feature file detected.\n")

if foldx_df is not None:
    df = foldx_df
    report_lines.append(f"## FoldX file\n")
    report_lines.append(f"- Source file: `{foldx_path}`")
    report_lines.append(f"- Rows: {len(df)}")
    report_lines.append(f"- Columns: {list(df.columns)}\n")

    organism_col = find_col(df, ["organism_type", "query_organism_type", "organism"])
    gh_col = find_col(df, ["gh_family", "query_gh_family", "family"])
    energy_col = find_col(df, ["foldx_energy", "total_energy", "stability", "foldx_total_energy", "energy"])
    energy_res_col = find_col(df, ["foldx_energy_per_residue", "energy_per_residue", "foldx_per_residue", "stability_per_residue"])

    report_lines.append("### Detected FoldX columns\n")
    for name, col in [
        ("organism_type", organism_col),
        ("gh_family", gh_col),
        ("foldx_energy", energy_col),
        ("foldx_energy_per_residue", energy_res_col),
    ]:
        report_lines.append(f"- {name}: `{col}`")
    report_lines.append("")

    foldx_rows = []
    for label, col in [
        ("FoldX total energy", energy_col),
        ("FoldX energy per residue", energy_res_col),
    ]:
        if col and col in df.columns:
            x = pd.to_numeric(df[col], errors="coerce").dropna()
            foldx_rows.append({
                "feature": label,
                "valid_records": int(x.shape[0]),
                "minimum": x.min(),
                "mean": x.mean(),
                "median": x.median(),
                "maximum": x.max(),
            })

    foldx_summary = pd.DataFrame(foldx_rows)
    foldx_summary.to_csv(OUTDIR / "foldx_wild_type_summary.csv", index=False)

    report_lines.append("### FoldX wild-type summary\n")
    report_lines.append(foldx_summary.to_markdown(index=False))
    report_lines.append("")

    if organism_col and gh_col and energy_res_col:
        group = df.copy()
        group["group"] = group[organism_col].astype(str) + " " + group[gh_col].astype(str)
        group[energy_res_col] = pd.to_numeric(group[energy_res_col], errors="coerce")
        group_summary = (
            group.groupby("group")
            .agg(
                records=("group", "size"),
                mean_foldx_energy_per_residue=(energy_res_col, "mean"),
                median_foldx_energy_per_residue=(energy_res_col, "median"),
                min_foldx_energy_per_residue=(energy_res_col, "min"),
                max_foldx_energy_per_residue=(energy_res_col, "max"),
            )
            .reset_index()
        )
        group_summary.to_csv(OUTDIR / "foldx_energy_per_residue_by_group.csv", index=False)

        report_lines.append("### FoldX energy per residue by group\n")
        report_lines.append(group_summary.to_markdown(index=False))
        report_lines.append("")

else:
    report_lines.append("## FoldX file\n")
    report_lines.append("- No FoldX file detected.\n")

with open(OUTDIR / "STRUCTURAL_FEATURES_FOLDX_SUMMARY_REPORT.md", "w") as h:
    h.write("\n".join(report_lines))

print("Done.")
print("Report:", OUTDIR / "STRUCTURAL_FEATURES_FOLDX_SUMMARY_REPORT.md")
