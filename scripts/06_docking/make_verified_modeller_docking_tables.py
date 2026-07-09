#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

OUT = Path("results/verified_docking_results")
OUT.mkdir(parents=True, exist_ok=True)

# Confirmed full MODELLER/homology file
infile = Path("results/integration/stratified_pdb_vs_homology_new_ml_pdb_foldx_docking/homology_branch_integrated_annotated_metadata_filled.csv")

if not infile.exists():
    infile = Path("results/integration/homology_full_docking_candidate_level.csv")

df = pd.read_csv(infile)

print("INPUT FILE:", infile)
print("shape:", df.shape)
print("columns:", list(df.columns))

# Standardise organism/GH columns
org_col = None
gh_col = None

for c in ["candidate_organism_type", "candidate_organism_type_source", "organism_type"]:
    if c in df.columns:
        org_col = c
        break

for c in ["candidate_gh_family", "candidate_gh_family_source", "gh_family"]:
    if c in df.columns:
        gh_col = c
        break

if org_col is None or gh_col is None:
    raise SystemExit("Could not find organism type or GH family columns.")

# Standardise key docking columns
needed = [
    "delta_binding_mut_minus_wt_mean",
    "delta_binding_mut_minus_wt_best",
    "delta_binding_mut_minus_wt_worst",
    "wt_binding_energy_mean",
    "mut_binding_energy_mean",
    "foldx_ddg"
]

missing = [c for c in needed if c not in df.columns]
if missing:
    raise SystemExit(f"Missing required columns: {missing}")

for c in needed:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Full group summary using mean paired mutant-WT docking change
group = (
    df.groupby([org_col, gh_col], dropna=False)
      .agg(
          candidate_count=("uniprot_accession", "count"),
          mean_foldx_ddg=("foldx_ddg", "mean"),
          mean_docking_change=("delta_binding_mut_minus_wt_mean", "mean"),
          median_docking_change=("delta_binding_mut_minus_wt_mean", "median"),
          min_docking_change=("delta_binding_mut_minus_wt_mean", "min"),
          max_docking_change=("delta_binding_mut_minus_wt_mean", "max"),
          mean_wt_binding=("wt_binding_energy_mean", "mean"),
          mean_mut_binding=("mut_binding_energy_mean", "mean"),
          improved_or_retained_count=("mean_docking_improved_vs_wt", "sum") if "mean_docking_improved_vs_wt" in df.columns else ("delta_binding_mut_minus_wt_mean", lambda x: (x <= 0).sum())
      )
      .reset_index()
)

group["improved_or_retained_percent"] = (
    group["improved_or_retained_count"] / group["candidate_count"] * 100
)

group = group.rename(columns={
    org_col: "organism_type",
    gh_col: "gh_family"
})

# Round only for output readability
group_out = group.copy()
for c in group_out.columns:
    if pd.api.types.is_float_dtype(group_out[c]):
        group_out[c] = group_out[c].round(4)

group_out.to_csv(OUT / "verified_modeller_docking_group_summary.csv", index=False)

# Full top/worst candidate tables from confirmed 150-row file
id_col = "uniprot_accession" if "uniprot_accession" in df.columns else "protein"
mut_col = "mutation"

candidate_cols = [
    id_col, org_col, gh_col, mut_col, "foldx_ddg",
    "delta_binding_mut_minus_wt_mean",
    "delta_binding_mut_minus_wt_best",
    "delta_binding_mut_minus_wt_worst",
    "wt_binding_energy_mean",
    "mut_binding_energy_mean"
]

candidate_cols = [c for c in candidate_cols if c in df.columns]

best = df.sort_values("delta_binding_mut_minus_wt_mean", ascending=True)[candidate_cols].head(20)
worst = df.sort_values("delta_binding_mut_minus_wt_mean", ascending=False)[candidate_cols].head(20)

best = best.rename(columns={
    id_col: "candidate",
    org_col: "organism_type",
    gh_col: "gh_family",
    "delta_binding_mut_minus_wt_mean": "mean_docking_change",
    "delta_binding_mut_minus_wt_best": "best_ligand_change",
    "delta_binding_mut_minus_wt_worst": "worst_ligand_change"
})

worst = worst.rename(columns={
    id_col: "candidate",
    org_col: "organism_type",
    gh_col: "gh_family",
    "delta_binding_mut_minus_wt_mean": "mean_docking_change",
    "delta_binding_mut_minus_wt_best": "best_ligand_change",
    "delta_binding_mut_minus_wt_worst": "worst_ligand_change"
})

for table in [best, worst]:
    for c in table.columns:
        if pd.api.types.is_float_dtype(table[c]):
            table[c] = table[c].round(4)

best.to_csv(OUT / "verified_modeller_top20_favourable_candidates.csv", index=False)
worst.to_csv(OUT / "verified_modeller_top20_weaker_candidates.csv", index=False)

# Write report
report = []
report.append("# Verified MODELLER docking result tables\n")
report.append(f"Input file: `{infile}`")
report.append(f"Input shape: {df.shape}\n")

report.append("## MODELLER group-level docking summary\n")
report.append(group_out.to_string(index=False))

report.append("\n\n## Top 20 favourable MODELLER candidates\n")
report.append(best.to_string(index=False))

report.append("\n\n## Top 20 weaker MODELLER candidates\n")
report.append(worst.to_string(index=False))

(OUT / "VERIFIED_MODELLER_DOCKING_RESULTS_REPORT.md").write_text("\n".join(report))

print("\n=== MODELLER GROUP SUMMARY ===")
print(group_out.to_string(index=False))

print("\n=== TOP 20 FAVOURABLE MODELLER CANDIDATES ===")
print(best.to_string(index=False))

print("\n=== TOP 20 WEAKER MODELLER CANDIDATES ===")
print(worst.to_string(index=False))

print("\nSaved:")
for p in sorted(OUT.glob("*modeller*")):
    print(" -", p)
