#!/usr/bin/env python3
# Purpose: Verify PDB ligand docking.

from pathlib import Path
import pandas as pd

OUT = Path("results/verified_docking_results")
OUT.mkdir(parents=True, exist_ok=True)

print("Searching for full PDB paired docking delta wide candidate table...\n")

full_candidates = []
for p in Path("results").rglob("*.csv"):
    s = str(p).lower()
    if "pdb60_paired_docking_delta_wide_candidate_table" in s and "top30" not in s:
        full_candidates.append(p)

if not full_candidates:
    raise SystemExit("Full PDB paired wide table not found.")

infile = full_candidates[0]
df = pd.read_csv(infile)

print("INPUT USED:", infile)
print("shape:", df.shape)
print("columns:", list(df.columns))
print(df.head(10).to_string(index=False))

# Accept either mean column name
if "mean_docking_change_calc" not in df.columns and "mean_delta_affinity" in df.columns:
    df["mean_docking_change_calc"] = df["mean_delta_affinity"]

required = [
    "uniprot_accession",
    "organism_type",
    "gh_family",
    "foldx_mutation_code",
    "foldx_ddg",
    "delta_affinity_mut_minus_wt_xylobiose",
    "delta_affinity_mut_minus_wt_xylotriose",
    "mean_docking_change_calc",
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing required columns: {missing}")

df["foldx_mutation_code"] = df["foldx_mutation_code"].astype(str).str.replace(";", "", regex=False)

for c in [
    "foldx_ddg",
    "delta_affinity_mut_minus_wt_xylobiose",
    "delta_affinity_mut_minus_wt_xylotriose",
    "mean_docking_change_calc",
]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Convert wide table to ligand-long table
long = pd.concat([
    df.assign(
        ligand="xylobiose",
        delta_binding_mut_minus_wt=df["delta_affinity_mut_minus_wt_xylobiose"]
    ),
    df.assign(
        ligand="xylotriose",
        delta_binding_mut_minus_wt=df["delta_affinity_mut_minus_wt_xylotriose"]
    )
], ignore_index=True)

# Ligand/group summaries
group = (
    long.groupby(["ligand", "organism_type", "gh_family"], dropna=False)
        .agg(
            docking_records=("uniprot_accession", "count"),
            mean_delta=("delta_binding_mut_minus_wt", "mean"),
            median_delta=("delta_binding_mut_minus_wt", "median"),
            min_delta=("delta_binding_mut_minus_wt", "min"),
            max_delta=("delta_binding_mut_minus_wt", "max"),
            improved_or_retained_count=("delta_binding_mut_minus_wt", lambda x: (x <= 0).sum())
        )
        .reset_index()
)

group["improved_or_retained_percent"] = (
    group["improved_or_retained_count"] / group["docking_records"] * 100
)

for c in group.columns:
    if pd.api.types.is_float_dtype(group[c]):
        group[c] = group[c].round(4)

group.to_csv(OUT / "verified_pdb_ligand_group_summary.csv", index=False)

candidate_cols = [
    "uniprot_accession",
    "organism_type",
    "gh_family",
    "foldx_mutation_code",
    "foldx_ddg",
    "ligand",
    "delta_binding_mut_minus_wt",
]

cand = long[candidate_cols].rename(columns={
    "uniprot_accession": "candidate",
    "foldx_mutation_code": "mutation",
})

for ligand in ["xylobiose", "xylotriose"]:
    sub = cand[cand["ligand"] == ligand].copy()

    best = sub.sort_values("delta_binding_mut_minus_wt", ascending=True).head(15)
    worst = sub.sort_values("delta_binding_mut_minus_wt", ascending=False).head(15)

    for table in [best, worst]:
        for c in table.columns:
            if pd.api.types.is_float_dtype(table[c]):
                table[c] = table[c].round(4)

    best.to_csv(OUT / f"verified_pdb_{ligand}_top15_favourable.csv", index=False)
    worst.to_csv(OUT / f"verified_pdb_{ligand}_top15_weaker.csv", index=False)

overall = df[[
    "uniprot_accession",
    "organism_type",
    "gh_family",
    "foldx_mutation_code",
    "foldx_ddg",
    "delta_affinity_mut_minus_wt_xylobiose",
    "delta_affinity_mut_minus_wt_xylotriose",
    "mean_docking_change_calc",
]].copy()

overall = overall.rename(columns={
    "uniprot_accession": "candidate",
    "foldx_mutation_code": "mutation",
    "delta_affinity_mut_minus_wt_xylobiose": "xylobiose_docking_change",
    "delta_affinity_mut_minus_wt_xylotriose": "xylotriose_docking_change",
    "mean_docking_change_calc": "mean_docking_change",
})

overall["retention_class"] = overall.apply(
    lambda r: "improved_both_ligands" if r["xylobiose_docking_change"] <= 0 and r["xylotriose_docking_change"] <= 0
    else ("mixed_retention" if r["xylobiose_docking_change"] <= 0 or r["xylotriose_docking_change"] <= 0
          else "docking_loss_both_ligands"),
    axis=1
)

for c in overall.columns:
    if pd.api.types.is_float_dtype(overall[c]):
        overall[c] = overall[c].round(4)

overall_sorted = overall.sort_values("mean_docking_change", ascending=True)
overall_sorted.to_csv(
    OUT / "verified_pdb_overall_candidates_by_mean_docking_change.csv",
    index=False
)

retention = (
    overall.groupby(["retention_class"], dropna=False)
           .agg(candidate_count=("candidate", "count"))
           .reset_index()
)
retention["percent"] = retention["candidate_count"] / retention["candidate_count"].sum() * 100
retention["percent"] = retention["percent"].round(1)
retention.to_csv(OUT / "verified_pdb_retention_classification.csv", index=False)

print("\n=== PDB LIGAND GROUP SUMMARY ===")
print(group.to_string(index=False))

for ligand in ["xylobiose", "xylotriose"]:
    print(f"\n=== PDB {ligand.upper()} TOP 15 FAVOURABLE ===")
    print(pd.read_csv(OUT / f"verified_pdb_{ligand}_top15_favourable.csv").to_string(index=False))

    print(f"\n=== PDB {ligand.upper()} TOP 15 WEAKER ===")
    print(pd.read_csv(OUT / f"verified_pdb_{ligand}_top15_weaker.csv").to_string(index=False))

print("\n=== PDB RETENTION CLASSIFICATION ===")
print(retention.to_string(index=False))

print("\n=== PDB OVERALL TOP 15 FAVOURABLE BY MEAN DOCKING CHANGE ===")
print(overall_sorted.head(15).to_string(index=False))

print("\n=== PDB OVERALL TOP 15 WEAKER BY MEAN DOCKING CHANGE ===")
print(overall.sort_values("mean_docking_change", ascending=False).head(15).to_string(index=False))

print("\nSaved files in:", OUT)
