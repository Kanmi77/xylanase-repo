#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

OUT = Path("results/verified_docking_results")
OUT.mkdir(parents=True, exist_ok=True)

# Find ligand-level MODELLER docking rows
candidate_files = []
for p in Path("results").rglob("*.csv"):
    s = str(p).lower()
    if any(k in s for k in ["docking", "vina", "ligand", "xylobiose", "xylotriose"]):
        candidate_files.append(p)

print("Searching for MODELLER ligand-level docking files...\n")

valid = []
for p in sorted(candidate_files):
    try:
        df = pd.read_csv(p)
    except Exception:
        continue

    cols = list(df.columns)
    lowcols = [c.lower() for c in cols]

    has_ligand = any(c in lowcols for c in ["ligand", "ligand_name"])
    has_delta = any("delta" in c and ("binding" in c or "affinity" in c) for c in lowcols)
    has_model = not any(k in str(p).lower() for k in ["pdb60", "pdb_branch", "top_pdb60"])

    if has_ligand and has_delta and has_model:
        valid.append(p)
        print("POSSIBLE MODELLER LIGAND FILE:", p)
        print("shape:", df.shape)
        print("columns:", cols)
        print(df.head(5).to_string(index=False))
        print()

if not valid:
    raise SystemExit("No ligand-level MODELLER docking file found. Paste the search output above.")

# Prefer the largest valid file
best_file = None
best_rows = -1
for p in valid:
    df = pd.read_csv(p)
    if df.shape[0] > best_rows:
        best_file = p
        best_rows = df.shape[0]

df = pd.read_csv(best_file)
print("USING:", best_file)
print("shape:", df.shape)

# Standardise columns
cols = list(df.columns)
lower_map = {c.lower(): c for c in cols}

lig_col = lower_map.get("ligand", lower_map.get("ligand_name"))

delta_col = None
for c in cols:
    lc = c.lower()
    if "delta" in lc and ("binding" in lc or "affinity" in lc):
        delta_col = c
        break

protein_col = None
for c in ["protein", "uniprot_accession", "candidate", "candidate_id"]:
    if c in cols:
        protein_col = c
        break

mut_col = None
for c in ["mutation", "foldx_mutation_code"]:
    if c in cols:
        mut_col = c
        break

org_col = None
for c in ["organism_type", "candidate_organism_type", "candidate_organism_type_source"]:
    if c in cols:
        org_col = c
        break

gh_col = None
for c in ["gh_family", "candidate_gh_family", "candidate_gh_family_source"]:
    if c in cols:
        gh_col = c
        break

if not all([lig_col, delta_col, protein_col, mut_col]):
    raise SystemExit(f"Missing required columns. ligand={lig_col}, delta={delta_col}, protein={protein_col}, mutation={mut_col}")

df[delta_col] = pd.to_numeric(df[delta_col], errors="coerce")
df[lig_col] = df[lig_col].astype(str).str.lower()

# Keep only xylobiose/xylotriose
df = df[df[lig_col].isin(["xylobiose", "xylotriose"])].copy()

# Group summaries per ligand
group_cols = [lig_col]
if org_col and gh_col:
    group_cols += [org_col, gh_col]

group = (
    df.groupby(group_cols, dropna=False)
      .agg(
          docking_rows=(protein_col, "count"),
          mean_delta=(delta_col, "mean"),
          median_delta=(delta_col, "median"),
          min_delta=(delta_col, "min"),
          max_delta=(delta_col, "max"),
          improved_or_retained_count=(delta_col, lambda x: (x <= 0).sum())
      )
      .reset_index()
)

group["improved_or_retained_percent"] = group["improved_or_retained_count"] / group["docking_rows"] * 100

rename = {
    lig_col: "ligand",
    delta_col: "delta_binding_mut_minus_wt"
}
if org_col:
    rename[org_col] = "organism_type"
if gh_col:
    rename[gh_col] = "gh_family"
group = group.rename(columns=rename)

for c in group.columns:
    if pd.api.types.is_float_dtype(group[c]):
        group[c] = group[c].round(4)

group.to_csv(OUT / "verified_modeller_ligand_group_summary.csv", index=False)

# Candidate-level ligand-specific best/worst tables
base_cols = [protein_col, mut_col, lig_col, delta_col]
if org_col:
    base_cols.insert(2, org_col)
if gh_col:
    base_cols.insert(3 if org_col else 2, gh_col)
if "foldx_ddg" in cols:
    base_cols.append("foldx_ddg")

cand = df[base_cols].copy()
cand = cand.rename(columns={
    protein_col: "candidate",
    mut_col: "mutation",
    lig_col: "ligand",
    delta_col: "delta_binding_mut_minus_wt",
    org_col if org_col else "": "organism_type",
    gh_col if gh_col else "": "gh_family"
})

for ligand in ["xylobiose", "xylotriose"]:
    sub = cand[cand["ligand"] == ligand].copy()

    best = sub.sort_values("delta_binding_mut_minus_wt", ascending=True).head(15)
    worst = sub.sort_values("delta_binding_mut_minus_wt", ascending=False).head(15)

    for table in [best, worst]:
        for c in table.columns:
            if pd.api.types.is_float_dtype(table[c]):
                table[c] = table[c].round(4)

    best.to_csv(OUT / f"verified_modeller_{ligand}_top15_favourable.csv", index=False)
    worst.to_csv(OUT / f"verified_modeller_{ligand}_top15_weaker.csv", index=False)

print("\n=== MODELLER LIGAND GROUP SUMMARY ===")
print(group.to_string(index=False))

for ligand in ["xylobiose", "xylotriose"]:
    print(f"\n=== MODELLER {ligand.upper()} TOP 15 FAVOURABLE ===")
    print(pd.read_csv(OUT / f"verified_modeller_{ligand}_top15_favourable.csv").to_string(index=False))

    print(f"\n=== MODELLER {ligand.upper()} TOP 15 WEAKER ===")
    print(pd.read_csv(OUT / f"verified_modeller_{ligand}_top15_weaker.csv").to_string(index=False))

print("\nSaved ligand-specific files in:", OUT)
