#!/usr/bin/env python3
# Purpose: Summarise MODELLER docking.

from pathlib import Path
import pandas as pd
import numpy as np

OUTDIR = Path("results/modeller_docking_results_summary")
OUTDIR.mkdir(parents=True, exist_ok=True)

candidate_files = [
    Path("results/foldx_clean/tier2_ddg_ranked_annotated.csv"),
    Path("results/foldx/tier2_ddg_ranked_annotated.csv"),
    Path("results/docking/top15_candidates_with_docking.csv"),
    Path("results/docking/top15_final_ranked.csv"),
    Path("results/docking/modeller_docking_delta_summary.csv"),
    Path("results/docking/modeller_docking_paired_summary.csv"),
    Path("results/docking/modeller_docking_all_scores.csv"),
    Path("results/docking/modeller_docking_best_scores.csv"),
    Path("results/docking/docking_delta_summary.csv"),
    Path("results/docking/docking_paired_summary.csv"),
    Path("results/docking/all_scores.csv"),
    Path("results/docking/best_scores.csv"),
    Path("results/final_candidate_selection/top15_candidates_with_docking.csv"),
    Path("results/final_candidate_selection/top15_final_ranked.csv"),
]

# Also search likely files automatically
auto_files = []
for root in ["results"]:
    for p in Path(root).rglob("*.csv"):
        s = str(p).lower()
        if any(k in s for k in ["docking", "vina", "candidate", "ddg", "delta"]):
            auto_files.append(p)

files = []
seen = set()
for p in candidate_files + auto_files:
    if p.exists() and p not in seen:
        files.append(p)
        seen.add(p)

report = []
report.append("# MODELLER docking result summary\n")
report.append("## Candidate CSV files found\n")
for p in files:
    report.append(f"- {p}")

def read_csv(p):
    try:
        return pd.read_csv(p)
    except Exception:
        return None

def find_col(cols, patterns):
    for pat in patterns:
        for c in cols:
            if pat in c.lower():
                return c
    return None

summary_rows = []

for p in files:
    df = read_csv(p)
    if df is None or df.empty:
        continue

    cols = list(df.columns)
    low = [c.lower() for c in cols]

    interesting = any(
        any(k in c for k in ["ddg", "dock", "vina", "xylobiose", "xylotriose", "candidate", "mutation", "delta"])
        for c in low
    )
    if not interesting:
        continue

    report.append("\n" + "="*100)
    report.append(f"FILE: {p}")
    report.append(f"shape: {df.shape}")
    report.append(f"columns: {cols}")
    report.append("\nPreview:")
    report.append(df.head(10).to_string(index=False))

    # Save a simple metadata row
    summary_rows.append({
        "file": str(p),
        "rows": df.shape[0],
        "columns": df.shape[1],
        "column_names": "; ".join(cols)
    })

    # Identify useful columns
    ddg_col = find_col(cols, ["ddg"])
    mutation_col = find_col(cols, ["mutation"])
    accession_col = find_col(cols, ["uniprot_accession", "accession", "protein"])
    organism_col = find_col(cols, ["organism_type"])
    gh_col = find_col(cols, ["gh_family"])
    x2_col = find_col(cols, ["xylobiose"])
    x3_col = find_col(cols, ["xylotriose"])
    delta_col = find_col(cols, ["mean_docking_delta", "mean_delta", "docking_delta", "delta"])

    # FoldX effect summaries
    if ddg_col:
        df[ddg_col] = pd.to_numeric(df[ddg_col], errors="coerce")
        report.append(f"\nDDG summary using column: {ddg_col}")
        report.append(df[ddg_col].describe().to_string())

        tmp = df.copy()
        tmp["effect_calc"] = np.where(tmp[ddg_col] < 0, "stabilising",
                              np.where(tmp[ddg_col] > 0, "destabilising", "zero"))
        report.append("\nEffect counts:")
        report.append(tmp["effect_calc"].value_counts().to_string())

        if organism_col and gh_col:
            by_group = pd.crosstab([tmp[organism_col], tmp[gh_col]], tmp["effect_calc"])
            report.append("\nEffect by organism type and GH family:")
            report.append(by_group.to_string())

        useful = [c for c in [accession_col, organism_col, gh_col, mutation_col, ddg_col] if c]
        if useful:
            top_stab = df.sort_values(ddg_col).head(30)[useful]
            top_stab.to_csv(OUTDIR / f"{p.stem}_top30_stabilising.csv", index=False)
            report.append("\nTop 30 stabilising by FoldX DDG:")
            report.append(top_stab.to_string(index=False))

    # Docking summaries if xylobiose/xylotriose columns exist
    numeric_cols = []
    for c in cols:
        lc = c.lower()
        if any(k in lc for k in ["xylobiose", "xylotriose", "delta", "dock", "vina", "score", "affinity"]):
            df[c] = pd.to_numeric(df[c], errors="coerce")
            if df[c].notna().sum() > 0:
                numeric_cols.append(c)

    if numeric_cols:
        report.append("\nDocking-related numeric summaries:")
        report.append(df[numeric_cols].describe().to_string())

    # If paired xylobiose/xylotriose deltas exist, rank candidates
    # Prefer columns containing delta and ligand names
    x2_delta_cols = [c for c in cols if "xylobiose" in c.lower() and ("delta" in c.lower() or "change" in c.lower())]
    x3_delta_cols = [c for c in cols if "xylotriose" in c.lower() and ("delta" in c.lower() or "change" in c.lower())]

    if x2_delta_cols and x3_delta_cols:
        x2d = x2_delta_cols[0]
        x3d = x3_delta_cols[0]
        df[x2d] = pd.to_numeric(df[x2d], errors="coerce")
        df[x3d] = pd.to_numeric(df[x3d], errors="coerce")
        df["mean_docking_change_calc"] = df[[x2d, x3d]].mean(axis=1)

        useful = [c for c in [accession_col, organism_col, gh_col, mutation_col, ddg_col, x2d, x3d, "mean_docking_change_calc"] if c and c in df.columns]

        best = df.sort_values("mean_docking_change_calc").head(30)[useful]
        worst = df.sort_values("mean_docking_change_calc", ascending=False).head(30)[useful]

        best.to_csv(OUTDIR / f"{p.stem}_top30_best_docking_delta.csv", index=False)
        worst.to_csv(OUTDIR / f"{p.stem}_top30_worst_docking_delta.csv", index=False)

        report.append(f"\nDetected paired ligand delta columns: {x2d}, {x3d}")
        report.append("\nTop 30 best docking changes:")
        report.append(best.to_string(index=False))
        report.append("\nTop 30 worst docking changes:")
        report.append(worst.to_string(index=False))

    # If a generic delta column exists, rank candidates
    elif delta_col:
        df[delta_col] = pd.to_numeric(df[delta_col], errors="coerce")
        useful = [c for c in [accession_col, organism_col, gh_col, mutation_col, ddg_col, delta_col] if c and c in df.columns]

        best = df.sort_values(delta_col).head(30)[useful]
        worst = df.sort_values(delta_col, ascending=False).head(30)[useful]

        best.to_csv(OUTDIR / f"{p.stem}_top30_best_generic_delta.csv", index=False)
        worst.to_csv(OUTDIR / f"{p.stem}_top30_worst_generic_delta.csv", index=False)

        report.append(f"\nDetected generic delta column: {delta_col}")
        report.append("\nTop 30 best generic docking delta:")
        report.append(best.to_string(index=False))
        report.append("\nTop 30 worst generic docking delta:")
        report.append(worst.to_string(index=False))

    # If raw xylobiose/xylotriose scores exist but no deltas
    elif x2_col and x3_col:
        df[x2_col] = pd.to_numeric(df[x2_col], errors="coerce")
        df[x3_col] = pd.to_numeric(df[x3_col], errors="coerce")
        df["mean_ligand_score_calc"] = df[[x2_col, x3_col]].mean(axis=1)

        useful = [c for c in [accession_col, organism_col, gh_col, mutation_col, ddg_col, x2_col, x3_col, "mean_ligand_score_calc"] if c and c in df.columns]

        best = df.sort_values("mean_ligand_score_calc").head(30)[useful]
        worst = df.sort_values("mean_ligand_score_calc", ascending=False).head(30)[useful]

        best.to_csv(OUTDIR / f"{p.stem}_top30_best_raw_scores.csv", index=False)
        worst.to_csv(OUTDIR / f"{p.stem}_top30_worst_raw_scores.csv", index=False)

        report.append("\nTop 30 best raw docking scores:")
        report.append(best.to_string(index=False))
        report.append("\nTop 30 worst raw docking scores:")
        report.append(worst.to_string(index=False))

# Save file inventory
pd.DataFrame(summary_rows).to_csv(OUTDIR / "modeller_docking_candidate_files_inventory.csv", index=False)

# Save report
report_path = OUTDIR / "MODELLER_DOCKING_RESULTS_SUMMARY.md"
report_path.write_text("\n".join(report))

print("Saved:", report_path)
print("\nGenerated files:")
for p in sorted(OUTDIR.glob("*")):
    print(" -", p)

print("\nQuick grep command:")
print(f"grep -nEi 'FILE:|shape:|DDG summary|Effect counts|Effect by organism|Top 30 best|Top 30 worst|Detected' {report_path}")
