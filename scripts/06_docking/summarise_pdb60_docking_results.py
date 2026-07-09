#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path("results/optionC_original_only/pdb60_foldx_mutations")
OUTDIR = BASE / "pdb60_results_summary"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Main expected files
ddg_all = BASE / "parsed_ddg_and_docking/pdb60_foldx_ddg_parsed_all.csv"
stab_all = BASE / "parsed_ddg_and_docking/pdb60_foldx_stabilizing_mutants.csv"
top_for_docking = BASE / "parsed_ddg_and_docking/pdb60_top_stabilizing_mutants_for_docking.csv"

vina_dir = BASE / "parsed_ddg_and_docking/mutant_docking_vina_compatible_pdbqt/parsed_scores"
vina_best = vina_dir / "pdb60_mutant_vina_best_scores.csv"
vina_all = vina_dir / "pdb60_mutant_vina_all_modes.csv"
vina_group = vina_dir / "pdb60_mutant_vina_summary_by_group.csv"
vina_ligand = vina_dir / "pdb60_mutant_vina_summary_by_ligand.csv"
top_x2 = vina_dir / "top20_pdb60_mutant_xylobiose_scores.csv"
top_x3 = vina_dir / "top20_pdb60_mutant_xylotriose_scores.csv"

report_lines = []
report_lines.append("# PDB60 FoldX and docking result summary\n")

def read_csv(path):
    if not path.exists():
        report_lines.append(f"Missing file: {path}")
        return None
    df = pd.read_csv(path)
    report_lines.append(f"Loaded {path}: shape={df.shape}")
    return df

ddg = read_csv(ddg_all)
stab = read_csv(stab_all)
top = read_csv(top_for_docking)
best = read_csv(vina_best)
allm = read_csv(vina_all)
grp = read_csv(vina_group)
lig = read_csv(vina_ligand)
x2 = read_csv(top_x2)
x3 = read_csv(top_x3)

# ---------- FoldX PDB60 summaries ----------
if ddg is not None:
    # Find ddg column
    ddg_col = None
    for c in ddg.columns:
        if c.lower() in ["ddg", "Δδg", "delta_delta_g"] or "ddg" in c.lower():
            ddg_col = c
            break

    if ddg_col:
        ddg[ddg_col] = pd.to_numeric(ddg[ddg_col], errors="coerce")
        ddg["effect_calc"] = np.where(ddg[ddg_col] < 0, "stabilising",
                              np.where(ddg[ddg_col] > 0, "destabilising", "zero"))

        ddg[ddg_col].describe().to_csv(OUTDIR / "pdb60_foldx_ddg_distribution.csv")

        effect_counts = ddg["effect_calc"].value_counts().rename_axis("effect").reset_index(name="count")
        effect_counts["percentage"] = effect_counts["count"] / len(ddg) * 100
        effect_counts.to_csv(OUTDIR / "pdb60_foldx_effect_counts.csv", index=False)

        if "organism_type" in ddg.columns and "gh_family" in ddg.columns:
            by_group = (
                ddg.groupby(["organism_type", "gh_family", "effect_calc"])
                   .size()
                   .unstack(fill_value=0)
                   .reset_index()
            )
            for col in ["stabilising", "destabilising", "zero"]:
                if col not in by_group.columns:
                    by_group[col] = 0
            by_group["total"] = by_group["stabilising"] + by_group["destabilising"] + by_group["zero"]
            by_group["stabilising_percentage"] = by_group["stabilising"] / by_group["total"] * 100
            by_group.to_csv(OUTDIR / "pdb60_foldx_effect_by_group.csv", index=False)

        useful_cols = [c for c in [
            "pdb_id", "pdb_tag", "structure_id", "uniprot_accession",
            "organism_type", "gh_family", "mutation", "mutation_code", ddg_col
        ] if c in ddg.columns]
        ddg.sort_values(ddg_col).head(30)[useful_cols].to_csv(
            OUTDIR / "pdb60_top30_stabilising_from_all_parsed_ddg.csv", index=False
        )

# Top stabilising docking subset distribution
if top is not None:
    ddg_col = None
    for c in top.columns:
        if c.lower() in ["ddg", "Δδg", "delta_delta_g"] or "ddg" in c.lower():
            ddg_col = c
            break
    if ddg_col:
        top[ddg_col] = pd.to_numeric(top[ddg_col], errors="coerce")
        top[ddg_col].describe().to_csv(OUTDIR / "pdb60_top_docking_subset_ddg_distribution.csv")

        useful_cols = [c for c in [
            "pdb_id", "pdb_tag", "structure_id", "uniprot_accession",
            "organism_type", "gh_family", "mutation", "mutation_code", ddg_col
        ] if c in top.columns]
        top.sort_values(ddg_col).head(30)[useful_cols].to_csv(
            OUTDIR / "pdb60_top30_stabilising_docking_subset.csv", index=False
        )

# ---------- Docking summaries ----------
if best is not None:
    # Guess key columns
    score_col = None
    ligand_col = None
    receptor_col = None

    for c in best.columns:
        lc = c.lower()
        if score_col is None and ("score" in lc or "affinity" in lc or "binding" in lc):
            score_col = c
        if ligand_col is None and "ligand" in lc:
            ligand_col = c
        if receptor_col is None and ("receptor" in lc or "pdb" in lc or "structure" in lc):
            receptor_col = c

    report_lines.append("\nBest-score columns:")
    report_lines.append(str(list(best.columns)))
    report_lines.append(f"Detected score_col={score_col}, ligand_col={ligand_col}, receptor_col={receptor_col}")

    if score_col:
        best[score_col] = pd.to_numeric(best[score_col], errors="coerce")

        # Overall best-score summary
        best[score_col].describe().to_csv(OUTDIR / "pdb60_mutant_vina_best_score_distribution.csv")

        # By ligand
        if ligand_col:
            by_ligand = best.groupby(ligand_col)[score_col].agg(
                count="count",
                mean="mean",
                std="std",
                min="min",
                median="median",
                max="max"
            ).reset_index()
            by_ligand.to_csv(OUTDIR / "pdb60_mutant_vina_best_score_by_ligand.csv", index=False)

        # By organism/GH if columns exist
        if "organism_type" in best.columns and "gh_family" in best.columns:
            by_group = best.groupby(["organism_type", "gh_family", ligand_col])[score_col].agg(
                count="count",
                mean="mean",
                std="std",
                min="min",
                median="median",
                max="max"
            ).reset_index()
            by_group.to_csv(OUTDIR / "pdb60_mutant_vina_best_score_by_group_and_ligand.csv", index=False)

        # Top scores by ligand
        if ligand_col:
            for lig_name, sub in best.groupby(ligand_col):
                safe = str(lig_name).replace("/", "_").replace(" ", "_")
                sub.sort_values(score_col).head(30).to_csv(
                    OUTDIR / f"pdb60_top30_mutant_vina_{safe}.csv", index=False
                )

# Copy existing summary files into output dir if available
for src in [vina_group, vina_ligand, top_x2, top_x3]:
    if src is not None and src.exists():
        try:
            pd.read_csv(src).to_csv(OUTDIR / src.name, index=False)
        except Exception:
            pass

# Save report
(OUTDIR / "PDB60_RESULTS_SUMMARY_REPORT.md").write_text("\n".join(report_lines))

print("\nSaved PDB60 result summaries to:")
print(OUTDIR)

print("\nGenerated files:")
for p in sorted(OUTDIR.glob("*")):
    print(" -", p)

print("\nQuick preview:")
for fname in [
    "pdb60_foldx_effect_counts.csv",
    "pdb60_foldx_effect_by_group.csv",
    "pdb60_top_docking_subset_ddg_distribution.csv",
    "pdb60_mutant_vina_best_score_by_ligand.csv",
    "pdb60_mutant_vina_best_score_by_group_and_ligand.csv",
]:
    p = OUTDIR / fname
    if p.exists():
        print("\n===", fname, "===")
        print(pd.read_csv(p).head(20).to_string(index=False))
