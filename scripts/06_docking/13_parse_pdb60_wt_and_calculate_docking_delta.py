#!/usr/bin/env python3
"""
Parse matching WT docking scores and calculate paired mutant-WT docking deltas.

Inputs:
1. Mutant best scores:
   results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/mutant_docking_vina_compatible_pdbqt/parsed_scores/pdb60_mutant_vina_best_scores.csv

2. WT matching docking manifest:
   results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/wt_matching_docking_vina_compatible_pdbqt/pdb60_matching_wt_docking_jobs_manifest.csv

3. WT Vina outputs:
   results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/wt_matching_docking_vina_compatible_pdbqt/vina_outputs/

Output:
   results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/paired_mutant_wt_docking_delta/
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd
import numpy as np


def parse_vina_output(out_file: Path):
    rows = []

    if not out_file.exists():
        return rows

    mode = 0

    with out_file.open("r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("REMARK VINA RESULT"):
                continue

            parts = line.strip().split()

            try:
                affinity = float(parts[3])
            except Exception:
                continue

            try:
                rmsd_lb = float(parts[4])
            except Exception:
                rmsd_lb = np.nan

            try:
                rmsd_ub = float(parts[5])
            except Exception:
                rmsd_ub = np.nan

            mode += 1

            rows.append({
                "vina_out_pdbqt": str(out_file),
                "mode": mode,
                "wt_affinity_kcal_mol": affinity,
                "wt_rmsd_lb": rmsd_lb,
                "wt_rmsd_ub": rmsd_ub,
            })

    return rows


def main():
    root = Path(".").resolve()

    mutant_best_file = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "mutant_docking_vina_compatible_pdbqt"
        / "parsed_scores"
        / "pdb60_mutant_vina_best_scores.csv"
    )

    wt_base = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "wt_matching_docking_vina_compatible_pdbqt"
    )

    wt_manifest_file = wt_base / "pdb60_matching_wt_docking_jobs_manifest.csv"

    outdir = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "paired_mutant_wt_docking_delta"
    )

    outdir.mkdir(parents=True, exist_ok=True)

    if not mutant_best_file.exists():
        raise FileNotFoundError(f"Missing mutant best score file: {mutant_best_file}")

    if not wt_manifest_file.exists():
        raise FileNotFoundError(f"Missing WT manifest file: {wt_manifest_file}")

    mutant = pd.read_csv(mutant_best_file, low_memory=False)
    wt_manifest = pd.read_csv(wt_manifest_file, low_memory=False)

    wt_rows = []

    for _, r in wt_manifest.iterrows():
        out_file = Path(str(r["vina_out_pdbqt"]))

        parsed = parse_vina_output(out_file)

        for row in parsed:
            row["wt_job_id"] = r.get("wt_job_id", "")
            row["paired_mutant_job_id"] = r.get("paired_mutant_job_id", "")
            row["structure_id"] = r.get("structure_id", "")
            row["uniprot_accession"] = r.get("uniprot_accession", "")
            row["pdb_id"] = r.get("pdb_id", "")
            row["organism_type"] = r.get("organism_type", "")
            row["gh_family"] = r.get("gh_family", "")
            row["foldx_mutation_code"] = r.get("foldx_mutation_code", "")
            row["mutant_ddg"] = r.get("mutant_ddg", np.nan)
            row["ligand_name"] = r.get("ligand_name", "")
            row["wt_receptor_pdbqt"] = r.get("wt_receptor_pdbqt", "")
            row["wt_pdb_path"] = r.get("wt_pdb_path", "")
            wt_rows.append(row)

    wt_all_modes = pd.DataFrame(wt_rows)

    if wt_all_modes.empty:
        raise RuntimeError("No WT Vina scores were parsed.")

    wt_all_modes.to_csv(outdir / "pdb60_wt_vina_all_modes.csv", index=False)

    wt_mode_counts = (
        wt_all_modes
        .groupby(["paired_mutant_job_id", "ligand_name"])
        .agg(wt_modes_returned=("mode", "count"))
        .reset_index()
    )

    wt_mode_counts.to_csv(outdir / "pdb60_wt_mode_counts_by_job.csv", index=False)

    wt_best = (
        wt_all_modes
        .sort_values("wt_affinity_kcal_mol", ascending=True)
        .groupby(["paired_mutant_job_id", "ligand_name"], as_index=False)
        .head(1)
        .copy()
    )

    wt_best.to_csv(outdir / "pdb60_wt_vina_best_scores.csv", index=False)

    # Prepare mutant table.
    mutant = mutant.rename(
        columns={
            "vina_affinity_kcal_mol": "mutant_affinity_kcal_mol",
            "rmsd_lb": "mutant_rmsd_lb",
            "rmsd_ub": "mutant_rmsd_ub",
            "job_id": "mutant_job_id",
            "ddg": "foldx_ddg",
        }
    )

    keep_mutant_cols = [
        "mutant_job_id",
        "structure_id",
        "uniprot_accession",
        "pdb_id",
        "organism_type",
        "gh_family",
        "foldx_mutation_code",
        "foldx_ddg",
        "ligand_name",
        "mutant_affinity_kcal_mol",
        "mutant_rmsd_lb",
        "mutant_rmsd_ub",
        "mutant_pdb_path",
        "receptor_pdbqt",
    ]

    keep_mutant_cols = [c for c in keep_mutant_cols if c in mutant.columns]

    mutant_sub = mutant[keep_mutant_cols].copy()

    keep_wt_cols = [
        "paired_mutant_job_id",
        "ligand_name",
        "wt_job_id",
        "wt_affinity_kcal_mol",
        "wt_rmsd_lb",
        "wt_rmsd_ub",
        "wt_receptor_pdbqt",
        "wt_pdb_path",
    ]

    keep_wt_cols = [c for c in keep_wt_cols if c in wt_best.columns]

    wt_sub = wt_best[keep_wt_cols].copy()

    paired = mutant_sub.merge(
        wt_sub,
        left_on=["mutant_job_id", "ligand_name"],
        right_on=["paired_mutant_job_id", "ligand_name"],
        how="inner",
    )

    paired["delta_affinity_mut_minus_wt"] = (
        paired["mutant_affinity_kcal_mol"] - paired["wt_affinity_kcal_mol"]
    )

    paired["docking_improved_vs_wt"] = paired["delta_affinity_mut_minus_wt"] < 0

    paired = paired.sort_values(
        ["ligand_name", "delta_affinity_mut_minus_wt"],
        ascending=[True, True],
    )

    paired.to_csv(outdir / "pdb60_paired_mutant_wt_docking_delta.csv", index=False)

    # Summary by ligand.
    summary = (
        paired
        .groupby("ligand_name")
        .agg(
            paired_jobs=("mutant_job_id", "count"),
            improved_jobs=("docking_improved_vs_wt", "sum"),
            mean_delta=("delta_affinity_mut_minus_wt", "mean"),
            median_delta=("delta_affinity_mut_minus_wt", "median"),
            best_delta=("delta_affinity_mut_minus_wt", "min"),
            worst_delta=("delta_affinity_mut_minus_wt", "max"),
            mean_mutant_affinity=("mutant_affinity_kcal_mol", "mean"),
            mean_wt_affinity=("wt_affinity_kcal_mol", "mean"),
            mean_foldx_ddg=("foldx_ddg", "mean"),
            median_foldx_ddg=("foldx_ddg", "median"),
        )
        .reset_index()
    )

    summary["improved_percent"] = 100 * summary["improved_jobs"] / summary["paired_jobs"]

    summary.to_csv(outdir / "pdb60_paired_docking_delta_summary_by_ligand.csv", index=False)

    group_summary = (
        paired
        .groupby(["organism_type", "gh_family", "ligand_name"])
        .agg(
            paired_jobs=("mutant_job_id", "count"),
            improved_jobs=("docking_improved_vs_wt", "sum"),
            mean_delta=("delta_affinity_mut_minus_wt", "mean"),
            median_delta=("delta_affinity_mut_minus_wt", "median"),
            best_delta=("delta_affinity_mut_minus_wt", "min"),
            mean_foldx_ddg=("foldx_ddg", "mean"),
        )
        .reset_index()
    )

    group_summary["improved_percent"] = 100 * group_summary["improved_jobs"] / group_summary["paired_jobs"]

    group_summary.to_csv(outdir / "pdb60_paired_docking_delta_summary_by_group.csv", index=False)

    # Make wide paired candidate table.
    wide = paired.pivot_table(
        index=[
            "mutant_job_id",
            "structure_id",
            "uniprot_accession",
            "pdb_id",
            "organism_type",
            "gh_family",
            "foldx_mutation_code",
            "foldx_ddg",
        ],
        columns="ligand_name",
        values=[
            "mutant_affinity_kcal_mol",
            "wt_affinity_kcal_mol",
            "delta_affinity_mut_minus_wt",
        ],
        aggfunc="first",
    )

    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()

    # Both-ligand improvement.
    delta_cols = [c for c in wide.columns if c.startswith("delta_affinity_mut_minus_wt_")]

    if "delta_affinity_mut_minus_wt_xylobiose" in wide.columns:
        wide["xylobiose_improved"] = wide["delta_affinity_mut_minus_wt_xylobiose"] < 0
    else:
        wide["xylobiose_improved"] = False

    if "delta_affinity_mut_minus_wt_xylotriose" in wide.columns:
        wide["xylotriose_improved"] = wide["delta_affinity_mut_minus_wt_xylotriose"] < 0
    else:
        wide["xylotriose_improved"] = False

    wide["both_ligands_improved"] = wide["xylobiose_improved"] & wide["xylotriose_improved"]

    if delta_cols:
        wide["mean_delta_affinity"] = wide[delta_cols].mean(axis=1)
    else:
        wide["mean_delta_affinity"] = np.nan

    wide = wide.sort_values(["both_ligands_improved", "mean_delta_affinity", "foldx_ddg"], ascending=[False, True, True])

    wide.to_csv(outdir / "pdb60_paired_docking_delta_wide_candidate_table.csv", index=False)

    top_both = wide[wide["both_ligands_improved"]].head(30).copy()
    top_both.to_csv(outdir / "top_pdb60_candidates_improved_both_ligands.csv", index=False)

    expected_wt_jobs = len(wt_manifest)
    successful_wt_best = wt_best.shape[0]
    paired_jobs = paired.shape[0]

    wt_less_than_10 = (wt_mode_counts["wt_modes_returned"] < 10).sum()

    report = outdir / "PDB60_PAIRED_DOCKING_DELTA_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 paired mutant-WT docking delta report\n\n")
        fh.write(f"- Expected WT docking jobs: {expected_wt_jobs}\n")
        fh.write(f"- WT jobs with best score parsed: {successful_wt_best}\n")
        fh.write(f"- WT jobs with fewer than 10 returned modes: {wt_less_than_10}\n")
        fh.write(f"- Total WT modes parsed: {wt_all_modes.shape[0]}\n")
        fh.write(f"- Paired mutant-WT comparisons: {paired_jobs}\n\n")

        fh.write("## Summary by ligand\n\n")
        fh.write(summary.to_string(index=False))
        fh.write("\n\n")

        fh.write("## Interpretation rule\n\n")
        fh.write(
            "delta_affinity_mut_minus_wt = mutant affinity - WT affinity. "
            "A negative value means the mutant docked more strongly than the paired WT structure under the same receptor-preparation and grid settings. "
            "A positive value means the mutant docked weaker than WT.\n\n"
        )

        fh.write("## Important limitation\n\n")
        fh.write(
            "Both mutant and WT docking were performed using the same Vina-compatible fallback receptor PDBQT preparation with zero charges. "
            "Therefore, the paired delta is more meaningful than the absolute docking score, but the result should still be treated as comparative screening evidence rather than precise binding-energy estimation.\n"
        )

    print("\n[DONE] WT docking parsed and paired mutant-WT docking deltas calculated.")
    print(f"Report: {report}")
    print(f"Paired delta table: {outdir / 'pdb60_paired_mutant_wt_docking_delta.csv'}")
    print(f"Wide candidate table: {outdir / 'pdb60_paired_docking_delta_wide_candidate_table.csv'}")
    print(f"Top both-ligand candidates: {outdir / 'top_pdb60_candidates_improved_both_ligands.csv'}")


if __name__ == "__main__":
    main()
