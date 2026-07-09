#!/usr/bin/env python3
"""
Create final PDB60 ranked candidate table without tiers.

No Tier A/B/C labels are used.

Ranking logic:
- FoldX ΔΔG remains the primary stability evidence.
- Docking deltas are secondary paired comparative evidence.
- Candidates are ranked by a transparent combined rank score:
    rank(FoldX ΔΔG) +
    rank(mean docking delta) +
    rank(xylobiose docking delta) +
    rank(xylotriose docking delta)

Lower score = better candidate.

The table keeps the actual evidence columns so the reader can see why each candidate ranks where it does.
"""

from pathlib import Path
import pandas as pd
import numpy as np


def main():
    root = Path(".").resolve()

    in_file = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "paired_mutant_wt_docking_delta"
        / "pdb60_paired_docking_delta_wide_candidate_table.csv"
    )

    if not in_file.exists():
        raise FileNotFoundError(f"Missing input: {in_file}")

    outdir = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_final_ranked_no_tiers"
    )

    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_file, low_memory=False)

    # Remove duplicated columns from previous pivot/merge operations.
    df = df.loc[:, ~df.columns.duplicated()].copy()

    required = [
        "mutant_job_id",
        "structure_id",
        "uniprot_accession",
        "pdb_id",
        "organism_type",
        "gh_family",
        "foldx_mutation_code",
        "foldx_ddg",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")

    # Recompute docking deltas directly from mutant and WT affinities.
    if (
        "mutant_affinity_kcal_mol_xylobiose" in df.columns
        and "wt_affinity_kcal_mol_xylobiose" in df.columns
    ):
        df["delta_affinity_mut_minus_wt_xylobiose"] = (
            pd.to_numeric(df["mutant_affinity_kcal_mol_xylobiose"], errors="coerce")
            - pd.to_numeric(df["wt_affinity_kcal_mol_xylobiose"], errors="coerce")
        )
    else:
        df["delta_affinity_mut_minus_wt_xylobiose"] = np.nan

    if (
        "mutant_affinity_kcal_mol_xylotriose" in df.columns
        and "wt_affinity_kcal_mol_xylotriose" in df.columns
    ):
        df["delta_affinity_mut_minus_wt_xylotriose"] = (
            pd.to_numeric(df["mutant_affinity_kcal_mol_xylotriose"], errors="coerce")
            - pd.to_numeric(df["wt_affinity_kcal_mol_xylotriose"], errors="coerce")
        )
    else:
        df["delta_affinity_mut_minus_wt_xylotriose"] = np.nan

    df["mean_delta_affinity"] = df[
        [
            "delta_affinity_mut_minus_wt_xylobiose",
            "delta_affinity_mut_minus_wt_xylotriose",
        ]
    ].mean(axis=1)

    df["xylobiose_improved"] = df["delta_affinity_mut_minus_wt_xylobiose"] < 0
    df["xylotriose_improved"] = df["delta_affinity_mut_minus_wt_xylotriose"] < 0
    df["both_ligands_improved"] = df["xylobiose_improved"] & df["xylotriose_improved"]
    df["any_ligand_improved"] = df["xylobiose_improved"] | df["xylotriose_improved"]

    # Ranking components.
    # More negative FoldX and docking deltas are better, so ascending rank is used.
    df["rank_foldx_ddg"] = df["foldx_ddg"].rank(method="min", ascending=True)
    df["rank_xylobiose_delta"] = df["delta_affinity_mut_minus_wt_xylobiose"].rank(method="min", ascending=True)
    df["rank_xylotriose_delta"] = df["delta_affinity_mut_minus_wt_xylotriose"].rank(method="min", ascending=True)
    df["rank_mean_delta"] = df["mean_delta_affinity"].rank(method="min", ascending=True)

    df["final_rank_score"] = (
        df["rank_foldx_ddg"] +
        df["rank_xylobiose_delta"] +
        df["rank_xylotriose_delta"] +
        df["rank_mean_delta"]
    )

    # No tiers. Just a transparent ranking.
    # Candidates with both-ligand improvement are placed first because they satisfy both substrate checks.
    df = df.sort_values(
        [
            "both_ligands_improved",
            "final_rank_score",
            "foldx_ddg",
            "mean_delta_affinity",
        ],
        ascending=[False, True, True, True],
    ).copy()

    df["final_rank"] = range(1, len(df) + 1)

    final_cols = [
        "final_rank",
        "final_rank_score",
        "mutant_job_id",
        "structure_id",
        "uniprot_accession",
        "pdb_id",
        "organism_type",
        "gh_family",
        "foldx_mutation_code",
        "foldx_ddg",
        "delta_affinity_mut_minus_wt_xylobiose",
        "delta_affinity_mut_minus_wt_xylotriose",
        "mean_delta_affinity",
        "mutant_affinity_kcal_mol_xylobiose",
        "wt_affinity_kcal_mol_xylobiose",
        "mutant_affinity_kcal_mol_xylotriose",
        "wt_affinity_kcal_mol_xylotriose",
        "xylobiose_improved",
        "xylotriose_improved",
        "both_ligands_improved",
        "any_ligand_improved",
        "rank_foldx_ddg",
        "rank_xylobiose_delta",
        "rank_xylotriose_delta",
        "rank_mean_delta",
    ]

    final_cols = [c for c in final_cols if c in df.columns]

    final = df[final_cols].copy()

    final.to_csv(outdir / "pdb60_final_ranked_no_tiers.csv", index=False)

    top15 = final.head(15).copy()
    top15.to_csv(outdir / "pdb60_top15_ranked_no_tiers.csv", index=False)

    top20 = final.head(20).copy()
    top20.to_csv(outdir / "pdb60_top20_ranked_no_tiers.csv", index=False)

    summary = {
        "total_candidates": len(final),
        "both_ligands_improved": int(final["both_ligands_improved"].sum()),
        "at_least_one_ligand_improved": int(final["any_ligand_improved"].sum()),
        "xylobiose_improved": int(final["xylobiose_improved"].sum()),
        "xylotriose_improved": int(final["xylotriose_improved"].sum()),
        "mean_foldx_ddg": float(final["foldx_ddg"].mean()),
        "median_foldx_ddg": float(final["foldx_ddg"].median()),
        "mean_delta_affinity": float(final["mean_delta_affinity"].mean()),
        "median_delta_affinity": float(final["mean_delta_affinity"].median()),
    }

    pd.DataFrame([summary]).to_csv(outdir / "pdb60_no_tiers_summary.csv", index=False)

    group_summary = (
        final
        .groupby(["organism_type", "gh_family"])
        .agg(
            candidates=("mutant_job_id", "count"),
            both_ligands_improved=("both_ligands_improved", "sum"),
            at_least_one_ligand_improved=("any_ligand_improved", "sum"),
            mean_foldx_ddg=("foldx_ddg", "mean"),
            mean_delta_affinity=("mean_delta_affinity", "mean"),
        )
        .reset_index()
    )

    group_summary.to_csv(outdir / "pdb60_no_tiers_summary_by_group.csv", index=False)

    report = outdir / "PDB60_FINAL_RANKED_NO_TIERS_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 final ranked candidate report without tiers\n\n")
        fh.write("- No Tier A/B/C classification was used.\n")
        fh.write("- Candidates are ranked using a transparent combined rank score.\n")
        fh.write("- FoldX ΔΔG remains the primary stability evidence.\n")
        fh.write("- Paired mutant-WT docking deltas are secondary comparative evidence.\n\n")

        fh.write("## Summary\n\n")
        for k, v in summary.items():
            fh.write(f"- {k}: {v}\n")

        fh.write("\n## Ranking formula\n\n")
        fh.write(
            "final_rank_score = rank(FoldX ΔΔG) + rank(xylobiose Δ docking) + "
            "rank(xylotriose Δ docking) + rank(mean Δ docking).\n\n"
        )

        fh.write("Lower final_rank_score indicates stronger combined computational support.\n\n")

        fh.write("## Top 15 ranked candidates\n\n")

        display_cols = [
            "final_rank",
            "final_rank_score",
            "mutant_job_id",
            "organism_type",
            "gh_family",
            "foldx_ddg",
            "delta_affinity_mut_minus_wt_xylobiose",
            "delta_affinity_mut_minus_wt_xylotriose",
            "mean_delta_affinity",
            "both_ligands_improved",
        ]

        fh.write(top15[display_cols].to_string(index=False))
        fh.write("\n\n")

        fh.write("## Interpretation note\n\n")
        fh.write(
            "The final rank should be interpreted as a prioritisation score, not as a biological class. "
            "A candidate with negative FoldX ΔΔG and negative docking deltas for both ligands is more balanced, "
            "but docking results remain comparative screening evidence because both WT and mutant receptor PDBQT files "
            "were generated using the same fallback zero-charge receptor preparation.\n"
        )

    print("[DONE] Final non-tiered ranked PDB60 candidate table created.")
    print(f"Report: {report}")
    print(f"Final table: {outdir / 'pdb60_final_ranked_no_tiers.csv'}")
    print(f"Top 15: {outdir / 'pdb60_top15_ranked_no_tiers.csv'}")


if __name__ == "__main__":
    main()
