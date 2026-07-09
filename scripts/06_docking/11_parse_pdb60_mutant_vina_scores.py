#!/usr/bin/env python3
"""
Parse PDB60 mutant Vina-compatible docking scores.

This script:
1. Parses all REMARK VINA RESULT lines from mutant docking output PDBQT files.
2. Extracts all 10 modes per docking job.
3. Selects the best binding score per mutant-ligand job.
4. Merges docking scores with FoldX ΔΔG and mutation metadata.
5. Writes docking summary tables.

Important:
These are mutant-only docking scores. They should not yet be interpreted as docking improvement
until matching WT receptors are docked with the same fallback/Vina-compatible receptor preparation.

Input:
    results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/mutant_docking_vina_compatible_pdbqt/

Output:
    results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/mutant_docking_vina_compatible_pdbqt/parsed_scores/
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

            # Expected:
            # REMARK VINA RESULT: -7.2 0.000 0.000
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
                "vina_affinity_kcal_mol": affinity,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub,
            })

    return rows


def main():
    root = Path(".").resolve()

    base = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "mutant_docking_vina_compatible_pdbqt"
    )

    manifest_file = base / "pdb60_mutant_docking_jobs_manifest_VINA_COMPATIBLE.csv"

    if not manifest_file.exists():
        raise FileNotFoundError(f"Missing manifest: {manifest_file}")

    outdir = base / "parsed_scores"
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(manifest_file, low_memory=False)

    all_rows = []

    for _, r in manifest.iterrows():
        out_file = Path(str(r["vina_out_pdbqt"]))

        parsed = parse_vina_output(out_file)

        for row in parsed:
            row["job_id"] = r.get("job_id", "")
            row["structure_id"] = r.get("structure_id", "")
            row["uniprot_accession"] = r.get("uniprot_accession", "")
            row["pdb_id"] = r.get("pdb_id", "")
            row["organism_type"] = r.get("organism_type", "")
            row["gh_family"] = r.get("gh_family", "")
            row["foldx_mutation_code"] = r.get("foldx_mutation_code", "")
            row["ddg"] = r.get("ddg", np.nan)
            row["ligand_name"] = r.get("ligand_name", "")
            row["ligand_pdbqt"] = r.get("ligand_pdbqt", "")
            row["mutant_pdb_path"] = r.get("mutant_pdb_path", "")
            row["receptor_pdbqt"] = r.get("receptor_pdbqt", "")
            all_rows.append(row)

    all_modes = pd.DataFrame(all_rows)

    if all_modes.empty:
        raise RuntimeError("No Vina scores were parsed. Check output PDBQT files.")

    all_modes = all_modes.sort_values(
        ["ligand_name", "structure_id", "foldx_mutation_code", "vina_affinity_kcal_mol"],
        ascending=[True, True, True, True],
    )

    all_modes.to_csv(outdir / "pdb60_mutant_vina_all_modes.csv", index=False)

    best = (
        all_modes
        .sort_values("vina_affinity_kcal_mol", ascending=True)
        .groupby(["job_id", "ligand_name"], as_index=False)
        .head(1)
        .copy()
    )

    best = best.sort_values(["ligand_name", "vina_affinity_kcal_mol"], ascending=[True, True])

    best.to_csv(outdir / "pdb60_mutant_vina_best_scores.csv", index=False)

    # Summary by ligand.
    ligand_summary = (
        best
        .groupby("ligand_name")
        .agg(
            docking_jobs=("job_id", "count"),
            mean_affinity=("vina_affinity_kcal_mol", "mean"),
            median_affinity=("vina_affinity_kcal_mol", "median"),
            min_affinity=("vina_affinity_kcal_mol", "min"),
            max_affinity=("vina_affinity_kcal_mol", "max"),
            mean_ddg=("ddg", "mean"),
            median_ddg=("ddg", "median"),
        )
        .reset_index()
    )

    ligand_summary.to_csv(outdir / "pdb60_mutant_vina_summary_by_ligand.csv", index=False)

    group_summary = (
        best
        .groupby(["organism_type", "gh_family", "ligand_name"])
        .agg(
            docking_jobs=("job_id", "count"),
            unique_structures=("structure_id", "nunique"),
            mean_affinity=("vina_affinity_kcal_mol", "mean"),
            median_affinity=("vina_affinity_kcal_mol", "median"),
            best_affinity=("vina_affinity_kcal_mol", "min"),
            mean_ddg=("ddg", "mean"),
            median_ddg=("ddg", "median"),
        )
        .reset_index()
    )

    group_summary.to_csv(outdir / "pdb60_mutant_vina_summary_by_group.csv", index=False)

    # Top candidates separately by ligand.
    top_xylobiose = best[best["ligand_name"] == "xylobiose"].head(20)
    top_xylotriose = best[best["ligand_name"] == "xylotriose"].head(20)

    top_xylobiose.to_csv(outdir / "top20_pdb60_mutant_xylobiose_scores.csv", index=False)
    top_xylotriose.to_csv(outdir / "top20_pdb60_mutant_xylotriose_scores.csv", index=False)

    expected_jobs = len(manifest)
    successful_jobs = best.shape[0]
    modes_parsed = all_modes.shape[0]

    missing_jobs = expected_jobs - successful_jobs

    report = outdir / "PDB60_MUTANT_VINA_SCORE_PARSE_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 mutant Vina score parsing report\n\n")
        fh.write(f"- Expected mutant-ligand docking jobs: {expected_jobs}\n")
        fh.write(f"- Successful best-score jobs parsed: {successful_jobs}\n")
        fh.write(f"- Missing/failed jobs: {missing_jobs}\n")
        fh.write(f"- Total Vina modes parsed: {modes_parsed}\n")
        fh.write(f"- Expected modes if 10 per job: {expected_jobs * 10}\n\n")

        fh.write("## Ligand summary\n\n")
        fh.write(ligand_summary.to_string(index=False))
        fh.write("\n\n")

        fh.write("## Important interpretation note\n\n")
        fh.write(
            "These are mutant-only docking scores. They should not yet be interpreted as docking improvement "
            "because matching WT structures still need to be docked using the same Vina-compatible receptor "
            "PDBQT preparation and the same grid-generation method. Only after mutant-WT paired comparison "
            "can docking deltas be calculated.\n"
        )

    print("\n[DONE] PDB60 mutant Vina scores parsed.")
    print(f"Report: {report}")
    print(f"Best scores: {outdir / 'pdb60_mutant_vina_best_scores.csv'}")
    print(f"All modes: {outdir / 'pdb60_mutant_vina_all_modes.csv'}")
    print(f"Summary by ligand: {outdir / 'pdb60_mutant_vina_summary_by_ligand.csv'}")
    print(f"Summary by group: {outdir / 'pdb60_mutant_vina_summary_by_group.csv'}")


if __name__ == "__main__":
    main()
