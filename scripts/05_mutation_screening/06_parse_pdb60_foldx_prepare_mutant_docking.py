#!/usr/bin/env python3
"""
Parse PDB60 FoldX BuildModel ΔΔG results and prepare mutant docking jobs.

This script:
1. Parses FoldX Dif_*.fxout files.
2. Maps FoldX model indices back to mutation codes.
3. Creates a ranked ΔΔG table.
4. Selects top stabilizing mutants per structure.
5. Prepares receptor PDBQT conversion script.
6. Prepares AutoDock Vina mutant docking script for xylobiose and xylotriose.

Run:
    cd ~/xylanase-thesis
    python scripts/10_optionC_original_only/06_parse_pdb60_foldx_prepare_mutant_docking.py

Optional:
    python scripts/10_optionC_original_only/06_parse_pdb60_foldx_prepare_mutant_docking.py --top-per-structure 3
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def safe_id(x):
    return (
        str(x)
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
        .replace(";", "")
        .replace(":", "_")
    )


def first_existing(paths):
    for f in paths:
        p = Path(f)
        if p.exists():
            return p
    return None


def find_ligand(root: Path, names: list[str]) -> str:
    search_dirs = [
        root / "data" / "ligands",
        root / "data" / "raw" / "ligands",
        root / "ligands",
        root / "docking" / "ligands",
        root / "results" / "docking" / "ligands",
    ]

    for d in search_dirs:
        if not d.exists():
            continue

        for name in names:
            p = d / name
            if p.exists():
                return str(p)

    return ""


def parse_first_float(tokens):
    for t in tokens:
        try:
            return float(t)
        except Exception:
            continue
    return np.nan


def parse_model_index(model_name: str):
    """
    Tries to extract FoldX mutant index from model name.

    Examples:
    1ABC_A_1.pdb -> 1
    1ABC_A_12.pdb -> 12
    1ABC_A_Repair_3.pdb -> 3
    """

    m = re.search(r"_(\d+)(?:_\d+)?\.pdb$", str(model_name))

    if m:
        return int(m.group(1))

    return None


def parse_dif_file(dif_file: Path) -> pd.DataFrame:
    rows = []

    with dif_file.open("r", errors="ignore") as fh:
        line_counter = 0

        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            if line.lower().startswith("pdb"):
                continue

            if line.startswith("---"):
                continue

            parts = line.replace("\t", " ").split()

            if len(parts) < 2:
                continue

            model_name = parts[0]

            ddg = parse_first_float(parts[1:])

            if pd.isna(ddg):
                continue

            line_counter += 1

            model_index = parse_model_index(model_name)

            if model_index is None:
                model_index = line_counter

            rows.append({
                "dif_file": str(dif_file),
                "foldx_model_name": model_name,
                "foldx_model_index": model_index,
                "ddg": ddg,
                "line_order": line_counter,
            })

    return pd.DataFrame(rows)


def pdb_center_and_box(pdb_path: str, padding: float = 8.0):
    p = Path(str(pdb_path))

    if not p.exists():
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    coords = []

    with p.open("r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except Exception:
                continue

    if not coords:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    arr = np.array(coords, dtype=float)

    mins = arr.min(axis=0)
    maxs = arr.max(axis=0)

    center = arr.mean(axis=0)
    size = (maxs - mins) + padding

    size = np.clip(size, 20.0, 60.0)

    return center[0], center[1], center[2], size[0], size[1], size[2]


def find_mutant_pdb(workdir: Path, model_name: str, model_index: int):
    direct = workdir / model_name

    if direct.exists() and direct.suffix.lower() == ".pdb":
        return str(direct)

    # Fallback: find PDBs ending with model index.
    candidates = []

    for p in workdir.glob("*.pdb"):
        name = p.name

        if name.endswith(f"_{model_index}.pdb"):
            candidates.append(p)

        elif re.search(rf"_{model_index}_\d+\.pdb$", name):
            candidates.append(p)

    if candidates:
        return str(sorted(candidates)[0])

    return ""


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--top-per-structure",
        type=int,
        default=1,
        help="Number of top stabilizing mutants per structure to send to docking.",
    )

    parser.add_argument(
        "--ddg-threshold",
        type=float,
        default=0.0,
        help="Only mutants with ddG below this threshold are selected for docking.",
    )

    args = parser.parse_args()

    root = Path(".").resolve()

    base = root / "results" / "optionC_original_only" / "pdb60_foldx_mutations"
    work_root = base / "foldx_buildmodel_runs"
    outdir = base / "parsed_ddg_and_docking"
    docking_dir = outdir / "mutant_docking"
    receptor_dir = docking_dir / "receptors_pdbqt"
    vina_out_dir = docking_dir / "vina_outputs"
    vina_log_dir = docking_dir / "vina_logs"

    for d in [outdir, docking_dir, receptor_dir, vina_out_dir, vina_log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    mutation_file = base / "pdb60_consensus_mutation_candidates.csv"

    if not mutation_file.exists():
        raise FileNotFoundError(f"Missing mutation candidate file: {mutation_file}")

    mutations = pd.read_csv(mutation_file, low_memory=False)

    mutations["structure_id"] = mutations["structure_id"].astype(str).map(safe_id)
    mutations["candidate_order"] = mutations.groupby("structure_id").cumcount() + 1

    parsed_parts = []

    for sid, sub in mutations.groupby("structure_id"):
        workdir = work_root / sid

        if not workdir.exists():
            continue

        dif_files = sorted(workdir.glob("Dif_*.fxout"))

        if not dif_files:
            continue

        # Usually one Dif file per structure.
        for dif in dif_files:
            parsed = parse_dif_file(dif)

            if parsed.empty:
                continue

            parsed["structure_id"] = sid
            parsed["workdir"] = str(workdir)
            parsed_parts.append(parsed)

    if not parsed_parts:
        raise RuntimeError("No FoldX Dif_*.fxout results could be parsed.")

    parsed_all = pd.concat(parsed_parts, ignore_index=True)

    merged = parsed_all.merge(
        mutations,
        left_on=["structure_id", "foldx_model_index"],
        right_on=["structure_id", "candidate_order"],
        how="left",
    )

    merged["mutant_pdb_path"] = merged.apply(
        lambda r: find_mutant_pdb(
            Path(r["workdir"]),
            str(r["foldx_model_name"]),
            int(r["foldx_model_index"]),
        ),
        axis=1,
    )

    merged["mutant_pdb_exists"] = merged["mutant_pdb_path"].map(
        lambda x: bool(x) and Path(str(x)).exists()
    )

    merged = merged.sort_values("ddg", ascending=True)

    merged.to_csv(outdir / "pdb60_foldx_ddg_parsed_all.csv", index=False)

    stabilizing = merged[merged["ddg"] < args.ddg_threshold].copy()

    stabilizing.to_csv(outdir / "pdb60_foldx_stabilizing_mutants.csv", index=False)

    top = (
        stabilizing
        .sort_values(["structure_id", "ddg"], ascending=[True, True])
        .groupby("structure_id")
        .head(args.top_per_structure)
        .copy()
    )

    top = top.sort_values("ddg", ascending=True)

    top.to_csv(outdir / "pdb60_top_stabilizing_mutants_for_docking.csv", index=False)

    ligand_xylobiose = find_ligand(
        root,
        [
            "xylobiose.pdbqt",
            "xylobiose_ligand.pdbqt",
            "XYLOBIOSE.pdbqt",
        ],
    )

    ligand_xylotriose = find_ligand(
        root,
        [
            "xylotriose.pdbqt",
            "xylotriose_ligand.pdbqt",
            "XYLOTRIOSE.pdbqt",
        ],
    )

    job_rows = []

    for _, r in top.iterrows():
        mutant_pdb = r["mutant_pdb_path"]

        if not mutant_pdb or not Path(str(mutant_pdb)).exists():
            continue

        mutation_label = str(r.get("foldx_mutation_code", f"model_{r['foldx_model_index']}")).replace(";", "")
        job_id = safe_id(f"{r['structure_id']}_{mutation_label}")

        receptor_pdbqt = receptor_dir / f"{job_id}.pdbqt"

        cx, cy, cz, sx, sy, sz = pdb_center_and_box(mutant_pdb)

        base_row = {
            "job_id": job_id,
            "structure_id": r["structure_id"],
            "uniprot_accession": r.get("uniprot_accession", ""),
            "pdb_id": r.get("pdb_id", ""),
            "organism_type": r.get("organism_type", ""),
            "gh_family": r.get("gh_family", ""),
            "foldx_mutation_code": r.get("foldx_mutation_code", ""),
            "ddg": r["ddg"],
            "mutant_pdb_path": mutant_pdb,
            "mutant_pdb_exists": Path(str(mutant_pdb)).exists(),
            "receptor_pdbqt": str(receptor_pdbqt),
            "center_x": cx,
            "center_y": cy,
            "center_z": cz,
            "size_x": sx,
            "size_y": sy,
            "size_z": sz,
        }

        for ligand_name, ligand_path in [
            ("xylobiose", ligand_xylobiose),
            ("xylotriose", ligand_xylotriose),
        ]:
            jr = base_row.copy()
            jr["ligand_name"] = ligand_name
            jr["ligand_pdbqt"] = ligand_path
            jr["ligand_pdbqt_exists"] = bool(ligand_path) and Path(ligand_path).exists()
            jr["vina_out_pdbqt"] = str(vina_out_dir / f"{job_id}_{ligand_name}_out.pdbqt")
            jr["vina_log"] = str(vina_log_dir / f"{job_id}_{ligand_name}.log")
            job_rows.append(jr)

    jobs = pd.DataFrame(job_rows)

    jobs.to_csv(docking_dir / "pdb60_mutant_docking_jobs_manifest.csv", index=False)

    prep_script = docking_dir / "01_prepare_pdb60_mutant_receptors_obabel.sh"

    with prep_script.open("w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("set -euo pipefail\n\n")
        fh.write("echo 'Preparing PDB60 mutant receptors using Open Babel...'\n\n")

        for mutant_pdb, receptor_pdbqt in (
            jobs[["mutant_pdb_path", "receptor_pdbqt"]]
            .drop_duplicates()
            .itertuples(index=False)
        ):
            fh.write(f"mkdir -p '{Path(receptor_pdbqt).parent}'\n")
            fh.write(
                f"if [ ! -f '{receptor_pdbqt}' ]; then\n"
                f"  obabel -ipdb '{mutant_pdb}' -opdbqt -O '{receptor_pdbqt}' --partialcharge gasteiger\n"
                f"fi\n\n"
            )

    prep_script.chmod(0o755)

    vina_script = docking_dir / "02_run_pdb60_mutant_vina.sh"

    with vina_script.open("w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("set -euo pipefail\n\n")
        fh.write("echo 'Running PDB60 mutant docking with AutoDock Vina...'\n\n")

        for _, r in jobs.iterrows():
            if not r["ligand_pdbqt"]:
                continue

            fh.write(f"mkdir -p '{Path(r['vina_out_pdbqt']).parent}' '{Path(r['vina_log']).parent}'\n")
            fh.write(
                "vina "
                f"--receptor '{r['receptor_pdbqt']}' "
                f"--ligand '{r['ligand_pdbqt']}' "
                f"--center_x {r['center_x']:.3f} "
                f"--center_y {r['center_y']:.3f} "
                f"--center_z {r['center_z']:.3f} "
                f"--size_x {r['size_x']:.3f} "
                f"--size_y {r['size_y']:.3f} "
                f"--size_z {r['size_z']:.3f} "
                "--exhaustiveness 8 "
                "--num_modes 10 "
                f"--out '{r['vina_out_pdbqt']}' "
                f"--log '{r['vina_log']}'\n\n"
            )

    vina_script.chmod(0o755)

    summary_rows = []

    summary_rows.append({
        "metric": "mutation_candidates_original",
        "value": len(mutations),
    })

    summary_rows.append({
        "metric": "parsed_foldx_ddg_rows",
        "value": len(merged),
    })

    summary_rows.append({
        "metric": "stabilizing_mutants_ddg_below_threshold",
        "value": len(stabilizing),
    })

    summary_rows.append({
        "metric": "structures_with_stabilizing_mutants",
        "value": stabilizing["structure_id"].nunique() if not stabilizing.empty else 0,
    })

    summary_rows.append({
        "metric": "top_mutants_selected_for_docking",
        "value": len(top),
    })

    summary_rows.append({
        "metric": "mutant_docking_jobs",
        "value": len(jobs),
    })

    summary_rows.append({
        "metric": "xylobiose_ligand_found",
        "value": bool(ligand_xylobiose),
    })

    summary_rows.append({
        "metric": "xylotriose_ligand_found",
        "value": bool(ligand_xylotriose),
    })

    pd.DataFrame(summary_rows).to_csv(outdir / "pdb60_foldx_parse_docking_summary.csv", index=False)

    report = outdir / "PDB60_FOLDX_DDG_AND_MUTANT_DOCKING_PREP_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 FoldX ΔΔG parsing and mutant docking preparation report\n\n")
        fh.write(f"- Original mutation candidates: {len(mutations)}\n")
        fh.write(f"- Parsed FoldX ΔΔG rows: {len(merged)}\n")
        fh.write(f"- Stabilizing mutants, ddG < {args.ddg_threshold}: {len(stabilizing)}\n")
        fh.write(f"- Structures with stabilizing mutants: {stabilizing['structure_id'].nunique() if not stabilizing.empty else 0}\n")
        fh.write(f"- Top mutants selected for docking: {len(top)}\n")
        fh.write(f"- Mutant docking jobs prepared: {len(jobs)}\n")
        fh.write(f"- Xylobiose ligand found: {bool(ligand_xylobiose)}\n")
        fh.write(f"- Xylotriose ligand found: {bool(ligand_xylotriose)}\n\n")
        fh.write("## Output files\n\n")
        fh.write("- `pdb60_foldx_ddg_parsed_all.csv`\n")
        fh.write("- `pdb60_foldx_stabilizing_mutants.csv`\n")
        fh.write("- `pdb60_top_stabilizing_mutants_for_docking.csv`\n")
        fh.write("- `mutant_docking/pdb60_mutant_docking_jobs_manifest.csv`\n")
        fh.write("- `mutant_docking/01_prepare_pdb60_mutant_receptors_obabel.sh`\n")
        fh.write("- `mutant_docking/02_run_pdb60_mutant_vina.sh`\n\n")
        fh.write("## Important note\n\n")
        fh.write(
            "Docking boxes are generated from receptor coordinates as whole-receptor fallback boxes. "
            "If the original docking workflow used catalytic-cleft grid boxes, replace these coordinates "
            "with the original grid-generation method before interpreting binding scores.\n"
        )

    print("\n[DONE] PDB60 FoldX ΔΔG parsed and mutant docking jobs prepared.")
    print(f"Report: {report}")
    print(f"Top mutants for docking: {outdir / 'pdb60_top_stabilizing_mutants_for_docking.csv'}")
    print(f"Docking manifest: {docking_dir / 'pdb60_mutant_docking_jobs_manifest.csv'}")
    print(f"Prepare receptors: {prep_script}")
    print(f"Run Vina: {vina_script}")


if __name__ == "__main__":
    main()
