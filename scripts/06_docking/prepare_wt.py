#!/usr/bin/env python3
# Purpose: Prepare wild-type docking.

from __future__ import annotations

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


def atom_type_from_protein_atom(atom_name: str):
    a = str(atom_name).strip().upper()

    if not a:
        return "C"

    if a[0].isdigit() and len(a) > 1:
        a = a[1:]

    first = a[0]

    if first == "C":
        return "C"

    if first == "N":
        return "N"

    if first == "O":
        return "OA"

    if first == "S":
        return "S"

    if first == "H":
        return "HD"

    if first == "P":
        return "P"

    return "C"


def pdb_to_receptor_pdbqt(in_pdb: str, out_pdbqt: Path):
    p = Path(str(in_pdb))

    if not p.exists():
        return 0

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    count = 0

    with p.open("r", errors="ignore") as inp, out_pdbqt.open("w") as out:
        out.write("REMARK Vina-compatible WT receptor PDBQT generated from WT PDB\n")
        out.write("REMARK Charges set to 0.000; atom types assigned from protein atom names\n")

        for line in inp:
            if not line.startswith("ATOM"):
                continue

            try:
                atom_name = line[12:16]
                resname = line[17:20]
                chain = line[21:22] or "A"
                resnum = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue

            ad_type = atom_type_from_protein_atom(atom_name)
            charge = 0.000

            count += 1

            out.write(
                f"ATOM  "
                f"{count:5d} "
                f"{atom_name:<4s}"
                f" "
                f"{resname:>3s} "
                f"{chain:1s}"
                f"{resnum:>4s}"
                f"    "
                f"{x:8.3f}"
                f"{y:8.3f}"
                f"{z:8.3f}"
                f"{1.00:6.2f}"
                f"{0.00:6.2f}"
                f"    "
                f"{charge:6.3f} "
                f"{ad_type:>2s}\n"
            )

    if count == 0:
        try:
            out_pdbqt.unlink()
        except Exception:
            pass

    return count


def main():
    root = Path(".").resolve()

    mutant_base = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "mutant_docking_vina_compatible_pdbqt"
    )

    mutant_manifest_file = mutant_base / "pdb60_mutant_docking_jobs_manifest_VINA_COMPATIBLE.csv"

    reps_file = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "clean_unique_pdb60_representatives.csv"
    )

    if not mutant_manifest_file.exists():
        raise FileNotFoundError(f"Missing mutant manifest: {mutant_manifest_file}")

    if not reps_file.exists():
        raise FileNotFoundError(f"Missing representative PDB file: {reps_file}")

    outdir = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "wt_matching_docking_vina_compatible_pdbqt"
    )

    receptor_dir = outdir / "receptors_pdbqt"
    vina_out_dir = outdir / "vina_outputs"
    vina_log_dir = outdir / "vina_logs"

    for d in [outdir, receptor_dir, vina_out_dir, vina_log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    mutant_jobs = pd.read_csv(mutant_manifest_file, low_memory=False)
    reps = pd.read_csv(reps_file, low_memory=False)

    reps["structure_id"] = reps["structure_id"].astype(str).map(safe_id)
    mutant_jobs["structure_id"] = mutant_jobs["structure_id"].astype(str).map(safe_id)

    if "foldx_ready_pdb_path" not in reps.columns:
        raise RuntimeError("Representative file does not contain foldx_ready_pdb_path.")

    rep_lookup = (
        reps[
            [
                "structure_id",
                "foldx_ready_pdb_path",
                "uniprot_accession",
                "pdb_id",
                "organism_type",
                "gh_family",
            ]
        ]
        .drop_duplicates("structure_id")
        .copy()
    )

    jobs = mutant_jobs.merge(
        rep_lookup,
        on="structure_id",
        how="left",
        suffixes=("", "_rep"),
    )

    missing_wt = jobs[
        jobs["foldx_ready_pdb_path"].isna() |
        ~jobs["foldx_ready_pdb_path"].astype(str).map(lambda x: Path(x).exists())
    ].copy()

    missing_wt.to_csv(outdir / "missing_wt_receptor_paths.csv", index=False)

    jobs = jobs[
        jobs["foldx_ready_pdb_path"].notna() &
        jobs["foldx_ready_pdb_path"].astype(str).map(lambda x: Path(x).exists())
    ].copy()

    receptor_rows = []

    for sid, sub in jobs.groupby("structure_id"):
        first = sub.iloc[0]

        receptor_pdbqt = receptor_dir / f"{safe_id(sid)}_WT.pdbqt"

        atom_count = pdb_to_receptor_pdbqt(first["foldx_ready_pdb_path"], receptor_pdbqt)

        receptor_rows.append({
            "structure_id": sid,
            "wt_pdb_path": first["foldx_ready_pdb_path"],
            "wt_receptor_pdbqt": str(receptor_pdbqt),
            "pdbqt_atom_count": atom_count,
            "pdbqt_ok": atom_count > 20 and receptor_pdbqt.exists() and receptor_pdbqt.stat().st_size > 0,
            "pdbqt_generation_method": "vina_compatible_python_receptor_pdbqt_zero_charge",
        })

    receptor_map = pd.DataFrame(receptor_rows)

    receptor_map.to_csv(outdir / "pdb60_matching_wt_receptor_pdbqt_map.csv", index=False)

    # Rename WT receptor-preparation columns before merging to avoid collision
    # with existing mutant PDBQT columns from the mutant docking manifest.
    receptor_map_for_merge = receptor_map.rename(
        columns={
            "pdbqt_atom_count": "wt_pdbqt_atom_count",
            "pdbqt_ok": "wt_pdbqt_ok",
            "pdbqt_generation_method": "wt_pdbqt_generation_method",
        }
    )

    jobs = jobs.merge(
        receptor_map_for_merge[
            [
                "structure_id",
                "wt_pdb_path",
                "wt_receptor_pdbqt",
                "wt_pdbqt_atom_count",
                "wt_pdbqt_ok",
                "wt_pdbqt_generation_method",
            ]
        ],
        on="structure_id",
        how="left",
    )

    jobs = jobs[jobs["wt_pdbqt_ok"]].copy()

    wt_rows = []

    for _, r in jobs.iterrows():
        # We deliberately reuse the mutant job grid values so WT and mutant are paired.
        wt_job_id = safe_id(f"{r['job_id']}_WT")

        rr = {
            "wt_job_id": wt_job_id,
            "paired_mutant_job_id": r["job_id"],
            "structure_id": r["structure_id"],
            "uniprot_accession": r.get("uniprot_accession", r.get("uniprot_accession_rep", "")),
            "pdb_id": r.get("pdb_id", r.get("pdb_id_rep", "")),
            "organism_type": r.get("organism_type", r.get("organism_type_rep", "")),
            "gh_family": r.get("gh_family", r.get("gh_family_rep", "")),
            "foldx_mutation_code": r.get("foldx_mutation_code", ""),
            "mutant_ddg": r.get("ddg", np.nan),
            "ligand_name": r["ligand_name"],
            "ligand_pdbqt": r["ligand_pdbqt"],
            "wt_pdb_path": r["wt_pdb_path"],
            "wt_receptor_pdbqt": r["wt_receptor_pdbqt"],
            "center_x": r["center_x"],
            "center_y": r["center_y"],
            "center_z": r["center_z"],
            "size_x": r["size_x"],
            "size_y": r["size_y"],
            "size_z": r["size_z"],
            "vina_out_pdbqt": str(vina_out_dir / f"{wt_job_id}_{r['ligand_name']}_out.pdbqt"),
            "vina_log": str(vina_log_dir / f"{wt_job_id}_{r['ligand_name']}.log"),
            "pdbqt_generation_method": r["wt_pdbqt_generation_method"],
        }

        wt_rows.append(rr)

    wt_jobs = pd.DataFrame(wt_rows)

    wt_jobs.to_csv(outdir / "pdb60_matching_wt_docking_jobs_manifest.csv", index=False)

    vina_script = outdir / "01_run_pdb60_matching_wt_vina.sh"

    with vina_script.open("w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("set -uo pipefail\n\n")
        fh.write("echo 'Running matching WT docking for PDB60 mutant pairs...'\n\n")

        if wt_jobs.empty:
            fh.write("echo 'No valid WT jobs found.'\n")
        else:
            for _, r in wt_jobs.iterrows():
                if not r["ligand_pdbqt"] or not Path(str(r["ligand_pdbqt"])).exists():
                    continue

                fh.write(f"mkdir -p '{Path(r['vina_out_pdbqt']).parent}' '{Path(r['vina_log']).parent}'\n")
                fh.write(
                    f"if [ -s '{r['wt_receptor_pdbqt']}' ]; then\n"
                    f"  vina "
                    f"--receptor '{r['wt_receptor_pdbqt']}' "
                    f"--ligand '{r['ligand_pdbqt']}' "
                    f"--center_x {r['center_x']:.3f} "
                    f"--center_y {r['center_y']:.3f} "
                    f"--center_z {r['center_z']:.3f} "
                    f"--size_x {r['size_x']:.3f} "
                    f"--size_y {r['size_y']:.3f} "
                    f"--size_z {r['size_z']:.3f} "
                    f"--exhaustiveness 8 "
                    f"--num_modes 10 "
                    f"--out '{r['vina_out_pdbqt']}' "
                    f"> '{r['vina_log']}' 2>&1\n"
                    f"else\n"
                    f"  echo 'SKIPPED_MISSING_OR_EMPTY_WT_RECEPTOR_PDBQT: {r['wt_receptor_pdbqt']}' > '{r['vina_log']}'\n"
                    f"fi\n\n"
                )

    vina_script.chmod(0o755)

    report = outdir / "PDB60_MATCHING_WT_DOCKING_PREP_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 matching WT docking preparation report\n\n")
        fh.write(f"- Mutant docking jobs used as pairs: {len(mutant_jobs)}\n")
        fh.write(f"- WT receptor structures found: {receptor_map['structure_id'].nunique()}\n")
        fh.write(f"- WT receptor PDBQT files created: {(receptor_map['pdbqt_ok']).sum()}\n")
        fh.write(f"- Matching WT docking jobs prepared: {len(wt_jobs)}\n")
        fh.write(f"- Missing WT receptor path rows: {len(missing_wt)}\n")
        fh.write(f"- Vina script: `{vina_script}`\n\n")
        fh.write("## Pairing rule\n\n")
        fh.write(
            "WT docking jobs reuse the same ligand, grid center, and grid size as the corresponding mutant docking job. "
            "This enables paired mutant-WT docking comparison.\n"
        )

    print("\n[DONE] Matching WT docking jobs prepared.")
    print(f"Report: {report}")
    print(f"WT jobs manifest: {outdir / 'pdb60_matching_wt_docking_jobs_manifest.csv'}")
    print(f"Run Vina: {vina_script}")


if __name__ == "__main__":
    main()
