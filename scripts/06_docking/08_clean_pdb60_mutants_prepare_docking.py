#!/usr/bin/env python3
"""
Clean FoldX mutant PDB files into strict PDB format and regenerate docking scripts.

Why this is needed:
- FoldX mutant PDB files can contain ATOM records but still fail Open Babel parsing.
- This script rewrites ATOM/HETATM lines into strict PDB column format.
- It then prepares Open Babel and Vina scripts using the cleaned PDB files.

Inputs:
    results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/mutant_docking_fixed/pdb60_mutant_docking_jobs_manifest_fixed.csv

Outputs:
    results/optionC_original_only/pdb60_foldx_mutations/parsed_ddg_and_docking/mutant_docking_cleaned/
"""

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


def guess_element(atom_name: str):
    a = atom_name.strip()

    if not a:
        return ""

    # Common protein atom names.
    if a[0].isdigit() and len(a) > 1:
        return a[1].upper()

    if len(a) >= 2 and a[:2].strip().upper() in {"FE", "ZN", "MG", "MN", "CA", "NA", "CL"}:
        return a[:2].strip().upper()

    return a[0].upper()


def parse_atom_line(line):
    """
    Parse ATOM/HETATM line using PDB columns first.
    Falls back to split parsing if needed.
    """

    rec = line[0:6].strip()

    try:
        serial = int(line[6:11])
    except Exception:
        serial = None

    atom_name = line[12:16].strip() if len(line) >= 16 else ""
    altloc = line[16:17].strip() if len(line) >= 17 else ""
    resname = line[17:20].strip() if len(line) >= 20 else "UNK"
    chain = line[21:22].strip() if len(line) >= 22 else "A"

    try:
        resnum = int(line[22:26])
    except Exception:
        resnum = None

    icode = line[26:27].strip() if len(line) >= 27 else ""

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        parts = line.split()

        # Very conservative fallback.
        # Typical split format:
        # ATOM serial atom res chain resnum x y z ...
        try:
            if len(parts) >= 9:
                if serial is None:
                    serial = int(parts[1])

                if not atom_name:
                    atom_name = parts[2]

                if not resname:
                    resname = parts[3]

                if not chain:
                    chain = parts[4]

                if resnum is None:
                    resnum = int(float(parts[5]))

                x = float(parts[6])
                y = float(parts[7])
                z = float(parts[8])
            else:
                return None
        except Exception:
            return None

    try:
        occ = float(line[54:60])
    except Exception:
        occ = 1.00

    try:
        bfac = float(line[60:66])
    except Exception:
        bfac = 0.00

    element = line[76:78].strip() if len(line) >= 78 else ""

    if not element:
        element = guess_element(atom_name)

    if not atom_name or not resname or resnum is None:
        return None

    if not chain:
        chain = "A"

    if rec not in {"ATOM", "HETATM"}:
        rec = "ATOM"

    return {
        "record": rec,
        "serial": serial,
        "atom_name": atom_name,
        "altloc": altloc[:1],
        "resname": resname[:3].upper(),
        "chain": chain[:1],
        "resnum": int(resnum),
        "icode": icode[:1],
        "x": x,
        "y": y,
        "z": z,
        "occ": occ,
        "bfac": bfac,
        "element": element[:2].upper(),
    }


def format_pdb_atom(atom, serial):
    atom_name = atom["atom_name"]

    # PDB atom name alignment.
    if len(atom_name) < 4 and atom["element"] and len(atom["element"]) == 1:
        atom_field = f" {atom_name:<3s}"
    else:
        atom_field = f"{atom_name:>4s}"

    return (
        f"{atom['record']:<6s}"
        f"{serial:5d} "
        f"{atom_field}"
        f"{atom['altloc']:<1s}"
        f"{atom['resname']:>3s} "
        f"{atom['chain']:<1s}"
        f"{atom['resnum']:4d}"
        f"{atom['icode']:<1s}"
        f"   "
        f"{atom['x']:8.3f}"
        f"{atom['y']:8.3f}"
        f"{atom['z']:8.3f}"
        f"{atom['occ']:6.2f}"
        f"{atom['bfac']:6.2f}"
        f"          "
        f"{atom['element']:>2s}"
    )


def clean_pdb(in_pdb: str, out_pdb: Path):
    atoms = []

    p = Path(str(in_pdb))

    if not p.exists():
        return 0

    with p.open("r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = parse_atom_line(line)

                if atom is not None:
                    atoms.append(atom)

    if not atoms:
        return 0

    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    with out_pdb.open("w") as out:
        for i, atom in enumerate(atoms, start=1):
            out.write(format_pdb_atom(atom, i) + "\n")

        out.write("TER\n")
        out.write("END\n")

    return len(atoms)


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


def main():
    root = Path(".").resolve()

    in_manifest = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "mutant_docking_fixed"
        / "pdb60_mutant_docking_jobs_manifest_fixed.csv"
    )

    if not in_manifest.exists():
        raise FileNotFoundError(f"Missing manifest: {in_manifest}")

    outdir = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "mutant_docking_cleaned"
    )

    clean_pdb_dir = outdir / "cleaned_mutant_pdb"
    receptor_dir = outdir / "receptors_pdbqt"
    vina_out_dir = outdir / "vina_outputs"
    vina_log_dir = outdir / "vina_logs"

    for d in [outdir, clean_pdb_dir, receptor_dir, vina_out_dir, vina_log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    jobs = pd.read_csv(in_manifest, low_memory=False)

    if jobs.empty:
        raise RuntimeError("Input manifest is empty.")

    unique_receptors = (
        jobs[
            [
                "job_id",
                "structure_id",
                "foldx_mutation_code",
                "ddg",
                "mutant_pdb_path",
            ]
        ]
        .drop_duplicates()
        .copy()
    )

    clean_map_rows = []

    for _, r in unique_receptors.iterrows():
        job_id = safe_id(r["job_id"])
        cleaned = clean_pdb_dir / f"{job_id}.pdb"

        atom_count = clean_pdb(r["mutant_pdb_path"], cleaned)

        clean_map_rows.append({
            "job_id": job_id,
            "structure_id": r.get("structure_id", ""),
            "foldx_mutation_code": r.get("foldx_mutation_code", ""),
            "ddg": r.get("ddg", np.nan),
            "original_mutant_pdb_path": r.get("mutant_pdb_path", ""),
            "cleaned_mutant_pdb_path": str(cleaned),
            "cleaned_atom_count": atom_count,
            "cleaned_pdb_ok": atom_count > 20,
        })

    clean_map = pd.DataFrame(clean_map_rows)

    clean_map.to_csv(outdir / "pdb60_cleaned_mutant_pdb_map.csv", index=False)

    jobs = jobs.merge(
        clean_map[
            [
                "job_id",
                "cleaned_mutant_pdb_path",
                "cleaned_atom_count",
                "cleaned_pdb_ok",
            ]
        ],
        on="job_id",
        how="left",
    )

    jobs = jobs[jobs["cleaned_pdb_ok"]].copy()

    fixed_rows = []

    for _, r in jobs.iterrows():
        receptor_pdbqt = receptor_dir / f"{safe_id(r['job_id'])}.pdbqt"

        cx, cy, cz, sx, sy, sz = pdb_center_and_box(r["cleaned_mutant_pdb_path"])

        rr = r.to_dict()
        rr["mutant_pdb_path_original"] = rr.get("mutant_pdb_path", "")
        rr["mutant_pdb_path"] = rr["cleaned_mutant_pdb_path"]
        rr["receptor_pdbqt"] = str(receptor_pdbqt)
        rr["center_x"] = cx
        rr["center_y"] = cy
        rr["center_z"] = cz
        rr["size_x"] = sx
        rr["size_y"] = sy
        rr["size_z"] = sz
        rr["vina_out_pdbqt"] = str(vina_out_dir / f"{safe_id(r['job_id'])}_{r['ligand_name']}_out.pdbqt")
        rr["vina_log"] = str(vina_log_dir / f"{safe_id(r['job_id'])}_{r['ligand_name']}.log")
        fixed_rows.append(rr)

    fixed_jobs = pd.DataFrame(fixed_rows)

    fixed_jobs.to_csv(outdir / "pdb60_mutant_docking_jobs_manifest_CLEANED.csv", index=False)

    prep_script = outdir / "01_prepare_pdb60_mutant_receptors_obabel_CLEANED.sh"

    with prep_script.open("w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("set -uo pipefail\n\n")
        fh.write("echo 'Preparing cleaned PDB60 mutant receptors using Open Babel...'\n\n")

        if fixed_jobs.empty:
            fh.write("echo 'No valid cleaned jobs found.'\n")
        else:
            for mutant_pdb, receptor_pdbqt in (
                fixed_jobs[["mutant_pdb_path", "receptor_pdbqt"]]
                .drop_duplicates()
                .itertuples(index=False)
            ):
                fh.write(f"mkdir -p '{Path(receptor_pdbqt).parent}'\n")
                fh.write(f"echo 'Preparing cleaned receptor: {Path(mutant_pdb).name}'\n")
                fh.write(
                    f"if [ ! -s '{receptor_pdbqt}' ]; then\n"
                    f"  obabel -ipdb '{mutant_pdb}' -opdbqt -O '{receptor_pdbqt}' --partialcharge gasteiger\n"
                    f"  if [ $? -ne 0 ]; then echo 'FAILED_OBABEL: {mutant_pdb}'; fi\n"
                    f"fi\n\n"
                )

    prep_script.chmod(0o755)

    vina_script = outdir / "02_run_pdb60_mutant_vina_CLEANED.sh"

    with vina_script.open("w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("set -uo pipefail\n\n")
        fh.write("echo 'Running cleaned PDB60 mutant docking with AutoDock Vina...'\n\n")

        if fixed_jobs.empty:
            fh.write("echo 'No valid cleaned jobs found.'\n")
        else:
            for _, r in fixed_jobs.iterrows():
                if not r["ligand_pdbqt"] or not Path(str(r["ligand_pdbqt"])).exists():
                    continue

                fh.write(f"mkdir -p '{Path(r['vina_out_pdbqt']).parent}' '{Path(r['vina_log']).parent}'\n")
                fh.write(
                    f"if [ -s '{r['receptor_pdbqt']}' ]; then\n"
                    f"  vina "
                    f"--receptor '{r['receptor_pdbqt']}' "
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
                    f"  echo 'SKIPPED_MISSING_OR_EMPTY_RECEPTOR_PDBQT: {r['receptor_pdbqt']}' > '{r['vina_log']}'\n"
                    f"fi\n\n"
                )

    vina_script.chmod(0o755)

    report = outdir / "PDB60_CLEANED_MUTANT_DOCKING_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 cleaned mutant docking preparation report\n\n")
        fh.write(f"- Input docking jobs: {len(pd.read_csv(in_manifest, low_memory=False))}\n")
        fh.write(f"- Unique mutant receptors: {len(unique_receptors)}\n")
        fh.write(f"- Cleaned mutant PDB files created: {(clean_map['cleaned_pdb_ok']).sum()}\n")
        fh.write(f"- Cleaned docking jobs retained: {len(fixed_jobs)}\n")
        fh.write(f"- Receptor-prep script: `{prep_script}`\n")
        fh.write(f"- Vina script: `{vina_script}`\n\n")
        fh.write("## Notes\n\n")
        fh.write(
            "FoldX mutant PDB files were rewritten into stricter PDB column format before Open Babel conversion. "
            "This fixes cases where FoldX files contain ATOM records but Open Babel rejects the file formatting.\n"
        )

    print("\n[DONE] Cleaned mutant PDB docking scripts prepared.")
    print(f"Report: {report}")
    print(f"Prepare receptors: {prep_script}")
    print(f"Run Vina: {vina_script}")


if __name__ == "__main__":
    main()
