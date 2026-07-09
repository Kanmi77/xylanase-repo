#!/usr/bin/env python3
# Purpose: Run Vina docking.

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
    """
    Protein-safe AutoDock/Vina atom typing.

    Important:
    CA in protein PDB = alpha carbon, not calcium.
    CB, CG, CD, CE, CZ etc. are carbon atoms.
    """

    a = str(atom_name).strip().upper()

    if not a:
        return "C"

    # Remove leading digit, e.g. 1H -> H
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

    # Conservative fallback for protein atoms.
    return "C"


def pdb_to_receptor_pdbqt(in_pdb: str, out_pdbqt: Path):
    p = Path(str(in_pdb))

    if not p.exists():
        return 0

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    count = 0

    with p.open("r", errors="ignore") as inp, out_pdbqt.open("w") as out:
        out.write("REMARK Vina-compatible rigid receptor PDBQT generated from cleaned FoldX mutant PDB\n")
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
        / "mutant_docking_cleaned"
        / "pdb60_mutant_docking_jobs_manifest_CLEANED.csv"
    )

    if not in_manifest.exists():
        raise FileNotFoundError(f"Missing input manifest: {in_manifest}")

    outdir = (
        root
        / "results"
        / "optionC_original_only"
        / "pdb60_foldx_mutations"
        / "parsed_ddg_and_docking"
        / "mutant_docking_vina_compatible_pdbqt"
    )

    receptor_dir = outdir / "receptors_pdbqt"
    vina_out_dir = outdir / "vina_outputs"
    vina_log_dir = outdir / "vina_logs"

    for d in [outdir, receptor_dir, vina_out_dir, vina_log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    jobs = pd.read_csv(in_manifest, low_memory=False)

    if jobs.empty:
        raise RuntimeError("Input manifest is empty.")

    receptor_rows = []

    for job_id, sub in jobs.groupby("job_id"):
        first = sub.iloc[0]

        receptor_pdbqt = receptor_dir / f"{safe_id(job_id)}.pdbqt"
        atom_count = pdb_to_receptor_pdbqt(first["mutant_pdb_path"], receptor_pdbqt)

        receptor_rows.append({
            "job_id": job_id,
            "structure_id": first.get("structure_id", ""),
            "foldx_mutation_code": first.get("foldx_mutation_code", ""),
            "ddg": first.get("ddg", np.nan),
            "mutant_pdb_path": first.get("mutant_pdb_path", ""),
            "receptor_pdbqt": str(receptor_pdbqt),
            "pdbqt_atom_count": atom_count,
            "pdbqt_ok": atom_count > 20 and receptor_pdbqt.exists() and receptor_pdbqt.stat().st_size > 0,
            "pdbqt_generation_method": "vina_compatible_python_receptor_pdbqt_zero_charge",
        })

    receptor_map = pd.DataFrame(receptor_rows)

    receptor_map.to_csv(outdir / "pdb60_vina_compatible_receptor_pdbqt_map.csv", index=False)

    jobs = jobs.drop(columns=[c for c in ["receptor_pdbqt"] if c in jobs.columns])

    jobs = jobs.merge(
        receptor_map[
            [
                "job_id",
                "receptor_pdbqt",
                "pdbqt_atom_count",
                "pdbqt_ok",
                "pdbqt_generation_method",
            ]
        ],
        on="job_id",
        how="left",
    )

    jobs = jobs[jobs["pdbqt_ok"]].copy()

    fixed_rows = []

    for _, r in jobs.iterrows():
        cx, cy, cz, sx, sy, sz = pdb_center_and_box(r["mutant_pdb_path"])

        rr = r.to_dict()
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

    fixed_jobs.to_csv(outdir / "pdb60_mutant_docking_jobs_manifest_VINA_COMPATIBLE.csv", index=False)

    vina_script = outdir / "01_run_pdb60_mutant_vina_VINA_COMPATIBLE.sh"

    with vina_script.open("w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write("set -uo pipefail\n\n")
        fh.write("echo 'Running PDB60 mutant docking using Vina-compatible receptor PDBQT files...'\n\n")

        if fixed_jobs.empty:
            fh.write("echo 'No valid Vina-compatible PDBQT jobs found.'\n")
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

    report = outdir / "PDB60_VINA_COMPATIBLE_PDBQT_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 Vina-compatible receptor PDBQT report\n\n")
        fh.write(f"- Input docking jobs: {len(pd.read_csv(in_manifest, low_memory=False))}\n")
        fh.write(f"- Unique mutant receptors: {receptor_map['job_id'].nunique()}\n")
        fh.write(f"- Vina-compatible receptor PDBQT files created: {(receptor_map['pdbqt_ok']).sum()}\n")
        fh.write(f"- Docking jobs retained: {len(fixed_jobs)}\n")
        fh.write(f"- Vina script: `{vina_script}`\n\n")
        fh.write("## Fix applied\n\n")
        fh.write(
            "Protein alpha-carbon atoms named CA were assigned atom type C instead of calcium. "
            "Ligand-style ROOT/ENDROOT/TORSDOF records were removed from receptor PDBQT files.\n"
        )
        fh.write("\n## Limitation\n\n")
        fh.write(
            "Charges were set to 0.000 because this is a fallback receptor preparation method. "
            "Docking scores from this correction layer should be treated as comparative screening values.\n"
        )

    print("\n[DONE] Vina-compatible receptor PDBQT files and Vina script prepared.")
    print(f"Report: {report}")
    print(f"Run Vina: {vina_script}")


if __name__ == "__main__":
    main()
