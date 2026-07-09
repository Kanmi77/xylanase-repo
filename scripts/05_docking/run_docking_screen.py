#!/usr/bin/env python3

from pathlib import Path
import argparse
import math
import subprocess

import pandas as pd
import yaml


TWO_LETTER_ELEMENTS = {
    "CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CO", "NI", "CU",
    "CD", "HG", "SR", "CS", "BA", "AL", "LI", "RB", "AG", "AU",
    "PB", "SN", "SE", "AS", "CR", "TI", "CA",
}

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL",
}


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def infer_element(atom_name, residue_name):
    atom = atom_name.strip()
    residue = residue_name.strip().upper()

    atom = atom.replace(" ", "")

    if not atom:
        return ""

    if atom[0].isdigit() and len(atom) > 1:
        atom = atom[1:]

    atom_upper = atom.upper()

    if residue in AA3:
        first = atom_upper[0]
        if first in {"C", "N", "O", "S", "H", "P"}:
            return first

    if len(atom_upper) >= 2:
        possible = atom_upper[:2]
        if possible in TWO_LETTER_ELEMENTS:
            if possible == "CA" and residue in AA3:
                return "C"
            return possible.title()

    first = atom_upper[0]

    if first in {"C", "N", "O", "S", "H", "P"}:
        return first

    return "C"


def fix_pdb_elements(input_pdb, output_pdb):
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL",
    }

    with open(input_pdb, "r", encoding="utf-8", errors="ignore") as src, open(
        output_pdb, "w", encoding="utf-8"
    ) as dst:
        for line in src:
            if line.startswith(("TER", "END")):
                dst.write(line)
                continue

            if not line.startswith("ATOM"):
                continue

            residue_name = line[17:20].strip().upper() if len(line) >= 20 else ""

            if residue_name not in standard_residues:
                continue

            atom_name = line[12:16] if len(line) >= 16 else ""
            element = line[76:78].strip() if len(line) >= 78 else ""

            if not element:
                element = infer_element(atom_name, residue_name)

            if element.upper() not in {"C", "N", "O", "S", "H", "P"}:
                element = infer_element(atom_name, residue_name)

            if element.upper() not in {"C", "N", "O", "S", "H", "P"}:
                continue

            line = line.rstrip("\n")
            if len(line) < 78:
                line = line.ljust(78)

            fixed = line[:76] + element.rjust(2) + line[78:]
            dst.write(fixed + "\n")


def parse_pdb_center(path):
    coords = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()

            if atom_name != "CA":
                continue

            try:
                coords.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
            except ValueError:
                continue

    if not coords:
        raise ValueError(f"No CA coordinates found in {path}")

    center = [
        sum(point[index] for point in coords) / len(coords)
        for index in range(3)
    ]

    max_distance = 0.0

    for point in coords:
        distance = math.sqrt(
            sum((point[index] - center[index]) ** 2 for index in range(3))
        )
        max_distance = max(max_distance, distance)

    box_size = max(22.0, min(44.0, max_distance * 1.8))

    return center, box_size


def prepare_receptor_with_meeko(meeko, pdb_path, receptor_path, log_path):
    receptor_path.parent.mkdir(parents=True, exist_ok=True)

    fixed_pdb = receptor_path.with_suffix(".fixed.pdb")
    fix_pdb_elements(pdb_path, fixed_pdb)

    output_base = receptor_path.with_suffix("")

    command = [
        meeko,
        "-i",
        str(fixed_pdb),
        "-o",
        str(output_base),
        "-p",
    ]

    with open(log_path, "w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(command) + "\n")
        result = subprocess.run(
            command,
            stdout=log,
            stderr=log,
            text=True,
        )

    expected_output = output_base.with_suffix(".pdbqt")

    if expected_output.exists() and expected_output != receptor_path:
        expected_output.rename(receptor_path)

    if not receptor_path.exists():
        return 1

    return result.returncode


def run_vina(vina, receptor, ligand, center, box_size, output, log_file, config):
    output.parent.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    command = [
        vina,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--center_x",
        str(round(center[0], 3)),
        "--center_y",
        str(round(center[1], 3)),
        "--center_z",
        str(round(center[2], 3)),
        "--size_x",
        str(round(box_size, 3)),
        "--size_y",
        str(round(box_size, 3)),
        "--size_z",
        str(round(box_size, 3)),
        "--exhaustiveness",
        str(config["exhaustiveness"]),
        "--num_modes",
        str(config["num_modes"]),
        "--energy_range",
        str(config["energy_range"]),
        "--out",
        str(output),
    ]

    with open(log_file, "w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(command) + "\n")
        result = subprocess.run(
            command,
            stdout=log,
            stderr=log,
            text=True,
        )

    return result.returncode


def parse_vina_score(log_file, output_file=None):
    if output_file is not None and output_file.exists():
        for line in output_file.read_text(errors="ignore").splitlines():
            line = line.strip()

            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()

                try:
                    return float(parts[3])
                except (ValueError, IndexError):
                    pass

    if not log_file.exists():
        return ""

    for line in log_file.read_text(errors="ignore").splitlines():
        line = line.strip()

        if not line:
            continue

        parts = line.split()

        if parts and parts[0].isdigit():
            try:
                return float(parts[1])
            except (ValueError, IndexError):
                continue

    return ""


def locate_repaired_pdb(row):
    work_dir = Path(str(row["work_dir"]))
    pdb_id = str(row["pdb_id"]).upper()

    candidates = sorted(work_dir.glob("*_Repair.pdb"))

    if candidates:
        return candidates[0]

    candidates = sorted(work_dir.glob(f"*{pdb_id}*.pdb"))

    if candidates:
        return candidates[0]

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Run AutoDock Vina docking for FoldX mutation records."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    docking_config = config["docking"]

    vina = docking_config["vina_binary"]
    meeko = docking_config["receptor_tool"]

    mutations = pd.read_csv(args.input, low_memory=False)
    mutations = mutations[mutations["status"] == "ready"].copy()

    if docking_config.get("test_mode", False):
        mutations = mutations.head(int(docking_config.get("max_mutations", 12)))

    receptor_dir = Path(docking_config["receptor_dir"])
    output_dir = Path(docking_config["output_dir"])

    rows = []

    for _, row in mutations.iterrows():
        accession = str(row["uniprot_accession"])
        pdb_id = str(row["pdb_id"]).upper()
        mutation_code = str(row["foldx_mutation_code"])
        receptor_tag = f"{accession}_{pdb_id}_{mutation_code}"

        pdb_path = locate_repaired_pdb(row)

        if pdb_path is None:
            rows.append(
                {
                    "uniprot_accession": accession,
                    "organism_type": row.get("organism_type", ""),
                    "gh_family": row.get("gh_family", ""),
                    "pdb_id": pdb_id,
                    "foldx_mutation_code": mutation_code,
                    "ligand": "",
                    "vina_affinity": "",
                    "status": "failed_no_receptor_pdb",
                }
            )
            continue

        center, auto_box_size = parse_pdb_center(pdb_path)
        box_size = float(docking_config.get("box_size", auto_box_size))

        receptor_pdbqt = receptor_dir / f"{receptor_tag}.pdbqt"
        receptor_log = receptor_dir / f"{receptor_tag}.prepare.log"

        prep_code = prepare_receptor_with_meeko(
            meeko,
            pdb_path,
            receptor_pdbqt,
            receptor_log,
        )

        for ligand_name, ligand_path_value in docking_config["ligands"].items():
            ligand_path = Path(str(ligand_path_value))

            vina_out = output_dir / f"{receptor_tag}_{ligand_name}.out.pdbqt"
            vina_log = output_dir / f"{receptor_tag}_{ligand_name}.log"

            vina_code = 1
            affinity = ""

            if prep_code == 0 and receptor_pdbqt.exists() and ligand_path.exists():
                vina_code = run_vina(
                    vina,
                    receptor_pdbqt,
                    ligand_path,
                    center,
                    box_size,
                    vina_out,
                    vina_log,
                    docking_config,
                )
                affinity = parse_vina_score(vina_log, vina_out)

            status = "ready" if affinity != "" else "failed"

            rows.append(
                {
                    "uniprot_accession": accession,
                    "organism_type": row.get("organism_type", ""),
                    "gh_family": row.get("gh_family", ""),
                    "pdb_id": pdb_id,
                    "foldx_mutation_code": mutation_code,
                    "foldx_ddg": row.get("foldx_ddg", ""),
                    "ligand": ligand_name,
                    "receptor_pdb": str(pdb_path),
                    "receptor_pdbqt": str(receptor_pdbqt),
                    "center_x": round(center[0], 3),
                    "center_y": round(center[1], 3),
                    "center_z": round(center[2], 3),
                    "box_size": round(box_size, 3),
                    "vina_affinity": affinity,
                    "prepare_exit_code": prep_code,
                    "vina_exit_code": vina_code,
                    "vina_output": str(vina_out),
                    "vina_log": str(vina_log),
                    "status": status,
                }
            )

    result = pd.DataFrame(rows)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    print(f"Docking rows: {len(result)}")

    if len(result):
        print(f"Ready rows: {(result['status'] == 'ready').sum()}")

    print(f"Saved: {output_path}")

    if result.empty or not (result["status"] == "ready").any():
        raise SystemExit("No docking result was successfully produced.")


if __name__ == "__main__":
    main()
