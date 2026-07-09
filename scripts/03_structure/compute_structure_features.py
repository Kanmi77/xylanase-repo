#!/usr/bin/env python3

from pathlib import Path
from urllib.request import Request, urlopen
import argparse
import math

import pandas as pd


STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def download_pdb(pdb_id, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and output_path.stat().st_size > 0:
        return True

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    request = Request(url, headers={"User-Agent": "enzyme-thermostability-workflow"})

    try:
        with urlopen(request, timeout=120) as response:
            text = response.read().decode("utf-8")
    except Exception as error:
        print(f"Download failed for {pdb_id}: {error}")
        return False

    if not text.startswith(("HEADER", "TITLE", "ATOM")):
        print(f"Unexpected PDB response for {pdb_id}")
        return False

    output_path.write_text(text, encoding="utf-8")
    return True


def parse_atom_line(line):
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip()
    chain_id = line[21].strip() or "_"
    residue_number = line[22:26].strip()
    insertion_code = line[26].strip()
    alternate_location = line[16].strip()

    if alternate_location not in {"", "A"}:
        return None

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    return {
        "atom_name": atom_name,
        "residue_name": residue_name,
        "chain_id": chain_id,
        "residue_id": f"{chain_id}:{residue_number}:{insertion_code}",
        "atom_key": f"{chain_id}:{residue_number}:{insertion_code}:{atom_name}",
        "x": x,
        "y": y,
        "z": z,
    }


def read_first_model_atoms(path):
    atoms = []
    seen_atom_keys = set()
    in_first_model = False
    model_seen = False

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if model_seen:
                    break
                model_seen = True
                in_first_model = True
                continue

            if line.startswith("ENDMDL") and in_first_model:
                break

            if model_seen and not in_first_model:
                continue

            if not line.startswith("ATOM"):
                continue

            atom = parse_atom_line(line)

            if atom is None:
                continue

            if atom["atom_key"] in seen_atom_keys:
                continue

            seen_atom_keys.add(atom["atom_key"])
            atoms.append(atom)

    return atoms


def compute_radius_of_gyration(coords):
    if not coords:
        return ""

    center = [
        sum(point[index] for point in coords) / len(coords)
        for index in range(3)
    ]

    squared_distances = []

    for point in coords:
        squared_distances.append(
            sum((point[index] - center[index]) ** 2 for index in range(3))
        )

    return math.sqrt(sum(squared_distances) / len(squared_distances))


def summarise_pdb(path):
    atoms = read_first_model_atoms(path)

    if not atoms:
        return {
            "structure_loaded": False,
            "atom_count": 0,
            "residue_count": 0,
            "chain_count": 0,
            "ca_atom_count": 0,
            "radius_of_gyration_ca": "",
            "positive_residue_count": "",
            "negative_residue_count": "",
            "polar_residue_count": "",
            "positive_fraction": "",
            "negative_fraction": "",
            "polar_fraction": "",
        }

    residues = {
        atom["residue_id"]
        for atom in atoms
        if atom["residue_name"] in STANDARD_AA
    }

    chains = {atom["chain_id"] for atom in atoms}

    ca_atoms = [
        atom for atom in atoms
        if atom["atom_name"] == "CA" and atom["residue_name"] in STANDARD_AA
    ]

    coords = [(atom["x"], atom["y"], atom["z"]) for atom in ca_atoms]

    positive_residues = {
        atom["residue_id"]
        for atom in atoms
        if atom["residue_name"] in {"LYS", "ARG", "HIS"}
    }

    negative_residues = {
        atom["residue_id"]
        for atom in atoms
        if atom["residue_name"] in {"ASP", "GLU"}
    }

    polar_residues = {
        atom["residue_id"]
        for atom in atoms
        if atom["residue_name"] in {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
    }

    residue_count = len(residues)

    return {
        "structure_loaded": True,
        "atom_count": len(atoms),
        "residue_count": residue_count,
        "chain_count": len(chains),
        "ca_atom_count": len(ca_atoms),
        "radius_of_gyration_ca": compute_radius_of_gyration(coords),
        "positive_residue_count": len(positive_residues),
        "negative_residue_count": len(negative_residues),
        "polar_residue_count": len(polar_residues),
        "positive_fraction": len(positive_residues) / residue_count if residue_count else "",
        "negative_fraction": len(negative_residues) / residue_count if residue_count else "",
        "polar_fraction": len(polar_residues) / residue_count if residue_count else "",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Download PDB files and compute structural features."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    inventory = pd.read_csv(args.input, low_memory=False)
    rows = []

    for _, row in inventory.iterrows():
        pdb_id = str(row["pdb_id"]).upper()
        structure_file = Path(row["structure_file"])

        downloaded = download_pdb(pdb_id, structure_file)

        if downloaded:
            features = summarise_pdb(structure_file)
        else:
            features = {
                "structure_loaded": False,
                "atom_count": 0,
                "residue_count": 0,
                "chain_count": 0,
                "ca_atom_count": 0,
                "radius_of_gyration_ca": "",
                "positive_residue_count": "",
                "negative_residue_count": "",
                "polar_residue_count": "",
                "positive_fraction": "",
                "negative_fraction": "",
                "polar_fraction": "",
            }

        result = row.to_dict()
        result.update(features)
        rows.append(result)

    output_table = pd.DataFrame(rows)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output_table.to_csv(output, index=False)

    print(f"Structure-feature rows: {len(output_table)}")

    if len(output_table):
        print(
            "Loaded structures:",
            int(output_table["structure_loaded"].astype(bool).sum()),
        )

    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
