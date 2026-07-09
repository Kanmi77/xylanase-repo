#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess

import pandas as pd
import yaml


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

MUTATION_TARGETS = ["A", "G", "V", "S", "T", "L", "N", "K", "P", "F"]


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def find_foldx(binary_name):
    direct = shutil.which(binary_name)

    if direct:
        return direct

    for candidate in ["FoldX", "foldx", "foldx_20231231"]:
        found = shutil.which(candidate)
        if found:
            return found

    raise SystemExit(
        "FoldX binary was not found. Set foldx.foldx_binary in config/workflow_config.yaml."
    )


def parse_residues_from_pdb(path):
    residues = []
    seen = set()

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            alt_loc = line[16].strip()

            if atom_name != "CA":
                continue

            if alt_loc not in {"", "A"}:
                continue

            residue_name = line[17:20].strip()
            chain_id = line[21].strip() or "A"
            residue_number = line[22:26].strip()

            if residue_name not in THREE_TO_ONE:
                continue

            key = (chain_id, residue_number)

            if key in seen:
                continue

            seen.add(key)

            residues.append(
                {
                    "wild_type": THREE_TO_ONE[residue_name],
                    "chain": chain_id,
                    "residue_number": residue_number,
                }
            )

    return residues


def choose_mutant(wild_type, index):
    for shift in range(len(MUTATION_TARGETS)):
        candidate = MUTATION_TARGETS[(index + shift) % len(MUTATION_TARGETS)]
        if candidate != wild_type:
            return candidate
    return "A" if wild_type != "A" else "G"


def select_mutations(residues, count):
    if not residues:
        return []

    if len(residues) <= count:
        selected = residues
    else:
        selected = []
        for index in range(count):
            position = round(index * (len(residues) - 1) / max(count - 1, 1))
            selected.append(residues[position])

    mutations = []

    for index, residue in enumerate(selected):
        mutant = choose_mutant(residue["wild_type"], index)
        mutation_code = (
            f"{residue['wild_type']}"
            f"{residue['chain']}"
            f"{residue['residue_number']}"
            f"{mutant}"
        )

        mutations.append(
            {
                "foldx_mutation_code": mutation_code,
                "wild_type": residue["wild_type"],
                "chain": residue["chain"],
                "residue_number": residue["residue_number"],
                "mutant": mutant,
            }
        )

    return mutations


def run_command(command, cwd, log_path):
    with open(log_path, "a", encoding="utf-8") as log:
        log.write("\nCOMMAND: " + " ".join(command) + "\n")
        result = subprocess.run(
            command,
            cwd=cwd,
            stdout=log,
            stderr=log,
            text=True,
        )

    return result.returncode


def parse_buildmodel_dif(work_dir):
    dif_files = sorted(work_dir.glob("Dif_*.fxout"))

    values = []

    for path in dif_files:
        lines = [
            line.strip()
            for line in path.read_text(errors="ignore").splitlines()
            if line.strip() and not line.startswith("#")
        ]

        for line in lines:
            parts = line.replace("\t", " ").split()
            numeric_values = []

            for part in parts:
                try:
                    numeric_values.append(float(part))
                except ValueError:
                    pass

            if numeric_values:
                values.append(numeric_values[0])

    return values


def main():
    parser = argparse.ArgumentParser(
        description="Run FoldX BuildModel mutation screening."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    foldx_config = config["foldx"]
    foldx_binary = find_foldx(foldx_config["foldx_binary"])

    wt = pd.read_csv(args.input, low_memory=False)
    wt = wt[wt["status"] == "ready"].copy()

    mutation_count = int(foldx_config["mutations_per_structure"])
    work_root = Path(foldx_config["work_dir"]) / "mutation_screen"
    work_root.mkdir(parents=True, exist_ok=True)

    rows = []

    for _, row in wt.iterrows():
        repaired_pdb = Path(row["repaired_pdb"])

        if not repaired_pdb.exists():
            continue

        accession = str(row["uniprot_accession"])
        pdb_id = str(row["pdb_id"]).upper()
        tag = f"{accession}_{pdb_id}"

        work_dir = work_root / tag
        work_dir.mkdir(parents=True, exist_ok=True)

        pdb_copy = work_dir / repaired_pdb.name
        shutil.copy2(repaired_pdb, pdb_copy)

        rotabase = Path(foldx_config.get("rotabase", ""))
        if rotabase.exists():
            shutil.copy2(rotabase, work_dir / "rotabase.txt")

        residues = parse_residues_from_pdb(pdb_copy)
        mutations = select_mutations(residues, mutation_count)

        if not mutations:
            continue

        mutation_file = work_dir / "individual_list.txt"

        with open(mutation_file, "w", encoding="utf-8") as handle:
            for mutation in mutations:
                handle.write(mutation["foldx_mutation_code"] + ";\n")

        log_path = work_dir / "foldx_buildmodel.log"

        exit_code = run_command(
            [
                foldx_binary,
                "--command=BuildModel",
                f"--pdb={pdb_copy.name}",
                f"--mutant-file={mutation_file.name}",
                "--numberOfRuns=1",
            ],
            cwd=work_dir,
            log_path=log_path,
        )

        ddg_values = parse_buildmodel_dif(work_dir)

        for index, mutation in enumerate(mutations):
            ddg = ddg_values[index] if index < len(ddg_values) else ""

            rows.append(
                {
                    "uniprot_accession": accession,
                    "organism_type": row.get("organism_type", ""),
                    "gh_family": row.get("gh_family", ""),
                    "pdb_id": pdb_id,
                    "foldx_mutation_code": mutation["foldx_mutation_code"],
                    "wild_type": mutation["wild_type"],
                    "chain": mutation["chain"],
                    "residue_number": mutation["residue_number"],
                    "mutant": mutation["mutant"],
                    "foldx_ddg": ddg,
                    "buildmodel_exit_code": exit_code,
                    "status": "ready" if ddg != "" else "failed",
                    "work_dir": str(work_dir),
                }
            )

    result = pd.DataFrame(rows)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    print(f"FoldX mutation rows: {len(result)}")

    if len(result):
        print(f"Successful mutation rows: {(result['status'] == 'ready').sum()}")

    print(f"Saved: {output_path}")

    if result.empty or not (result["status"] == "ready").any():
        raise SystemExit("No FoldX mutation result was successfully parsed.")


if __name__ == "__main__":
    main()
