#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess

import pandas as pd
import yaml


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


def copy_input_structure(row, work_dir):
    source = Path(row["structure_file"])
    pdb_id = str(row["pdb_id"]).upper()
    accession = str(row["uniprot_accession"])
    target = work_dir / f"{accession}_{pdb_id}.pdb"

    if not source.exists():
        raise FileNotFoundError(f"Missing structure file: {source}")

    shutil.copy2(source, target)
    return target


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


def parse_stability_output(work_dir):
    candidates = sorted(work_dir.glob("*_ST.fxout")) + sorted(work_dir.glob("*.fxout"))

    for path in candidates:
        try:
            lines = [
                line.strip()
                for line in path.read_text(errors="ignore").splitlines()
                if line.strip()
            ]
        except Exception:
            continue

        for line in reversed(lines):
            parts = line.replace("\t", " ").split()
            numeric_values = []

            for part in parts:
                try:
                    numeric_values.append(float(part))
                except ValueError:
                    pass

            if numeric_values:
                return numeric_values[0], path.name

    return None, ""


def main():
    parser = argparse.ArgumentParser(
        description="Run FoldX RepairPDB and Stability on available structures."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    foldx_config = config["foldx"]
    foldx_binary = find_foldx(foldx_config["foldx_binary"])

    work_root = Path(foldx_config["work_dir"]) / "wt_stability"
    repaired_dir = Path(foldx_config["repaired_dir"])
    output_path = Path(args.output)

    work_root.mkdir(parents=True, exist_ok=True)
    repaired_dir.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    structures = pd.read_csv(args.input, low_memory=False)
    structures = structures[structures["structure_loaded"].astype(bool)].copy()

    if foldx_config.get("test_mode", False):
        per_group = int(foldx_config.get("max_structures_per_group", 1))
        structures = (
            structures
            .sort_values(["organism_type", "gh_family", "uniprot_accession", "pdb_id"])
            .groupby(["organism_type", "gh_family"], as_index=False, group_keys=False)
            .head(per_group)
        )
    else:
        max_structures = foldx_config.get("max_structures")
        if max_structures:
            structures = structures.head(int(max_structures))

    rows = []

    for _, row in structures.iterrows():
        pdb_id = str(row["pdb_id"]).upper()
        accession = str(row["uniprot_accession"])
        tag = f"{accession}_{pdb_id}"
        work_dir = work_root / tag
        work_dir.mkdir(parents=True, exist_ok=True)

        log_path = work_dir / "foldx_wt.log"
        input_pdb = copy_input_structure(row, work_dir)

        rotabase = Path(foldx_config.get("rotabase", ""))
        if rotabase.exists():
            shutil.copy2(rotabase, work_dir / "rotabase.txt")

        repair_code = run_command(
            [
                foldx_binary,
                "--command=RepairPDB",
                f"--pdb={input_pdb.name}",
            ],
            cwd=work_dir,
            log_path=log_path,
        )

        repaired_name = input_pdb.name.replace(".pdb", "_Repair.pdb")
        repaired_pdb = work_dir / repaired_name

        if repaired_pdb.exists():
            shutil.copy2(repaired_pdb, repaired_dir / repaired_name)

        stability_code = 1
        foldx_total_energy = None
        parsed_file = ""

        if repaired_pdb.exists():
            stability_code = run_command(
                [
                    foldx_binary,
                    "--command=Stability",
                    f"--pdb={repaired_pdb.name}",
                ],
                cwd=work_dir,
                log_path=log_path,
            )
            foldx_total_energy, parsed_file = parse_stability_output(work_dir)

        residue_count = pd.to_numeric(row.get("residue_count", ""), errors="coerce")

        if foldx_total_energy is not None and residue_count and residue_count > 0:
            foldx_energy_per_residue = foldx_total_energy / residue_count
        else:
            foldx_energy_per_residue = ""

        rows.append(
            {
                "uniprot_accession": accession,
                "organism_type": row.get("organism_type", ""),
                "gh_family": row.get("gh_family", ""),
                "pdb_id": pdb_id,
                "structure_file": str(input_pdb),
                "repaired_pdb": str(repaired_pdb) if repaired_pdb.exists() else "",
                "residue_count": row.get("residue_count", ""),
                "repair_exit_code": repair_code,
                "stability_exit_code": stability_code,
                "foldx_wt_total_energy": foldx_total_energy if foldx_total_energy is not None else "",
                "foldx_energy_per_residue": foldx_energy_per_residue,
                "parsed_foldx_file": parsed_file,
                "status": "ready" if foldx_total_energy is not None else "failed",
            }
        )

    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)

    print(f"FoldX WT rows: {len(result)}")
    print(f"Successful WT rows: {(result['status'] == 'ready').sum()}")
    print(f"Saved: {output_path}")

    if not (result["status"] == "ready").any():
        raise SystemExit("No FoldX WT stability result was successfully parsed.")


if __name__ == "__main__":
    main()
