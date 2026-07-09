#!/usr/bin/env python3
# Purpose: Merge salt-bridge features.

from pathlib import Path
import pandas as pd
import math
from collections import defaultdict

BASE = Path.home() / "xylanase-thesis"

MASTER_IN = BASE / "data/curated/xylanase_thesis_master.csv"
MASTER_OUT = BASE / "data/curated/xylanase_thesis_master_final_v3.csv"

SALT_OUT = BASE / "results/structures/salt_bridge_counts_all_structures.csv"
SUMMARY_OUT = BASE / "results/reports/xylanase_thesis_master_salt_bridge_summary.txt"

DIST_CUTOFF = 4.0

POSITIVE_ATOMS = {
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    "HIS": {"ND1", "NE2"},  # included as potentially protonated histidine
}

NEGATIVE_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

def parse_pdb_atoms(pdb_path):
    atoms = []

    try:
        with open(pdb_path, "r", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip() or "_"
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue

                icode = line[26].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                atoms.append({
                    "atom_name": atom_name,
                    "resname": resname,
                    "chain": chain,
                    "resseq": resseq,
                    "icode": icode,
                    "x": x,
                    "y": y,
                    "z": z,
                    "residue_id": (chain, resseq, icode, resname),
                })
    except FileNotFoundError:
        return None

    return atoms

def dist(a, b):
    return math.sqrt(
        (a["x"] - b["x"]) ** 2 +
        (a["y"] - b["y"]) ** 2 +
        (a["z"] - b["z"]) ** 2
    )

def salt_bridge_count(pdb_path):
    atoms = parse_pdb_atoms(pdb_path)
    if atoms is None:
        return None, "missing_pdb"
    if not atoms:
        return None, "no_atoms"

    positive = []
    negative = []

    for a in atoms:
        res = a["resname"]
        atom = a["atom_name"]

        if res in POSITIVE_ATOMS and atom in POSITIVE_ATOMS[res]:
            positive.append(a)

        if res in NEGATIVE_ATOMS and atom in NEGATIVE_ATOMS[res]:
            negative.append(a)

    residue_pairs = set()

    for p in positive:
        for n in negative:
            if p["residue_id"] == n["residue_id"]:
                continue

            if dist(p, n) <= DIST_CUTOFF:
                # Store residue-level pair only once even if multiple atom pairs satisfy cutoff
                pair = tuple(sorted([p["residue_id"], n["residue_id"]]))
                residue_pairs.add(pair)

    return len(residue_pairs), "ok"

def main():
    df = pd.read_csv(MASTER_IN)

    required_cols = ["uniprot_accession", "foldx_input_pdb", "has_any_structure"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column in final master: {c}")

    rows = []

    for _, r in df.iterrows():
        acc = str(r.get("uniprot_accession", "")).strip()
        pdb_path_raw = str(r.get("foldx_input_pdb", "")).strip()
        source_layer = str(r.get("foldx_source_layer", "")).strip()
        structure_source = str(r.get("structure_sources", "")).strip()
        structure_id = str(r.get("foldx_structure_id", "")).strip()

        if not pdb_path_raw or pdb_path_raw.lower() == "nan":
            rows.append({
                "uniprot_accession": acc,
                "structure_sources": structure_source,
                "foldx_source_layer": source_layer,
                "foldx_structure_id": structure_id,
                "foldx_input_pdb": pdb_path_raw,
                "salt_bridge_count_recomputed": None,
                "salt_bridge_status": "no_foldx_input_pdb",
            })
            continue

        pdb_path = Path(pdb_path_raw)

        count, status = salt_bridge_count(pdb_path)

        rows.append({
            "uniprot_accession": acc,
            "structure_sources": structure_source,
            "foldx_source_layer": source_layer,
            "foldx_structure_id": structure_id,
            "foldx_input_pdb": pdb_path_raw,
            "salt_bridge_count_recomputed": count,
            "salt_bridge_status": status,
        })

    salt = pd.DataFrame(rows)
    salt.to_csv(SALT_OUT, index=False)

    # Merge recomputed values back into final master
    df = df.merge(
        salt[["uniprot_accession", "salt_bridge_count_recomputed", "salt_bridge_status"]],
        on="uniprot_accession",
        how="left"
    )

    # Preserve old column for summary
    if "salt_bridge_count" in df.columns:
        df["salt_bridge_count_previous"] = df["salt_bridge_count"]
    else:
        df["salt_bridge_count_previous"] = ""

    df["salt_bridge_count"] = pd.to_numeric(df["salt_bridge_count_recomputed"], errors="coerce")

    # Per-residue normalized salt bridge count
    df["chain_length"] = pd.to_numeric(df["chain_length"], errors="coerce")
    df["salt_bridge_per_res"] = df["salt_bridge_count"] / df["chain_length"]

    df.to_csv(MASTER_OUT, index=False)

    summary_lines = []
    summary_lines.append("Xylanase Thesis Master Final v3 Salt-Bridge Summary")
    summary_lines.append("=" * 70)
    summary_lines.append(f"Input master: {MASTER_IN}")
    summary_lines.append(f"Salt-bridge output: {SALT_OUT}")
    summary_lines.append(f"Output master: {MASTER_OUT}")
    summary_lines.append("")
    summary_lines.append(f"Rows in master: {len(df)}")
    summary_lines.append(f"Unique UniProt accessions: {df['uniprot_accession'].nunique()}")
    summary_lines.append("")
    summary_lines.append("Salt-bridge recomputation status:")
    summary_lines.append(str(salt["salt_bridge_status"].value_counts(dropna=False)))
    summary_lines.append("")
    summary_lines.append("Salt-bridge non-null coverage after merge:")
    summary_lines.append(str(df["salt_bridge_count"].notna().value_counts(dropna=False)))
    summary_lines.append("")
    summary_lines.append("Salt-bridge coverage by structure source:")
    summary_lines.append(str(pd.crosstab(df["structure_sources"].fillna("none"), df["salt_bridge_count"].notna())))
    summary_lines.append("")
    summary_lines.append("Salt-bridge summary:")
    summary_lines.append(str(df["salt_bridge_count"].describe()))
    summary_lines.append("")
    summary_lines.append("Salt-bridge per residue summary:")
    summary_lines.append(str(df["salt_bridge_per_res"].describe()))
    summary_lines.append("")
    summary_lines.append("Top 20 salt-bridge counts:")
    top_cols = [
        "uniprot_accession", "organism", "organism_type", "gh_family",
        "structure_sources", "foldx_source_layer", "foldx_structure_id",
        "chain_length", "salt_bridge_count", "salt_bridge_per_res",
        "foldx_energy_per_residue"
    ]
    available_top_cols = [c for c in top_cols if c in df.columns]
    summary_lines.append(
        df.sort_values("salt_bridge_count", ascending=False)[available_top_cols]
        .head(20)
        .to_string(index=False)
    )

    SUMMARY_OUT.write_text("\n".join(summary_lines))

    print(f"Saved: {SALT_OUT}")
    print(f"Saved: {MASTER_OUT}")
    print(f"Saved: {SUMMARY_OUT}")
    print("")
    print("Salt-bridge status:")
    print(salt["salt_bridge_status"].value_counts(dropna=False))
    print("")
    print("Salt-bridge non-null coverage after merge:")
    print(df["salt_bridge_count"].notna().value_counts(dropna=False))
    print("")
    print("Coverage by structure source:")
    print(pd.crosstab(df["structure_sources"].fillna("none"), df["salt_bridge_count"].notna()))
    print("")
    print("Top 20 salt-bridge counts:")
    print(
        df.sort_values("salt_bridge_count", ascending=False)[available_top_cols]
        .head(20)
        .to_string(index=False)
    )

if __name__ == "__main__":
    main()
