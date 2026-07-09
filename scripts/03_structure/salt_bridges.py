#!/usr/bin/env python3
# Purpose: Calculate salt-bridge features.

from pathlib import Path
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

INPUT = Path("results/foldx/foldx_normalized.csv")
OUTDIR = Path("results/structural_features_saltbridges")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUTDIR / "foldx_normalized_with_salt_bridges.csv"
OUT_REPORT = OUTDIR / "SALT_BRIDGE_SUMMARY_REPORT.md"

ACIDIC = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

BASIC = {
    "LYS": ["NZ"],
    "ARG": ["NH1", "NH2", "NE"],
    "HIS": ["ND1", "NE2"],
    "HSD": ["ND1", "NE2"],
    "HSE": ["ND1", "NE2"],
    "HSP": ["ND1", "NE2"],
}

DIST_CUTOFF = 4.0

parser = PDBParser(QUIET=True)

df = pd.read_csv(INPUT)

def count_salt_bridges(pdb_path):
    pdb_path = Path(str(pdb_path))
    if not pdb_path.exists():
        return np.nan

    try:
        structure = parser.get_structure("protein", pdb_path)
    except Exception:
        return np.nan

    acidic_atoms = []
    basic_atoms = []

    for atom in structure.get_atoms():
        residue = atom.get_parent()
        resname = residue.get_resname().strip().upper()
        atom_name = atom.get_name().strip().upper()

        if resname in ACIDIC and atom_name in ACIDIC[resname]:
            acidic_atoms.append(atom)

        if resname in BASIC and atom_name in BASIC[resname]:
            basic_atoms.append(atom)

    if not acidic_atoms or not basic_atoms:
        return 0

    ns = NeighborSearch(basic_atoms)
    pairs = set()

    for a in acidic_atoms:
        nearby = ns.search(a.coord, DIST_CUTOFF, level="A")
        for b in nearby:
            res_a = a.get_parent()
            res_b = b.get_parent()

            # avoid counting contacts within same residue, though acidic/basic same residue is unlikely
            if res_a == res_b:
                continue

            chain_a = res_a.get_parent().id
            chain_b = res_b.get_parent().id
            id_a = res_a.get_id()
            id_b = res_b.get_id()

            pair_id = (
                chain_a,
                id_a[1],
                res_a.get_resname(),
                chain_b,
                id_b[1],
                res_b.get_resname(),
            )
            pairs.add(pair_id)

    return len(pairs)

salt_counts = []

for i, row in df.iterrows():
    pdb_path = row.get("foldx_input_pdb", None)

    # Prefer repaired PDB if available in structure_id_y naming
    # Otherwise use foldx_input_pdb from the table.
    count = count_salt_bridges(pdb_path)
    salt_counts.append(count)

df["salt_bridge_count"] = salt_counts
df["salt_bridge_per_res"] = df["salt_bridge_count"] / df["chain_length"]

df.to_csv(OUT_TABLE, index=False)

summary = df["salt_bridge_count"].describe()
summary_per_res = df["salt_bridge_per_res"].describe()

group_summary = (
    df.groupby(["organism_type", "gh_family"])
    .agg(
        records=("uniprot_accession", "count"),
        mean_chain_length=("chain_length", "mean"),
        mean_salt_bridge_count=("salt_bridge_count", "mean"),
        median_salt_bridge_count=("salt_bridge_count", "median"),
        mean_salt_bridge_per_res=("salt_bridge_per_res", "mean"),
        median_salt_bridge_per_res=("salt_bridge_per_res", "median"),
    )
    .reset_index()
)

group_summary.to_csv(OUTDIR / "salt_bridge_group_summary.csv", index=False)

with open(OUT_REPORT, "w") as h:
    h.write("# Salt-bridge structural feature summary\n\n")
    h.write(f"- Input file: `{INPUT}`\n")
    h.write(f"- Output table: `{OUT_TABLE}`\n")
    h.write(f"- Distance cutoff: {DIST_CUTOFF} Å\n")
    h.write(f"- Records analysed: {len(df)}\n")
    h.write(f"- Records with salt-bridge values: {df['salt_bridge_count'].notna().sum()}\n\n")

    h.write("## Salt-bridge count summary\n\n")
    h.write(summary.to_frame("salt_bridge_count").to_markdown())
    h.write("\n\n")

    h.write("## Salt-bridge per-residue summary\n\n")
    h.write(summary_per_res.to_frame("salt_bridge_per_res").to_markdown())
    h.write("\n\n")

    h.write("## Salt-bridge group summary\n\n")
    h.write(group_summary.to_markdown(index=False))
    h.write("\n")

print("Done")
print("Output:", OUT_TABLE)
print("Report:", OUT_REPORT)
