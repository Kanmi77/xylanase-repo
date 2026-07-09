#!/usr/bin/env python3
"""
Prepare FoldX mutation screening for clean unique PDB xylanases.

This script:
1. Uses original PDB structured subset only.
2. Keeps bacterial/fungal GH10/GH11 only.
3. Selects one best PDB structure per UniProt accession.
4. Extracts PDB chain sequences from FoldX-ready PDB files.
5. Runs MAFFT per group: bacterial GH10, bacterial GH11, fungal GH10, fungal GH11.
6. Uses group consensus to propose candidate mutations.
7. Creates per-structure FoldX individual_list.txt files.
8. Creates a shell script to run FoldX BuildModel.

Important:
- This uses PDB chain residue numbering, not UniProt numbering.
- Mutations involving E/C residues are excluded by default to reduce risk of catalytic/disulfide disruption.
- Frozen/Snakemake files are ignored.

Run:
    cd ~/xylanase-thesis
    python scripts/10_optionC_original_only/05_prepare_pdb60_consensus_foldx_mutations.py

Then run:
    bash results/optionC_original_only/pdb60_foldx_mutations/01_run_pdb60_foldx_buildmodel.sh
"""

from __future__ import annotations

import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_csv_first(paths):
    for f in paths:
        p = Path(f)
        if p.exists() and "frozen" not in str(p).lower():
            return p, pd.read_csv(p, low_memory=False)
    raise FileNotFoundError("None of the expected PDB structured subset files was found.")


def safe_id(x):
    return (
        str(x)
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
        .replace(";", "")
        .replace(":", "_")
    )


def infer_structure_id(row):
    if "structure_id" in row.index and pd.notna(row["structure_id"]):
        return str(row["structure_id"])

    if "pdb_tag" in row.index and pd.notna(row["pdb_tag"]):
        return str(row["pdb_tag"])

    pdb_id = str(row.get("pdb_id", "")).strip()
    chain = str(row.get("chosen_chain", "")).strip()

    if pdb_id and chain:
        return f"{pdb_id}_{chain}_Repair"

    return ""


def resolve_foldx_pdb_path(root: Path, row: pd.Series) -> str:
    for c in ["foldx_pdb", "foldx_input_pdb", "selected_structure_path", "structure_path", "file_path"]:
        if c in row.index and pd.notna(row[c]):
            p = Path(str(row[c]))
            if p.exists() and p.suffix.lower() == ".pdb":
                return str(p)

    pdb_id = str(row.get("pdb_id", "")).strip()
    chain = str(row.get("chosen_chain", row.get("target_chain", ""))).strip()

    if pdb_id and chain:
        candidates = [
            root / "foldx" / "pdb_chains" / f"{pdb_id}_{chain}.pdb",
            root / "foldx" / "pdb_chains" / f"{pdb_id}_{chain}_Repair.pdb",
            root / "data" / "processed" / "pdb_chains" / f"{pdb_id}_{chain}.pdb",
        ]

        for p in candidates:
            if p.exists():
                return str(p)

    return ""


def parse_pdb_chain(pdb_path: str, chain_filter: str):
    p = Path(str(pdb_path))

    if not p.exists():
        return "", pd.DataFrame()

    residues = []
    seen = set()

    with p.open("r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            chain = line[21].strip() or "_"

            if chain_filter and chain != chain_filter:
                continue

            resname = line[17:20].strip().upper()
            resnum = line[22:26].strip()
            icode = line[26].strip()

            key = (chain, resnum, icode)

            if key in seen:
                continue

            seen.add(key)

            aa = AA3_TO_1.get(resname, "X")

            residues.append({
                "chain_id": chain,
                "pdb_residue_number": resnum,
                "insertion_code": icode,
                "residue_name_3letter": resname,
                "wt_residue": aa,
            })

    df = pd.DataFrame(residues)

    if df.empty:
        return "", df

    df["chain_sequence_index"] = np.arange(1, len(df) + 1)
    seq = "".join(df["wt_residue"].tolist())

    return seq, df


def write_fasta(records, out_fasta):
    with open(out_fasta, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


def read_fasta(path):
    records = {}
    name = None
    seqs = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if name:
                    records[name] = "".join(seqs)
                name = line[1:].strip().split()[0]
                seqs = []
            else:
                seqs.append(line.strip())

        if name:
            records[name] = "".join(seqs)

    return records


def run_mafft(in_fasta, out_fasta):
    mafft = shutil.which("mafft")

    if not mafft:
        raise RuntimeError("MAFFT not found. Please activate your environment or install mafft.")

    cmd = [mafft, "--auto", str(in_fasta)]

    with open(out_fasta, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)


def build_alignment_to_residue_map(aligned_seq, residue_df):
    """
    Maps alignment column index to PDB residue row.
    Returns dict: aln_index -> residue_df row dict
    """

    mapping = {}
    residue_idx = 0
    residues = residue_df.to_dict("records")

    for aln_idx, aa in enumerate(aligned_seq):
        if aa == "-":
            continue

        if residue_idx >= len(residues):
            break

        mapping[aln_idx] = residues[residue_idx]
        residue_idx += 1

    return mapping


def is_safe_mutation(wt, mut):
    if wt not in VALID_AA:
        return False

    if mut not in VALID_AA:
        return False

    if wt == mut:
        return False

    # Avoid cysteine to reduce disulfide disruption / accidental disulfide creation.
    if wt == "C" or mut == "C":
        return False

    # Avoid glutamate because GH10/GH11 catalytic residues are usually glutamates.
    # This is conservative; we can relax later after active-site annotation.
    if wt == "E" or mut == "E":
        return False

    return True


def main():
    root = Path(".").resolve()

    outdir = root / "results" / "optionC_original_only" / "pdb60_foldx_mutations"
    fasta_dir = outdir / "group_fastas"
    aln_dir = outdir / "group_alignments"
    input_dir = outdir / "foldx_inputs"
    work_dir = outdir / "foldx_buildmodel_runs"

    for d in [outdir, fasta_dir, aln_dir, input_dir, work_dir]:
        d.mkdir(parents=True, exist_ok=True)

    pdb_file, pdb = read_csv_first([
        "data/curated/xylanase_structured_subset_with_foldx_norm.csv",
        "data/curated/xylanase_structured_subset_with_foldx.csv",
        "data/curated/xylanase_structured_subset.csv",
    ])

    print(f"[INFO] PDB file used: {pdb_file}")

    clean = pdb[
        pdb["organism_type"].isin(["bacterial", "fungal"]) &
        pdb["gh_family"].isin(["GH10", "GH11"])
    ].copy()

    if "foldx_energy_per_residue" not in clean.columns:
        if "foldx_wt_total_energy" in clean.columns and "chain_length" in clean.columns:
            clean["foldx_energy_per_residue"] = (
                pd.to_numeric(clean["foldx_wt_total_energy"], errors="coerce") /
                pd.to_numeric(clean["chain_length"], errors="coerce")
            )
        else:
            clean["foldx_energy_per_residue"] = np.nan

    clean["structure_id"] = clean.apply(infer_structure_id, axis=1)
    clean["foldx_ready_pdb_path"] = clean.apply(lambda r: resolve_foldx_pdb_path(root, r), axis=1)
    clean["foldx_ready_pdb_exists"] = clean["foldx_ready_pdb_path"].map(lambda x: bool(x) and Path(str(x)).exists())

    clean = clean[clean["foldx_ready_pdb_exists"]].copy()

    # One best PDB structure per UniProt accession.
    # If multiple PDB structures represent the same accession, select the one with lowest FoldX energy/residue.
    reps = (
        clean.sort_values(["uniprot_accession", "foldx_energy_per_residue"], ascending=[True, True])
        .drop_duplicates("uniprot_accession", keep="first")
        .copy()
    )

    reps = reps.sort_values(["organism_type", "gh_family", "foldx_energy_per_residue"], ascending=[True, True, True])
    reps.to_csv(outdir / "clean_unique_pdb60_representatives.csv", index=False)

    print(f"[INFO] Clean PDB rows: {len(clean)}")
    print(f"[INFO] Clean unique PDB representatives: {len(reps)}")

    all_residue_maps = []
    sequence_records = []

    for _, row in reps.iterrows():
        sid = safe_id(row["structure_id"])
        chain = str(row.get("chosen_chain", row.get("target_chain", ""))).strip()
        pdb_path = row["foldx_ready_pdb_path"]

        seq, residues = parse_pdb_chain(pdb_path, chain)

        if not seq or residues.empty:
            continue

        residues["structure_id"] = sid
        residues["pdb_id"] = row.get("pdb_id", "")
        residues["uniprot_accession"] = row.get("uniprot_accession", "")
        residues["organism_type"] = row.get("organism_type", "")
        residues["gh_family"] = row.get("gh_family", "")
        residues["foldx_ready_pdb_path"] = pdb_path

        all_residue_maps.append(residues)

        sequence_records.append({
            "structure_id": sid,
            "sequence": seq,
            "organism_type": row.get("organism_type", ""),
            "gh_family": row.get("gh_family", ""),
            "uniprot_accession": row.get("uniprot_accession", ""),
            "pdb_id": row.get("pdb_id", ""),
            "chosen_chain": chain,
            "foldx_ready_pdb_path": pdb_path,
        })

    residue_map = pd.concat(all_residue_maps, ignore_index=True)
    seq_df = pd.DataFrame(sequence_records)

    residue_map.to_csv(outdir / "pdb60_chain_residue_map.csv", index=False)
    seq_df.to_csv(outdir / "pdb60_chain_sequences.csv", index=False)

    mutation_rows = []

    for (organism_type, gh_family), group in seq_df.groupby(["organism_type", "gh_family"]):
        group_name = f"{organism_type}_{gh_family}"
        fasta = fasta_dir / f"{group_name}.fasta"
        aligned = aln_dir / f"{group_name}.aligned.fasta"

        records = [(r["structure_id"], r["sequence"]) for _, r in group.iterrows()]
        write_fasta(records, fasta)

        if len(records) < 2:
            print(f"[WARN] Skipping {group_name}: fewer than 2 sequences.")
            continue

        print(f"[INFO] Running MAFFT for {group_name}: {len(records)} sequences")
        run_mafft(fasta, aligned)

        aln_records = read_fasta(aligned)

        aln_len = max(len(s) for s in aln_records.values())

        # Consensus per alignment column.
        consensus = {}

        for i in range(aln_len):
            chars = []

            for seq in aln_records.values():
                if i < len(seq):
                    aa = seq[i]
                    if aa != "-" and aa in VALID_AA:
                        chars.append(aa)

            if not chars:
                continue

            counts = Counter(chars)
            aa, n = counts.most_common(1)[0]
            freq = n / len(chars)

            consensus[i] = {
                "consensus_residue": aa,
                "consensus_count": n,
                "non_gap_count": len(chars),
                "consensus_frequency": freq,
            }

        # Build per-sequence alignment-to-PDB mapping.
        for _, seqrow in group.iterrows():
            sid = seqrow["structure_id"]
            aligned_seq = aln_records.get(sid)

            if not aligned_seq:
                continue

            resdf = residue_map[residue_map["structure_id"] == sid].copy()
            aln_to_res = build_alignment_to_residue_map(aligned_seq, resdf)

            candidates = []

            for aln_idx, res in aln_to_res.items():
                if aln_idx not in consensus:
                    continue

                con = consensus[aln_idx]
                wt = res["wt_residue"]
                mut = con["consensus_residue"]

                if con["consensus_frequency"] < 0.60:
                    continue

                if not is_safe_mutation(wt, mut):
                    continue

                # Avoid extreme termini.
                seq_idx = int(res["chain_sequence_index"])
                seq_len = len(resdf)

                if seq_idx <= 5 or seq_idx >= seq_len - 5:
                    continue

                pdb_num = str(res["pdb_residue_number"])
                chain = str(res["chain_id"])

                if not pdb_num.replace("-", "").isdigit():
                    continue

                foldx_code = f"{wt}{chain}{pdb_num}{mut};"

                candidates.append({
                    "structure_id": sid,
                    "uniprot_accession": seqrow["uniprot_accession"],
                    "pdb_id": seqrow["pdb_id"],
                    "chain_id": chain,
                    "pdb_residue_number": pdb_num,
                    "chain_sequence_index": seq_idx,
                    "organism_type": organism_type,
                    "gh_family": gh_family,
                    "wt_residue": wt,
                    "mutant_residue": mut,
                    "foldx_mutation_code": foldx_code,
                    "alignment_column": aln_idx + 1,
                    "consensus_frequency": con["consensus_frequency"],
                    "consensus_count": con["consensus_count"],
                    "non_gap_count": con["non_gap_count"],
                    "foldx_ready_pdb_path": seqrow["foldx_ready_pdb_path"],
                    "mutation_reason": "PDB-group consensus-supported substitution",
                })

            # Keep top 15 consensus-supported mutations per structure.
            candidates = sorted(
                candidates,
                key=lambda x: (
                    x["consensus_frequency"],
                    x["consensus_count"],
                    -abs(x["chain_sequence_index"] - (len(resdf) / 2)),
                ),
                reverse=True,
            )[:15]

            mutation_rows.extend(candidates)

    mutations = pd.DataFrame(mutation_rows)

    if mutations.empty:
        raise RuntimeError("No mutation candidates were generated. Check MAFFT/group sequences.")

    mutations.to_csv(outdir / "pdb60_consensus_mutation_candidates.csv", index=False)

    # Create per-structure FoldX input files and BuildModel run script.
    run_script = outdir / "01_run_pdb60_foldx_buildmodel.sh"

    with run_script.open("w") as sh:
        sh.write("#!/usr/bin/env bash\n")
        sh.write("set -euo pipefail\n\n")
        sh.write("FOLDX_BIN=${FOLDX_BIN:-foldx}\n\n")
        sh.write("echo \"Using FoldX binary: ${FOLDX_BIN}\"\n\n")

        for sid, sub in mutations.groupby("structure_id"):
            sid_safe = safe_id(sid)

            pdb_path = sub["foldx_ready_pdb_path"].iloc[0]
            pdb_name = Path(pdb_path).name

            target_dir = work_dir / sid_safe
            target_dir.mkdir(parents=True, exist_ok=True)

            individual_list = target_dir / "individual_list.txt"

            with individual_list.open("w") as fh:
                for code in sub["foldx_mutation_code"].drop_duplicates():
                    fh.write(code + "\n")

            sh.write(f"echo 'Running FoldX BuildModel for {sid_safe}'\n")
            sh.write(f"mkdir -p '{target_dir}'\n")
            sh.write(f"cp '{pdb_path}' '{target_dir}/{pdb_name}'\n")
            sh.write(f"cd '{target_dir}'\n")
            sh.write(
                f"${{FOLDX_BIN}} "
                f"--command=BuildModel "
                f"--pdb='{pdb_name}' "
                f"--mutant-file='individual_list.txt' "
                f"--numberOfRuns=1 "
                f"--out-pdb=true\n"
            )
            sh.write(f"cd '{root}'\n\n")

    run_script.chmod(0o755)

    # Summary by group.
    summary = (
        mutations
        .groupby(["organism_type", "gh_family"])
        .agg(
            structures=("structure_id", "nunique"),
            mutation_candidates=("foldx_mutation_code", "count"),
            median_consensus_frequency=("consensus_frequency", "median"),
        )
        .reset_index()
    )

    summary.to_csv(outdir / "pdb60_mutation_candidate_summary_by_group.csv", index=False)

    report = outdir / "PDB60_FOLDX_MUTATION_PREP_REPORT.md"

    with report.open("w") as fh:
        fh.write("# PDB60 FoldX mutation preparation report\n\n")
        fh.write(f"- Original PDB file used: `{pdb_file}`\n")
        fh.write(f"- Clean PDB rows: {len(clean)}\n")
        fh.write(f"- Clean unique PDB representatives: {len(reps)}\n")
        fh.write(f"- PDB chain sequences extracted: {len(seq_df)}\n")
        fh.write(f"- Residue-map rows: {len(residue_map)}\n")
        fh.write(f"- Mutation candidates generated: {len(mutations)}\n")
        fh.write(f"- Structures with mutation candidates: {mutations['structure_id'].nunique()}\n")
        fh.write(f"- FoldX run script: `{run_script}`\n\n")
        fh.write("## Mutation design rule\n\n")
        fh.write(
            "Mutations were generated from group-specific PDB-chain consensus alignments. "
            "Only bacterial/fungal GH10/GH11 PDB structures were used. One representative "
            "PDB structure was selected per UniProt accession. Mutations involving cysteine "
            "or glutamate were excluded conservatively to reduce risk of disulfide or catalytic "
            "residue disruption.\n"
        )

    print("\n[DONE] PDB60 FoldX mutation preparation complete.")
    print(f"Report: {report}")
    print(f"Mutation candidates: {outdir / 'pdb60_consensus_mutation_candidates.csv'}")
    print(f"FoldX run script: {run_script}")


if __name__ == "__main__":
    main()
