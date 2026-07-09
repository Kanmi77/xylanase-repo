#!/usr/bin/env python3
"""
Reconcile PDB structures for mutation screening and docking.

This script ignores all files with 'frozen' in the name because those belong
to the Snakemake pipeline workstream, not the original thesis dataset.

It will:
1. Load PDB candidates from the original Option C/PDB files.
2. Resolve FoldX-ready PDB chain paths.
3. Reconcile existing WT docking results.
4. Identify missing WT docking jobs.
5. Extract PDB-chain residue maps for FoldX mutation numbering.
6. Create a manual mutation design template.
7. Optionally validate filled mutations and convert them to FoldX syntax.

Run:
    cd ~/xylanase-thesis
    python scripts/10_optionC_original_only/03_reconcile_pdb_mutation_docking.py

After filling the mutation template:
    python scripts/10_optionC_original_only/03_reconcile_pdb_mutation_docking.py \
      --manual-mutations results/optionC_original_only/pdb_reconciliation/pdb_manual_mutation_design_filled.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


def is_frozen(path: Path) -> bool:
    return "frozen" in str(path).lower()


def read_csv(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists() or is_frozen(path):
        return None
    try:
        return pd.read_csv(path, low_memory=False)
    except Exception as exc:
        print(f"[WARN] Could not read {path}: {exc}")
        return None


def first_existing(root: Path, candidates: list[str]) -> Optional[Path]:
    for rel in candidates:
        p = root / rel
        if p.exists() and not is_frozen(p):
            return p
    return None


def infer_structure_id(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    if "structure_id" not in df.columns:
        if "pdb_tag" in df.columns:
            df["structure_id"] = df["pdb_tag"]
        elif "pdb_id" in df.columns and "chosen_chain" in df.columns:
            df["structure_id"] = (
                df["pdb_id"].astype(str)
                + "_"
                + df["chosen_chain"].astype(str)
                + "_Repair"
            )

    if "pdb_tag" not in df.columns and "structure_id" in df.columns:
        df["pdb_tag"] = df["structure_id"]

    return df


def load_pdb_candidates(root: Path, max_pdb: int) -> pd.DataFrame:
    candidate_files = [
        "results/optionC_original_only/optionC_original_top50_pdb_reference_structures.csv",
        "results/optionC_original_only/optionC_original_ranked_structures_source_aware.csv",
        "results/optionC_original_only/optionC_original_mutation_and_docking_targets.csv",

        "results/foldx/top15_candidates_with_docking.csv",
        "results/foldx/top15_candidates.csv",
        "top15_candidates_with_docking.csv",
        "top15_candidates.csv",

        "data/curated/xylanase_structured_subset_with_foldx_norm.csv",
        "data/curated/xylanase_structured_subset_with_foldx.csv",
        "data/curated/xylanase_structured_subset.csv",
    ]

    path = first_existing(root, candidate_files)

    if path is None:
        raise FileNotFoundError(
            "No original PDB candidate table found. "
            "Run Option C original-only stage 00 first or check original PDB files."
        )

    df = read_csv(path)

    if df is None or df.empty:
        raise ValueError(f"Could not read candidate table: {path}")

    df = infer_structure_id(df)

    # If source column exists, keep only PDB rows.
    # If no source column exists, we assume this is an original PDB-specific file.
    if "structure_source_norm" in df.columns:
        df = df[df["structure_source_norm"].astype(str).str.lower().eq("pdb")]
    elif "structure_source" in df.columns:
        df = df[df["structure_source"].astype(str).str.lower().eq("pdb")]

    df = df.copy()
    df["pdb_candidate_source_file"] = str(path.relative_to(root))

    sort_cols = [
        c for c in [
            "optionC_structural_score_no_docking",
            "final_score",
            "foldx_energy_per_residue",
        ]
        if c in df.columns
    ]

    if sort_cols:
        ascending = []
        for c in sort_cols:
            if c == "foldx_energy_per_residue":
                ascending.append(True)
            else:
                ascending.append(False)
        df = df.sort_values(sort_cols, ascending=ascending)

    return df.head(max_pdb).reset_index(drop=True)


def resolve_foldx_pdb_path(root: Path, row: pd.Series) -> str:
    # Prefer existing FoldX-ready PDB paths if present.
    for c in [
        "foldx_pdb",
        "foldx_input_pdb",
        "selected_structure_path",
        "structure_path",
        "file_path",
    ]:
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


def parse_pdb_residues(pdb_path: str, chain_filter: Optional[str] = None) -> pd.DataFrame:
    p = Path(str(pdb_path))

    if not p.exists() or p.suffix.lower() != ".pdb":
        return pd.DataFrame()

    residues = []
    seen = set()

    with p.open("r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            chain = line[21].strip() or "_"

            if chain_filter and chain != str(chain_filter):
                continue

            resname = line[17:20].strip().upper()
            resseq = line[22:26].strip()
            icode = line[26].strip()

            key = (chain, resseq, icode)

            if key in seen:
                continue

            seen.add(key)

            residues.append({
                "chain_id": chain,
                "pdb_residue_number": resseq,
                "insertion_code": icode,
                "residue_name_3letter": resname,
                "wt_residue": AA3_TO_1.get(resname, "X"),
            })

    out = pd.DataFrame(residues)

    if not out.empty:
        out["chain_sequence_index"] = out.groupby("chain_id").cumcount() + 1

    return out


def find_docking_tables(root: Path) -> list[Path]:
    paths = []

    for base in [root / "results", root]:
        if not base.exists():
            continue

        for p in base.rglob("*.csv"):
            if is_frozen(p):
                continue

            name = str(p).lower()

            if "dock" not in name and "vina" not in name and "top15" not in name:
                continue

            try:
                cols = pd.read_csv(p, nrows=0).columns.tolist()
            except Exception:
                continue

            vina_cols = [
                "vina_best_xylobiose",
                "vina_best_xylotriose",
                "vina_best_min",
                "vina_best_mean",
            ]

            if any(c in cols for c in vina_cols):
                paths.append(p)

    return sorted(set(paths))


def load_existing_docking(root: Path) -> pd.DataFrame:
    parts = []

    for p in find_docking_tables(root):
        df = read_csv(p)

        if df is None or df.empty:
            continue

        df = infer_structure_id(df)
        df["docking_source_file"] = str(p.relative_to(root))
        parts.append(df)

    if not parts:
        return pd.DataFrame()

    out = pd.concat(parts, ignore_index=True, sort=False)

    keys = [
        c for c in [
            "structure_id",
            "pdb_tag",
            "pdb_id",
            "chosen_chain",
            "uniprot_accession",
        ]
        if c in out.columns
    ]

    if keys:
        out = out.drop_duplicates(keys, keep="first")

    return out


def reconcile_docking(pdb: pd.DataFrame, dock: pd.DataFrame) -> pd.DataFrame:
    pdb = infer_structure_id(pdb)

    if dock.empty:
        out = pdb.copy()
        out["has_existing_wt_docking"] = False
        return out

    dock = infer_structure_id(dock)

    docking_cols = [
        c for c in [
            "structure_id",
            "pdb_tag",
            "pdb_id",
            "chosen_chain",
            "uniprot_accession",
            "vina_best_xylobiose",
            "vina_best_xylotriose",
            "vina_best_min",
            "vina_best_mean",
            "docking_source_file",
        ]
        if c in dock.columns
    ]

    if "structure_id" in pdb.columns and "structure_id" in dock.columns:
        out = pdb.merge(
            dock[docking_cols].drop_duplicates("structure_id"),
            on="structure_id",
            how="left",
            suffixes=("", "_dock"),
        )
    elif "pdb_tag" in pdb.columns and "pdb_tag" in dock.columns:
        out = pdb.merge(
            dock[docking_cols].drop_duplicates("pdb_tag"),
            on="pdb_tag",
            how="left",
            suffixes=("", "_dock"),
        )
    else:
        out = pdb.copy()

    vina_cols = [
        c for c in [
            "vina_best_xylobiose",
            "vina_best_xylotriose",
            "vina_best_min",
            "vina_best_mean",
        ]
        if c in out.columns
    ]

    if vina_cols:
        out["has_existing_wt_docking"] = out[vina_cols].notna().any(axis=1)
    else:
        out["has_existing_wt_docking"] = False

    return out


def validate_manual_mutations(
    outdir: Path,
    manual_path: Path,
    residue_map: pd.DataFrame,
) -> pd.DataFrame:
    manual = read_csv(manual_path)

    if manual is None or manual.empty:
        raise ValueError(f"Manual mutation file could not be read: {manual_path}")

    required = [
        "structure_id",
        "chain_id",
        "pdb_residue_number",
        "wt_residue",
        "mutant_residue",
    ]

    missing = [c for c in required if c not in manual.columns]

    if missing:
        raise KeyError(f"Manual mutation file missing columns: {missing}")

    if residue_map.empty:
        raise ValueError(
            "Residue map is empty, so mutations cannot be validated. "
            "Check whether FoldX-ready PDB chain files exist."
        )

    m = manual.copy()
    m["pdb_residue_number"] = m["pdb_residue_number"].astype(str)
    residue_map["pdb_residue_number"] = residue_map["pdb_residue_number"].astype(str)

    check = m.merge(
        residue_map[
            [
                "structure_id",
                "chain_id",
                "pdb_residue_number",
                "wt_residue",
            ]
        ].rename(columns={"wt_residue": "wt_residue_from_pdb"}),
        on=[
            "structure_id",
            "chain_id",
            "pdb_residue_number",
        ],
        how="left",
    )

    check["wt_residue"] = check["wt_residue"].astype(str).str.upper().str.strip()
    check["mutant_residue"] = check["mutant_residue"].astype(str).str.upper().str.strip()
    check["wt_residue_from_pdb"] = check["wt_residue_from_pdb"].astype(str).str.upper().str.strip()

    check["wt_matches_pdb"] = check["wt_residue"] == check["wt_residue_from_pdb"]

    def make_foldx_code(row):
        if not row["wt_matches_pdb"]:
            return "CHECK_WT_RESIDUE"
        return (
            f"{row['wt_residue']}"
            f"{row['chain_id']}"
            f"{row['pdb_residue_number']}"
            f"{row['mutant_residue']};"
        )

    check["foldx_mutation_code"] = check.apply(make_foldx_code, axis=1)

    check.to_csv(
        outdir / "pdb_manual_mutations_validated_with_foldx_codes.csv",
        index=False,
    )

    with (outdir / "pdb_foldx_individual_list_template.txt").open("w") as fh:
        for sid, sub in check.groupby("structure_id"):
            valid_codes = (
                sub.loc[sub["wt_matches_pdb"] == True, "foldx_mutation_code"]
                .dropna()
                .tolist()
            )

            if valid_codes:
                fh.write(f"# {sid}\n")
                for code in valid_codes:
                    fh.write(code + "\n")

    return check


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root. Default: current directory.",
    )

    parser.add_argument(
        "--max-pdb",
        type=int,
        default=40,
        help="Maximum number of PDB targets to reconcile.",
    )

    parser.add_argument(
        "--manual-mutations",
        default=None,
        help="Optional filled manual mutation CSV to validate and convert to FoldX syntax.",
    )

    args = parser.parse_args()

    root = Path(args.project_root).resolve()

    outdir = root / "results" / "optionC_original_only" / "pdb_reconciliation"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Project root: {root}")
    print("[INFO] Frozen/Snakemake files will be ignored.")
    print(f"[INFO] Output folder: {outdir}")

    pdb = load_pdb_candidates(root, args.max_pdb)

    pdb["foldx_ready_pdb_path"] = pdb.apply(
        lambda row: resolve_foldx_pdb_path(root, row),
        axis=1,
    )

    pdb["foldx_ready_pdb_exists"] = pdb["foldx_ready_pdb_path"].map(
        lambda x: bool(x) and Path(str(x)).exists()
    )

    pdb.to_csv(
        outdir / "pdb_targets_resolved_for_mutation_docking.csv",
        index=False,
    )

    docking = load_existing_docking(root)

    if not docking.empty:
        docking.to_csv(
            outdir / "existing_pdb_docking_tables_combined.csv",
            index=False,
        )

    reconciled = reconcile_docking(pdb, docking)

    reconciled.to_csv(
        outdir / "pdb_wt_docking_reconciled.csv",
        index=False,
    )

    missing_docking = reconciled[
        ~reconciled["has_existing_wt_docking"].astype(bool)
    ].copy()

    missing_docking.to_csv(
        outdir / "pdb_missing_wt_docking_jobs.csv",
        index=False,
    )

    residue_maps = []

    for _, row in pdb.iterrows():
        pdb_path = row.get("foldx_ready_pdb_path", "")
        chain = row.get("chosen_chain", row.get("target_chain", None))

        if pd.notna(chain):
            chain_filter = str(chain).strip()
        else:
            chain_filter = None

        residues = parse_pdb_residues(
            pdb_path,
            chain_filter=chain_filter,
        )

        if residues.empty:
            continue

        residues["structure_id"] = row.get("structure_id", row.get("pdb_tag", ""))
        residues["pdb_id"] = row.get("pdb_id", "")
        residues["uniprot_accession"] = row.get("uniprot_accession", "")

        residue_maps.append(residues)

    if residue_maps:
        residue_map = pd.concat(residue_maps, ignore_index=True)
    else:
        residue_map = pd.DataFrame()

    residue_map.to_csv(
        outdir / "pdb_chain_residue_map_for_foldx_numbering.csv",
        index=False,
    )

    template_cols = [
        "structure_id",
        "pdb_id",
        "uniprot_accession",
        "chain_id",
        "pdb_residue_number",
        "wt_residue",
        "mutant_residue",
        "mutation_reason",
        "avoid_if_catalytic_or_binding_site",
    ]

    if not pdb.empty:
        if not residue_map.empty:
            template = residue_map[
                [
                    "structure_id",
                    "pdb_id",
                    "uniprot_accession",
                    "chain_id",
                    "pdb_residue_number",
                    "wt_residue",
                ]
            ].copy()

            template["mutant_residue"] = ""
            template["mutation_reason"] = ""
            template["avoid_if_catalytic_or_binding_site"] = ""

            template = template[template_cols]
        else:
            template = pdb[
                [
                    c for c in [
                        "structure_id",
                        "pdb_id",
                        "uniprot_accession",
                        "chosen_chain",
                    ]
                    if c in pdb.columns
                ]
            ].copy()

            if "chosen_chain" in template.columns:
                template = template.rename(columns={"chosen_chain": "chain_id"})

            for c in template_cols:
                if c not in template.columns:
                    template[c] = ""

            template = template[template_cols]

        template.to_csv(
            outdir / "pdb_manual_mutation_design_template.csv",
            index=False,
        )

    if args.manual_mutations:
        validated = validate_manual_mutations(
            outdir=outdir,
            manual_path=Path(args.manual_mutations),
            residue_map=residue_map,
        )

        print(f"[INFO] Validated manual mutations: {len(validated)}")

    report_lines = []

    report_lines.append("# PDB mutation/docking reconciliation report\n\n")
    report_lines.append("Frozen/Snakemake files were ignored.\n\n")
    report_lines.append(f"- PDB targets reconciled: {len(pdb)}\n")
    report_lines.append(
        f"- FoldX-ready PDB paths found: "
        f"{int(pdb['foldx_ready_pdb_exists'].sum()) if 'foldx_ready_pdb_exists' in pdb else 0}\n"
    )
    report_lines.append(
        f"- Existing WT docking hits found: "
        f"{int(reconciled['has_existing_wt_docking'].sum()) if 'has_existing_wt_docking' in reconciled else 0}\n"
    )
    report_lines.append(f"- Missing WT docking jobs: {len(missing_docking)}\n")
    report_lines.append(f"- PDB residue-map rows: {len(residue_map)}\n\n")

    report_lines.append("## Interpretation\n\n")
    report_lines.append(
        "PDB mutation screening must use PDB chain residue numbering, not UniProt numbering, "
        "unless an explicit UniProt-to-PDB mapping is created. The generated residue map is the "
        "safe basis for FoldX mutation codes.\n"
    )

    (outdir / "PDB_MUTATION_DOCKING_RECONCILIATION_REPORT.md").write_text(
        "".join(report_lines),
        encoding="utf-8",
    )

    print("\n[DONE] PDB mutation/docking reconciliation complete.")
    print(f"Report: {outdir / 'PDB_MUTATION_DOCKING_RECONCILIATION_REPORT.md'}")


if __name__ == "__main__":
    main()
