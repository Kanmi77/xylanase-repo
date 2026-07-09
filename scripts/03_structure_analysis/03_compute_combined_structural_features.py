#!/usr/bin/env python3
"""
Compute structural descriptors for standardized protein structures.

The calculated descriptors are geometry-based proxies intended for comparative
screening across a large structural dataset.

Features:
    - residue/chain length
    - atom count
    - polar atom count
    - H-bond proxy count
    - salt-bridge proxy count
    - disulfide count
    - SASA total using Biopython ShrakeRupley, when available
    - per-residue normalised features
"""

from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time
import warnings

import numpy as np
import pandas as pd


try:
    from scipy.spatial import cKDTree
except Exception:
    cKDTree = None


try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley
    BIOPYTHON_SASA_AVAILABLE = True
except Exception:
    PDBParser = None
    MMCIFParser = None
    ShrakeRupley = None
    BIOPYTHON_SASA_AVAILABLE = False


AA3_STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

ACIDIC_ATOMS = {
    ("ASP", "OD1"), ("ASP", "OD2"),
    ("GLU", "OE1"), ("GLU", "OE2"),
}

BASIC_ATOMS = {
    ("LYS", "NZ"),
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("HIS", "ND1"), ("HIS", "NE2"),
}

POLAR_ELEMENTS = {"N", "O", "S"}


def load_yaml_config(path: str | None) -> dict:
    if not path:
        return {}

    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml") from exc

    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def deep_get(data: dict, keys: list[str], default=None):
    cur = data
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def resolve_path(path_value: str | None, project_root: str | None = None) -> Path | None:
    if not path_value:
        return None

    path = Path(os.path.expanduser(path_value))

    if path.is_absolute():
        return path

    if project_root:
        return Path(os.path.expanduser(project_root)) / path

    return Path.cwd() / path


class Logger:
    def __init__(self, log_path: Path | None):
        self.log_path = log_path
        if self.log_path:
            self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, msg: str) -> None:
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        if self.log_path:
            with self.log_path.open("a", encoding="utf-8") as handle:
                handle.write(line + "\n")


def clean_text(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "na", "n/a", "null"}:
        return ""
    return text



def extract_target_chain(row: pd.Series) -> str:
    """
    Extract selected chain from PDB-style structure_id.

    Examples:
        1AXK_B_Repair -> B
        7NWN_AAA_Repair -> AAA
        3NIY_A_Repair -> A

    MODELLER rows are not chain-filtered.
    """
    source = clean_text(row.get("structure_source_normalized", row.get("structure_source", ""))).lower()

    if source != "pdb":
        return ""

    structure_id = clean_text(row.get("structure_id", ""))

    if not structure_id:
        return ""

    core = structure_id
    if core.endswith("_Repair"):
        core = core[:-7]

    # Expected pattern: PDBID_CHAIN
    parts = core.split("_")

    if len(parts) >= 2 and len(parts[0]) == 4:
        return parts[1]

    return ""



def select_longest_standard_chain(atoms: list[dict]) -> tuple[str, list[dict]]:
    """
    Select the chain with the largest number of standard amino-acid residues.
    Used only as a fallback for PDB-derived structures without encoded chain ID.
    """
    chain_to_residues = {}

    for atom in atoms:
        if atom["record"] != "ATOM" or not atom["is_standard_aa"]:
            continue

        chain = atom["chain_id"]
        chain_to_residues.setdefault(chain, set()).add(atom["residue_key"])

    if not chain_to_residues:
        return "", atoms

    best_chain = sorted(
        chain_to_residues,
        key=lambda c: (-len(chain_to_residues[c]), c)
    )[0]

    return best_chain, [a for a in atoms if a["chain_id"] == best_chain]


def parse_float_slice(line: str, start: int, end: int):
    try:
        return float(line[start:end])
    except Exception:
        return None


def infer_element(line: str, atom_name: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element

    atom = atom_name.strip().upper()
    atom = atom.lstrip("0123456789")
    return atom[0] if atom else ""



def is_probably_mmcif(path: Path) -> bool:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for _ in range(20):
                line = handle.readline()
                if not line:
                    break
                stripped = line.strip()
                if stripped.startswith("data_") or stripped.startswith("_atom_site."):
                    return True
                if stripped.startswith("ATOM") or stripped.startswith("HETATM"):
                    return False
    except Exception:
        return False

    return path.suffix.lower() == ".cif"


def parse_atoms_with_biopython(path: Path) -> tuple[list[dict], str]:
    if not BIOPYTHON_SASA_AVAILABLE:
        return [], "Biopython parser unavailable"

    try:
        parser = MMCIFParser(QUIET=True) if is_probably_mmcif(path) else PDBParser(QUIET=True)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure(path.stem, str(path))

        atoms = []

        # Use only the first model to avoid inflated counts in multi-model structures.
        models = list(structure.get_models())
        if not models:
            return [], "No model found in structure"

        first_model = models[0]

        for atom in first_model.get_atoms():
            residue = atom.get_parent()
            chain = residue.get_parent()

            hetfield, resseq, icode = residue.id
            resname = residue.get_resname().strip().upper()
            atom_name = atom.get_name().strip().upper()
            element = clean_text(getattr(atom, "element", "")).upper()

            if not element:
                element = atom_name.lstrip("0123456789")[0] if atom_name else ""

            record = "ATOM" if hetfield == " " else "HETATM"

            coord = atom.get_coord()

            atoms.append({
                "record": record,
                "atom_name": atom_name,
                "resname": resname,
                "chain_id": clean_text(chain.id) or "_",
                "resseq": str(resseq),
                "icode": clean_text(icode),
                "residue_key": (clean_text(chain.id) or "_", str(resseq), clean_text(icode), resname),
                "coord": (float(coord[0]), float(coord[1]), float(coord[2])),
                "element": element,
                "is_standard_aa": resname in AA3_STANDARD,
            })

        if atoms:
            return atoms, ""

        return parse_atoms_with_biopython(path)

    except Exception as exc:
        return [], str(exc)


def parse_pdb_atoms(path: Path) -> tuple[list[dict], str]:
    if path.suffix.lower() == ".cif" or is_probably_mmcif(path):
        return parse_atoms_with_biopython(path)

    atoms = []
    current_model = 0
    seen_model_record = False

    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("MODEL"):
                    seen_model_record = True
                    current_model += 1
                    continue

                if seen_model_record and current_model > 1:
                    continue

                if not line.startswith(("ATOM", "HETATM")):
                    continue

                record = line[0:6].strip()
                atom_name = line[12:16].strip()
                altloc = line[16:17].strip()
                resname = line[17:20].strip().upper()
                chain_id = line[21:22].strip() or "_"
                resseq = line[22:26].strip()
                icode = line[26:27].strip()
                x = parse_float_slice(line, 30, 38)
                y = parse_float_slice(line, 38, 46)
                z = parse_float_slice(line, 46, 54)

                if x is None or y is None or z is None:
                    continue

                # Keep blank or A conformer; skip alternate B/C/etc. to avoid double counting.
                if altloc not in {"", "A"}:
                    continue

                element = infer_element(line, atom_name)

                atoms.append({
                    "record": record,
                    "atom_name": atom_name.upper(),
                    "resname": resname,
                    "chain_id": chain_id,
                    "resseq": resseq,
                    "icode": icode,
                    "residue_key": (chain_id, resseq, icode, resname),
                    "coord": (x, y, z),
                    "element": element,
                    "is_standard_aa": resname in AA3_STANDARD,
                })

        return atoms, ""

    except Exception as exc:
        return [], str(exc)


def count_hbond_proxy(atoms: list[dict], cutoff: float = 3.5) -> int:
    polar_indices = [
        i for i, atom in enumerate(atoms)
        if atom["record"] == "ATOM" and atom["element"] in POLAR_ELEMENTS
    ]

    if len(polar_indices) < 2:
        return 0

    coords = np.array([atoms[i]["coord"] for i in polar_indices], dtype=float)

    count = 0

    if cKDTree is not None:
        tree = cKDTree(coords)
        pairs = tree.query_pairs(cutoff)

        for i, j in pairs:
            ai = atoms[polar_indices[i]]
            aj = atoms[polar_indices[j]]

            if ai["residue_key"] == aj["residue_key"]:
                continue

            count += 1

        return int(count)

    # fallback without scipy; chunked but slower
    for i in range(len(polar_indices)):
        ai = atoms[polar_indices[i]]
        ci = np.array(ai["coord"], dtype=float)

        for j in range(i + 1, len(polar_indices)):
            aj = atoms[polar_indices[j]]

            if ai["residue_key"] == aj["residue_key"]:
                continue

            cj = np.array(aj["coord"], dtype=float)

            if np.linalg.norm(ci - cj) <= cutoff:
                count += 1

    return int(count)


def count_salt_bridges(atoms: list[dict], cutoff: float = 4.0) -> int:
    acidic = [
        atom for atom in atoms
        if atom["record"] == "ATOM" and (atom["resname"], atom["atom_name"]) in ACIDIC_ATOMS
    ]

    basic = [
        atom for atom in atoms
        if atom["record"] == "ATOM" and (atom["resname"], atom["atom_name"]) in BASIC_ATOMS
    ]

    if not acidic or not basic:
        return 0

    acid_coords = np.array([a["coord"] for a in acidic], dtype=float)
    basic_coords = np.array([b["coord"] for b in basic], dtype=float)

    count = 0

    if cKDTree is not None:
        acid_tree = cKDTree(acid_coords)
        basic_tree = cKDTree(basic_coords)
        hits = acid_tree.query_ball_tree(basic_tree, cutoff)

        for acid_idx, basic_hits in enumerate(hits):
            acid_res = acidic[acid_idx]["residue_key"]

            for basic_idx in basic_hits:
                if acid_res == basic[basic_idx]["residue_key"]:
                    continue
                count += 1

        return int(count)

    for a in acidic:
        ca = np.array(a["coord"], dtype=float)

        for b in basic:
            if a["residue_key"] == b["residue_key"]:
                continue

            cb = np.array(b["coord"], dtype=float)

            if np.linalg.norm(ca - cb) <= cutoff:
                count += 1

    return int(count)


def count_disulfides(atoms: list[dict], cutoff: float = 2.5) -> int:
    sg_atoms = [
        atom for atom in atoms
        if atom["record"] == "ATOM"
        and atom["resname"] == "CYS"
        and atom["atom_name"] == "SG"
    ]

    if len(sg_atoms) < 2:
        return 0

    coords = np.array([a["coord"] for a in sg_atoms], dtype=float)
    count = 0

    for i, j in itertools.combinations(range(len(sg_atoms)), 2):
        if sg_atoms[i]["residue_key"] == sg_atoms[j]["residue_key"]:
            continue

        dist = np.linalg.norm(coords[i] - coords[j])

        if dist <= cutoff:
            count += 1

    return int(count)


def compute_sasa(path: Path) -> tuple[float | None, str, str]:
    if not BIOPYTHON_SASA_AVAILABLE:
        return None, "not_available", "Biopython ShrakeRupley unavailable"

    try:
        parser = MMCIFParser(QUIET=True) if is_probably_mmcif(path) else PDBParser(QUIET=True)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure(path.stem, str(path))
            sr = ShrakeRupley()
            sr.compute(structure, level="A")

        total = 0.0

        for atom in structure.get_atoms():
            sasa = getattr(atom, "sasa", None)
            if sasa is not None and not math.isnan(float(sasa)):
                total += float(sasa)

        return total, "computed", ""

    except Exception as exc:
        return None, "failed", str(exc)


def compute_features_for_structure(row: pd.Series, compute_sasa_flag: bool) -> dict:
    path_text = clean_text(row.get("stage3_standardized_structure_path", "")) or clean_text(row.get("usable_structure_path", ""))
    path = Path(path_text)

    out = row.to_dict()
    out.update({
        "feature_status": "",
        "feature_error": "",
        "structure_path_for_features": path_text,
        "atom_count": 0,
        "hetatm_count": 0,
        "standard_residue_count": 0,
        "all_residue_count": 0,
        "chain_count": 0,
        "chain_ids": "",
        "polar_atom_count": 0,
        "hbond_proxy_count": 0,
        "salt_bridge_count": 0,
        "disulfide_count": 0,
        "sasa_total": "",
        "sasa_status": "not_requested" if not compute_sasa_flag else "",
        "sasa_error": "",
        "hbond_proxy_per_residue": "",
        "salt_bridge_per_residue": "",
        "disulfide_per_residue": "",
        "sasa_per_residue": "",
    })

    if not path_text:
        out["feature_status"] = "failed_missing_path"
        out["feature_error"] = "No structure path supplied"
        return out

    if not path.exists():
        out["feature_status"] = "failed_missing_file"
        out["feature_error"] = f"Structure file not found: {path}"
        return out

    atoms, parse_error = parse_pdb_atoms(path)

    if parse_error:
        out["feature_status"] = "failed_parse_error"
        out["feature_error"] = parse_error
        return out

    if not atoms:
        out["feature_status"] = "failed_no_atoms"
        out["feature_error"] = "No ATOM/HETATM records parsed"
        return out

    target_chain = extract_target_chain(row)
    source_norm = clean_text(row.get("structure_source_normalized", row.get("structure_source", ""))).lower()

    out["target_chain"] = target_chain
    out["chain_filter_status"] = "not_required"

    if target_chain:
        exact_chain_atoms = [a for a in atoms if a["chain_id"] == target_chain]

        if exact_chain_atoms:
            atoms = exact_chain_atoms
            out["chain_filter_status"] = "exact_chain_selected"
        elif len(target_chain) > 1:
            # Some parsers reduce multi-character mmCIF chain IDs to the first character.
            first_char_chain_atoms = [a for a in atoms if a["chain_id"] == target_chain[0]]

            if first_char_chain_atoms:
                atoms = first_char_chain_atoms
                out["chain_filter_status"] = "first_character_chain_selected"
            else:
                available_chains = sorted({a["chain_id"] for a in atoms})
                out["feature_status"] = "failed_target_chain_not_found"
                out["feature_error"] = f"Target chain {target_chain} not found; available chains: {available_chains}"
                return out
        else:
            available_chains = sorted({a["chain_id"] for a in atoms})
            out["feature_status"] = "failed_target_chain_not_found"
            out["feature_error"] = f"Target chain {target_chain} not found; available chains: {available_chains}"
            return out

    elif source_norm == "pdb":
        selected_chain, selected_atoms = select_longest_standard_chain(atoms)

        if selected_chain:
            atoms = selected_atoms
            out["target_chain"] = selected_chain
            out["chain_filter_status"] = "longest_chain_selected_no_target"
        else:
            out["feature_status"] = "failed_no_standard_chain_for_pdb"
            out["feature_error"] = "PDB-derived structure had no target chain and no standard amino-acid chain could be selected"
            return out

    atom_records = [a for a in atoms if a["record"] == "ATOM"]
    hetatm_records = [a for a in atoms if a["record"] == "HETATM"]

    standard_residues = {
        a["residue_key"]
        for a in atom_records
        if a["is_standard_aa"]
    }

    all_residues = {
        a["residue_key"]
        for a in atom_records
    }

    chains = sorted({a["chain_id"] for a in atom_records})

    residue_count = len(standard_residues)

    hbond_count = count_hbond_proxy(atoms)
    salt_count = count_salt_bridges(atoms)
    disulfide_count = count_disulfides(atoms)

    polar_atom_count = sum(
        1 for a in atom_records
        if a["element"] in POLAR_ELEMENTS
    )

    out["atom_count"] = int(len(atom_records))
    out["hetatm_count"] = int(len(hetatm_records))
    out["standard_residue_count"] = int(residue_count)
    out["all_residue_count"] = int(len(all_residues))
    out["chain_count"] = int(len(chains))
    out["chain_ids"] = ";".join(chains)
    out["polar_atom_count"] = int(polar_atom_count)
    out["hbond_proxy_count"] = int(hbond_count)
    out["salt_bridge_count"] = int(salt_count)
    out["disulfide_count"] = int(disulfide_count)

    if residue_count > 0:
        out["hbond_proxy_per_residue"] = float(hbond_count / residue_count)
        out["salt_bridge_per_residue"] = float(salt_count / residue_count)
        out["disulfide_per_residue"] = float(disulfide_count / residue_count)

    if compute_sasa_flag:
        sasa_total, sasa_status, sasa_error = compute_sasa(path)
        out["sasa_status"] = sasa_status
        out["sasa_error"] = sasa_error

        if sasa_total is not None:
            out["sasa_total"] = float(sasa_total)
            if residue_count > 0:
                out["sasa_per_residue"] = float(sasa_total / residue_count)

    out["feature_status"] = "computed"

    return out


def build_group_summary(features: pd.DataFrame) -> pd.DataFrame:
    computed = features[features["feature_status"] == "computed"].copy()

    if computed.empty:
        return pd.DataFrame(columns=[
            "organism_type",
            "gh_family",
            "structure_source_normalized",
            "n",
        ])

    numeric_cols = [
        "atom_count",
        "hetatm_count",
        "standard_residue_count",
        "all_residue_count",
        "chain_count",
        "polar_atom_count",
        "hbond_proxy_count",
        "salt_bridge_count",
        "disulfide_count",
        "sasa_total",
        "hbond_proxy_per_residue",
        "salt_bridge_per_residue",
        "disulfide_per_residue",
        "sasa_per_residue",
    ]

    for col in numeric_cols:
        if col in computed.columns:
            computed[col] = pd.to_numeric(computed[col], errors="coerce")

    group_cols = ["organism_type", "gh_family", "structure_source_normalized"]

    summary = (
        computed
        .groupby(group_cols, dropna=False)
        .agg(
            n=("uniprot_accession", "size"),
            unique_proteins=("uniprot_accession", "nunique"),
            mean_residue_count=("standard_residue_count", "mean"),
            median_residue_count=("standard_residue_count", "median"),
            mean_hbond_proxy_per_residue=("hbond_proxy_per_residue", "mean"),
            mean_salt_bridge_per_residue=("salt_bridge_per_residue", "mean"),
            mean_disulfide_per_residue=("disulfide_per_residue", "mean"),
            mean_sasa_per_residue=("sasa_per_residue", "mean"),
        )
        .reset_index()
    )

    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute structural descriptors for standardized structures.")

    parser.add_argument("--config", default=None, help="Optional YAML config.")
    parser.add_argument("--manifest", required=True, help="Standardized structure manifest CSV.")
    parser.add_argument("--output", required=True, help="Output structural feature CSV.")
    parser.add_argument("--group-summary-output", required=True, help="Output group summary CSV.")
    parser.add_argument("--problem-output", required=True, help="Output failed/warning rows CSV.")
    parser.add_argument("--no-sasa", action="store_true", help="Disable SASA calculation.")
    parser.add_argument("--log", default=None, help="Log file path.")

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_yaml_config(args.config)

    project_root = deep_get(config, ["project", "root"], None)

    manifest_path = resolve_path(args.manifest, project_root)
    output_path = resolve_path(args.output, project_root)
    group_summary_path = resolve_path(args.group_summary_output, project_root)
    problem_output_path = resolve_path(args.problem_output, project_root)

    log_path = resolve_path(
        args.log
        or str(Path(deep_get(config, ["outputs", "logs_dir"], "logs")) / "04_structure_03_features.log"),
        project_root,
    )

    logger = Logger(log_path)

    logger.write("Starting structural feature calculation")
    logger.write(f"Input manifest: {manifest_path}")
    logger.write(f"Output: {output_path}")
    logger.write(f"SASA enabled: {not args.no_sasa}")
    logger.write(f"SciPy cKDTree available: {cKDTree is not None}")
    logger.write(f"Biopython ShrakeRupley available: {BIOPYTHON_SASA_AVAILABLE}")

    if manifest_path is None or not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False).fillna("")

    rows = []
    problem_rows = []

    for i, (_, row) in enumerate(manifest.iterrows(), start=1):
        result = compute_features_for_structure(row, compute_sasa_flag=not args.no_sasa)
        rows.append(result)

        status = result.get("feature_status", "")

        if status != "computed":
            problem_rows.append(result)

        if i % 50 == 0 or i == len(manifest):
            logger.write(f"Processed {i}/{len(manifest)} structures")

    features = pd.DataFrame(rows)
    problems = pd.DataFrame(problem_rows, columns=features.columns)
    group_summary = build_group_summary(features)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    group_summary_path.parent.mkdir(parents=True, exist_ok=True)
    problem_output_path.parent.mkdir(parents=True, exist_ok=True)

    features.to_csv(output_path, index=False)
    group_summary.to_csv(group_summary_path, index=False)
    problems.to_csv(problem_output_path, index=False)

    status_counts = dict(Counter(features["feature_status"].astype(str)))
    sasa_counts = dict(Counter(features["sasa_status"].astype(str))) if "sasa_status" in features.columns else {}

    metadata = {
        "script": Path(__file__).name,
        "manifest": str(manifest_path),
        "output": str(output_path),
        "group_summary_output": str(group_summary_path),
        "problem_output": str(problem_output_path),
        "input_rows": int(len(manifest)),
        "output_rows": int(len(features)),
        "problem_rows": int(len(problems)),
        "status_counts": status_counts,
        "sasa_status_counts": sasa_counts,
        "scipy_ckdtree_available": cKDTree is not None,
        "biopython_sasa_available": BIOPYTHON_SASA_AVAILABLE,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    logger.write(f"Input rows: {len(manifest)}")
    logger.write(f"Output rows: {len(features)}")
    logger.write(f"Problem rows: {len(problems)}")
    logger.write(f"Feature status counts: {status_counts}")
    logger.write(f"SASA status counts: {sasa_counts}")
    logger.write("Finished structural feature calculation")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
