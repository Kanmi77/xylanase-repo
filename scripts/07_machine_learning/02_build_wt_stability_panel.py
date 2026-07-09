#!/usr/bin/env python3
"""
Build the wild-type experimental stability ML panel.

Goal:
    Estimate wild-type xylanase stability from experimental Tm/Topt.

Inputs:
    - Deduplicated Tmxyl/Toptxyl experimental labels
    - Curated GH10/GH11 master xylanase table
    - PDB structural subset with FoldX WT stability

Output:
    results/stability_ml_goal/datasets/wt_experimental_stability_panel.csv
"""

from pathlib import Path
import hashlib
import math
import re

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis


BASE_DIR = Path.home() / "xylanase-thesis"

LABEL_FILE = (
    BASE_DIR
    / "results"
    / "experimental_ml_correction"
    / "tmxyl_toptxyl_sequence_deduplicated.csv"
)

MASTER_FILE = (
    BASE_DIR
    / "data"
    / "curated"
    / "xylanase_master_gh10_gh11_frozen_thesis_with_refseq_brenda.csv"
)

STRUCTURE_FILE = (
    BASE_DIR
    / "data"
    / "curated"
    / "xylanase_structured_subset_with_foldx_norm.csv"
)

OUTPUT_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "datasets"
    / "wt_experimental_stability_panel.csv"
)

REPORT_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "reports"
    / "wt_experimental_stability_panel_report.md"
)


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

AA_GROUPS = {
    "hydrophobic_fraction": set("AILMFWV"),
    "aromatic_fraction": set("FWY"),
    "polar_fraction": set("STNQCY"),
    "charged_fraction": set("DEKRH"),
    "acidic_fraction": set("DE"),
    "basic_fraction": set("KRH"),
    "tiny_fraction": set("ACGST"),
    "small_fraction": set("ACDGNPSTV"),
    "proline_fraction": set("P"),
    "glycine_fraction": set("G"),
    "cysteine_fraction": set("C"),
}


def clean_sequence(sequence):
    """Clean protein sequence."""
    if pd.isna(sequence):
        return ""

    sequence = str(sequence).upper()
    sequence = re.sub(r"[^A-Z]", "", sequence)

    return sequence


def sequence_md5(sequence):
    """Calculate MD5 hash for a cleaned sequence."""
    sequence = clean_sequence(sequence)
    return hashlib.md5(sequence.encode()).hexdigest()


def clean_accession(value):
    """Extract a likely UniProt accession from a text field."""
    if pd.isna(value):
        return ""

    text = str(value).strip().upper().replace("\xa0", "")

    patterns = [
        r"\b[A-NR-Z][0-9][A-Z0-9]{3}[0-9]\b",
        r"\b[A-Z0-9]{6,10}\b",
    ]

    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(0)

    return text


def safe_ratio(numerator, denominator):
    """Calculate a safe ratio."""
    if denominator == 0:
        return math.nan
    return numerator / denominator


def calculate_sequence_features(sequence):
    """Calculate sequence descriptors from a protein sequence."""
    sequence = clean_sequence(sequence)
    length = len(sequence)

    row = {
        "seq_length_calculated": length,
        "seq_invalid_residue_count": sum(
            1 for residue in sequence if residue not in STANDARD_AA
        ),
    }

    if length == 0:
        return row

    for aa in sorted(STANDARD_AA):
        row[f"aa_frac_{aa}"] = sequence.count(aa) / length

    for name, residues in AA_GROUPS.items():
        row[name] = sum(sequence.count(aa) for aa in residues) / length

    polar = row["polar_fraction"]
    basic = row["basic_fraction"]

    row["charged_to_polar_ratio"] = safe_ratio(row["charged_fraction"], polar)
    row["hydrophobic_to_polar_ratio"] = safe_ratio(
        row["hydrophobic_fraction"],
        polar,
    )
    row["acidic_to_basic_ratio"] = safe_ratio(row["acidic_fraction"], basic)

    if row["seq_invalid_residue_count"] == 0:
        analysis = ProteinAnalysis(sequence)

        row["molecular_weight"] = analysis.molecular_weight()
        row["isoelectric_point"] = analysis.isoelectric_point()
        row["aromaticity"] = analysis.aromaticity()
        row["instability_index"] = analysis.instability_index()
        row["gravy"] = analysis.gravy()

        helix, turn, sheet = analysis.secondary_structure_fraction()
        row["helix_fraction_protparam"] = helix
        row["turn_fraction_protparam"] = turn
        row["sheet_fraction_protparam"] = sheet
    else:
        row["molecular_weight"] = math.nan
        row["isoelectric_point"] = math.nan
        row["aromaticity"] = math.nan
        row["instability_index"] = math.nan
        row["gravy"] = math.nan
        row["helix_fraction_protparam"] = math.nan
        row["turn_fraction_protparam"] = math.nan
        row["sheet_fraction_protparam"] = math.nan

    return row


def aggregate_structure_features(structure):
    """Aggregate PDB/FoldX rows by accession."""
    numeric_columns = [
        "chain_length",
        "hbond_proxy_count",
        "salt_bridge_count",
        "disulfide_count",
        "sasa_total",
        "foldx_wt_total_energy",
        "foldx_energy_per_residue",
    ]

    rows = []

    for accession, group in structure.groupby("uniprot_accession_clean"):
        row = {
            "uniprot_accession_clean": accession,
            "has_pdb_structure": 1,
            "n_pdb_structures": len(group),
            "pdb_ids": ";".join(sorted(group["pdb_id"].dropna().astype(str).unique())),
            "pdb_tags": ";".join(sorted(group["pdb_tag"].dropna().astype(str).unique())),
        }

        for column in numeric_columns:
            if column in group.columns:
                row[f"{column}_mean"] = pd.to_numeric(
                    group[column],
                    errors="coerce",
                ).mean()
                row[f"{column}_min"] = pd.to_numeric(
                    group[column],
                    errors="coerce",
                ).min()
                row[f"{column}_max"] = pd.to_numeric(
                    group[column],
                    errors="coerce",
                ).max()

        rows.append(row)

    return pd.DataFrame(rows)


def write_report(panel):
    """Write markdown report."""
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Wild-Type Experimental Stability ML Panel",
        "",
        f"Rows: {len(panel)}",
        f"Rows with target temperature: {panel['target_temperature_c'].notna().sum()}",
        f"Rows matched to curated master: {panel['matched_master'].sum()}",
        f"Rows with PDB/FoldX structure: {panel['has_pdb_structure'].sum()}",
        "",
        "## Target summary",
        "",
        panel["target_temperature_c"].describe().to_string(),
        "",
        "## Thermal class counts",
        "",
        panel["thermal_label"].value_counts(dropna=False).to_string(),
        "",
        "## Target type counts",
        "",
        panel["target_type"].value_counts(dropna=False).to_string(),
        "",
        "## GH family counts after master matching",
        "",
        panel["gh_family"].value_counts(dropna=False).to_string(),
        "",
        "## Organism type counts after master matching",
        "",
        panel["organism_type"].value_counts(dropna=False).to_string(),
        "",
    ]

    REPORT_FILE.write_text("\n".join(lines) + "\n")


def main():
    """Build WT stability panel."""
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(LABEL_FILE)
    master = pd.read_csv(MASTER_FILE)
    structure = pd.read_csv(STRUCTURE_FILE)

    labels["sequence_clean"] = labels["sequence"].apply(clean_sequence)
    labels["sequence_hash_check"] = labels["sequence_clean"].apply(sequence_md5)
    labels["example_accession_clean"] = labels[
        "example_protein_name_or_accession"
    ].apply(clean_accession)

    master["sequence_clean_master"] = master["sequence"].apply(clean_sequence)
    master["sequence_hash_master"] = master["sequence_clean_master"].apply(
        sequence_md5
    )
    master["uniprot_accession_clean"] = master["uniprot_accession"].apply(
        clean_accession
    )

    structure["uniprot_accession_clean"] = structure[
        "uniprot_accession"
    ].apply(clean_accession)

    structure_summary = aggregate_structure_features(structure)

    sequence_feature_rows = [
        calculate_sequence_features(sequence)
        for sequence in labels["sequence_clean"]
    ]

    sequence_features = pd.DataFrame(sequence_feature_rows)

    panel = pd.concat(
        [
            labels.reset_index(drop=True),
            sequence_features.reset_index(drop=True),
        ],
        axis=1,
    )

    master_keep = [
        "sequence_hash_master",
        "uniprot_accession_clean",
        "organism_type",
        "organism",
        "gh_family",
        "xref_cazy",
        "has_cazy_xref",
        "cazy_status",
        "brenda_temperature_optimum",
        "brenda_temperature_range",
        "brenda_temperature_stability",
        "brenda_ph_optimum",
    ]

    master_by_sequence = master[master_keep].drop_duplicates(
        "sequence_hash_master"
    )

    panel = panel.merge(
        master_by_sequence,
        left_on="sequence_hash",
        right_on="sequence_hash_master",
        how="left",
    )

    panel["matched_master"] = panel["uniprot_accession_clean"].notna().astype(int)

    fallback = panel["matched_master"] == 0

    if fallback.any():
        master_by_accession = master.drop_duplicates("uniprot_accession_clean")
        master_by_accession = master_by_accession[
            [
                "uniprot_accession_clean",
                "organism_type",
                "organism",
                "gh_family",
                "xref_cazy",
                "has_cazy_xref",
                "cazy_status",
                "brenda_temperature_optimum",
                "brenda_temperature_range",
                "brenda_temperature_stability",
                "brenda_ph_optimum",
            ]
        ]

        fallback_panel = panel.loc[fallback].drop(
            columns=[
                "uniprot_accession_clean",
                "organism_type",
                "organism",
                "gh_family",
                "xref_cazy",
                "has_cazy_xref",
                "cazy_status",
                "brenda_temperature_optimum",
                "brenda_temperature_range",
                "brenda_temperature_stability",
                "brenda_ph_optimum",
            ],
            errors="ignore",
        )

        fallback_panel = fallback_panel.merge(
            master_by_accession,
            left_on="example_accession_clean",
            right_on="uniprot_accession_clean",
            how="left",
        )

        for column in fallback_panel.columns:
            panel.loc[fallback, column] = fallback_panel[column].values

        panel["matched_master"] = panel[
            "uniprot_accession_clean"
        ].notna().astype(int)

    panel = panel.merge(
        structure_summary,
        on="uniprot_accession_clean",
        how="left",
    )

    panel["has_pdb_structure"] = panel["has_pdb_structure"].fillna(0).astype(int)
    panel["n_pdb_structures"] = panel["n_pdb_structures"].fillna(0).astype(int)

    panel["target_temperature_c"] = panel["temperature_median_c"]
    panel["thermal_label"] = panel["thermal_label_60c_median"]

    panel.to_csv(OUTPUT_FILE, index=False)
    write_report(panel)

    print(f"Wrote: {OUTPUT_FILE}")
    print(f"Rows: {len(panel)}")
    print(f"Columns: {len(panel.columns)}")

    print("\nRows with target temperature:")
    print(panel["target_temperature_c"].notna().sum())

    print("\nMatched master:")
    print(panel["matched_master"].value_counts(dropna=False))

    print("\nPDB/FoldX structure:")
    print(panel["has_pdb_structure"].value_counts(dropna=False))

    print("\nThermal label counts:")
    print(panel["thermal_label"].value_counts(dropna=False))

    print("\nTarget temperature summary:")
    print(panel["target_temperature_c"].describe())


if __name__ == "__main__":
    main()
