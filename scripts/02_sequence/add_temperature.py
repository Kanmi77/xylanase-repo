#!/usr/bin/env python3
# Purpose: Add temperature metadata.

from pathlib import Path
import re

import pandas as pd


BASE_DIR = Path.home() / "xylanase-thesis"

MASTER_FILE = (
    BASE_DIR
    / "data/thermal_sequence_rerun/master_692_thermal_annotated/"
    / "master_692_thermal_annotated_blast_expanded.csv"
)

BRENDA_FILES = [
    (
        BASE_DIR / "data/brenda/xylanase_temperature_optimum_FIXED.csv",
        "brenda_temperature_optimum",
        1,
    ),
    (
        BASE_DIR / "data/brenda/xylanase_temperature_stability_FIXED.csv",
        "brenda_temperature_stability",
        2,
    ),
    (
        BASE_DIR / "data/brenda/xylanase_temperature_range_FIXED.csv",
        "brenda_temperature_range",
        3,
    ),
]

OUTPUT_DIR = BASE_DIR / "data/thermal_sequence_rerun/master_692_thermal_annotated"
FASTA_DIR = OUTPUT_DIR / "fasta"
SUMMARY_DIR = BASE_DIR / "results/sequence_thermal_rerun/master_692_thermal_annotated/with_brenda"

OUTPUT_TABLE = OUTPUT_DIR / "master_692_thermal_annotated_blast_expanded_with_brenda.csv"
OUTPUT_FASTA = FASTA_DIR / "master_692_thermal_annotated_blast_expanded_with_brenda_all.fasta"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FASTA_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)


def normalise_id(value):
    if pd.isna(value):
        return ""

    text = str(value).strip()
    text = text.replace("\u00a0", "")
    text = text.replace(";", "")
    text = re.sub(r"\s+", "", text)

    return text.upper()


def extract_temperatures(value):
    if pd.isna(value):
        return []

    numbers = re.findall(r"[-+]?\d+(?:\.\d+)?", str(value))

    temperatures = []
    for number in numbers:
        temp = float(number)

        # Keep realistic enzyme temperature values.
        if 0 <= temp <= 120:
            temperatures.append(temp)

    return temperatures


def classify_temperature(temp):
    if pd.isna(temp):
        return "unknown"

    if float(temp) >= 60:
        return "thermophilic"

    return "mesophilic_lower"


def load_brenda_records():
    rows = []

    for path, basis, priority in BRENDA_FILES:
        if not path.exists():
            print(f"Missing BRENDA file, skipped: {path}")
            continue

        df = pd.read_csv(path)

        accession_col = None
        temperature_col = None

        for col in df.columns:
            if str(col).upper() in ["UNIPROT", "UNIPROT_ACCESSION", "ACCESSION"]:
                accession_col = col

            if str(col).upper() in ["TEMPERATURE", "TEMP", "TEMPERATURE_C"]:
                temperature_col = col

        if accession_col is None or temperature_col is None:
            print(f"Could not identify columns in: {path}")
            print(list(df.columns))
            continue

        for _, row in df.iterrows():
            accession = normalise_id(row[accession_col])
            temperatures = extract_temperatures(row[temperature_col])

            if not accession or not temperatures:
                continue

            rows.append(
                {
                    "uniprot_key": accession,
                    "brenda_basis": basis,
                    "brenda_priority": priority,
                    "brenda_temperature_c": max(temperatures),
                    "brenda_original_temperature": row[temperature_col],
                }
            )

    brenda = pd.DataFrame(rows)

    if brenda.empty:
        return brenda

    brenda["brenda_thermal_class"] = brenda["brenda_temperature_c"].apply(
        classify_temperature
    )

    # Keep the best BRENDA evidence per accession.
    # Priority order: optimum > stability > range.
    # Within the same basis, keep the highest reported temperature.
    brenda = (
        brenda.sort_values(
            ["uniprot_key", "brenda_priority", "brenda_temperature_c"],
            ascending=[True, True, False],
        )
        .drop_duplicates("uniprot_key")
        .copy()
    )

    return brenda


def make_fasta_header(row):
    temp = row["thermal_temperature_c"]

    if pd.isna(temp):
        temp_text = "NA"
    else:
        temp_text = f"{float(temp):.1f}C"

    return (
        f"{row['uniprot_accession']}|{row['organism_type']}|{row['gh_family']}|"
        f"thermal={row['thermal_class']}|temp={temp_text}"
    )


def write_fasta(df, path):
    with path.open("w") as handle:
        for _, row in df.iterrows():
            sequence = str(row["sequence_clean"]).strip().upper()

            if not sequence:
                continue

            handle.write(f">{make_fasta_header(row)}\n")
            handle.write(f"{sequence}\n")


def write_group_fastas(master):
    write_fasta(master, OUTPUT_FASTA)

    for organism_type in ["bacterial", "fungal", "archaea"]:
        for gh_family in ["GH10", "GH11"]:
            subset = master[
                (master["organism_type"] == organism_type)
                & (master["gh_family"] == gh_family)
            ].copy()

            if subset.empty:
                continue

            filename = (
                f"master_692_{organism_type}_{gh_family.lower()}"
                "_thermal_annotated_blast_expanded_with_brenda.fasta"
            )

            write_fasta(subset, FASTA_DIR / filename)


def main():
    master = pd.read_csv(MASTER_FILE)
    master = master.copy()

    master["uniprot_key"] = master["uniprot_accession"].apply(normalise_id)

    brenda = load_brenda_records()

    if brenda.empty:
        raise SystemExit("No usable BRENDA temperature records were found.")

    brenda.to_csv(SUMMARY_DIR / "brenda_temperature_records_best_per_accession.csv", index=False)

    master = master.merge(
        brenda[
            [
                "uniprot_key",
                "brenda_basis",
                "brenda_temperature_c",
                "brenda_thermal_class",
                "brenda_original_temperature",
            ]
        ],
        on="uniprot_key",
        how="left",
    )

    # Add BRENDA only to records still marked unknown.
    brenda_fill_mask = (
        (master["thermal_class"] == "unknown")
        & master["brenda_thermal_class"].isin(["thermophilic", "mesophilic_lower"])
    )

    master.loc[brenda_fill_mask, "thermal_class"] = master.loc[
        brenda_fill_mask,
        "brenda_thermal_class",
    ]

    master.loc[brenda_fill_mask, "thermal_temperature_c"] = master.loc[
        brenda_fill_mask,
        "brenda_temperature_c",
    ]

    master.loc[brenda_fill_mask, "thermal_label_source"] = master.loc[
        brenda_fill_mask,
        "brenda_basis",
    ]

    master.loc[brenda_fill_mask, "thermal_match_detail"] = (
        "brenda_raw_temperature_by_uniprot_accession"
    )

    master.to_csv(OUTPUT_TABLE, index=False)

    label_counts = (
        master["thermal_class"]
        .value_counts(dropna=False)
        .reset_index()
    )
    label_counts.columns = ["thermal_class", "count"]

    source_counts = (
        master["thermal_label_source"]
        .value_counts(dropna=False)
        .reset_index()
    )
    source_counts.columns = ["thermal_label_source", "count"]

    group_summary = (
        master.groupby(["organism_type", "gh_family", "thermal_class"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["organism_type", "gh_family", "thermal_class"])
    )

    newly_added = master[brenda_fill_mask].copy()

    label_counts.to_csv(SUMMARY_DIR / "master_692_with_brenda_label_counts.csv", index=False)
    source_counts.to_csv(SUMMARY_DIR / "master_692_with_brenda_label_source_counts.csv", index=False)
    group_summary.to_csv(SUMMARY_DIR / "master_692_with_brenda_group_summary.csv", index=False)
    newly_added.to_csv(SUMMARY_DIR / "new_labels_added_from_raw_brenda.csv", index=False)

    write_group_fastas(master)

    print("BRENDA raw temperature evidence added.")
    print(f"Rows in master: {len(master)}")
    print(f"New labels added from raw BRENDA: {len(newly_added)}")
    print(f"Output table: {OUTPUT_TABLE}")
    print(f"Output FASTA: {OUTPUT_FASTA}")

    print("\nThermal label counts:")
    print(label_counts.to_string(index=False))

    print("\nThermal label source counts:")
    print(source_counts.to_string(index=False))

    print("\nGroup summary:")
    print(group_summary.to_string(index=False))


if __name__ == "__main__":
    main()
