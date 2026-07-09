#!/usr/bin/env python3
# Purpose: Prepare experimental thermal labels.

from pathlib import Path
import hashlib
import pandas as pd
import numpy as np
import re


AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(seq):
    if pd.isna(seq):
        return ""

    seq = str(seq).strip().upper()
    seq = re.sub(r"[^A-Z]", "", seq)

    return seq


def sequence_hash(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def invalid_residue_count(seq):
    return sum(1 for aa in seq if aa not in AA)


def binary_label(temp):
    if pd.isna(temp):
        return "unknown"

    return "thermophilic" if float(temp) >= 60.0 else "mesophilic_lower"


def main():
    root = Path(".").resolve()

    tm_file = root / "data/external/TmxylAndToptxyl/extracted/Tm/data_Tm.xlsx"
    topt_file = root / "data/external/TmxylAndToptxyl/extracted/Topt/data_Topt.xlsx"

    outdir = root / "results/experimental_ml_correction"
    outdir.mkdir(parents=True, exist_ok=True)

    if not tm_file.exists():
        raise FileNotFoundError(tm_file)

    if not topt_file.exists():
        raise FileNotFoundError(topt_file)

    tm_raw = pd.read_excel(tm_file)
    topt_raw = pd.read_excel(topt_file)

    # Standardize Tm dataset.
    tm = pd.DataFrame(index=tm_raw.index)
    tm["source_dataset"] = "Tmxyl"
    tm["target_type"] = "Tm"
    tm["original_no"] = tm_raw.get("NO.", pd.Series(range(1, len(tm_raw) + 1)))
    tm["paper_name"] = tm_raw.get("Paper Name", "")
    tm["pmid"] = tm_raw.get("PMID", "")
    tm["organism_or_biology"] = tm_raw.get("biology", "")
    tm["protein_name_or_accession"] = tm_raw.get("Name", "")
    tm["protein_family"] = np.nan
    tm["sequence"] = tm_raw["Sequence"].apply(clean_sequence)
    tm["temperature_c"] = pd.to_numeric(tm_raw["Tm"], errors="coerce")
    tm["thermal_label_60c"] = tm["temperature_c"].apply(binary_label)

    # Standardize Topt dataset.
    topt = pd.DataFrame(index=topt_raw.index)
    topt["source_dataset"] = "Toptxyl"
    topt["target_type"] = "Topt"
    topt["original_no"] = topt_raw.get("NO.", pd.Series(range(1, len(topt_raw) + 1)))
    topt["paper_name"] = np.nan
    topt["pmid"] = np.nan
    topt["organism_or_biology"] = np.nan
    topt["protein_name_or_accession"] = topt_raw.get("name", "")
    topt["protein_family"] = topt_raw.get("protein family", "")
    topt["sequence"] = topt_raw["Sequence"].apply(clean_sequence)
    topt["temperature_c"] = pd.to_numeric(topt_raw["Topt"], errors="coerce")
    topt["thermal_label_60c"] = topt["temperature_c"].apply(binary_label)

    combined = pd.concat([tm, topt], ignore_index=True)

    combined["sequence_length"] = combined["sequence"].str.len()
    combined["sequence_hash"] = combined["sequence"].apply(sequence_hash)
    combined["invalid_residue_count"] = combined["sequence"].apply(invalid_residue_count)
    combined["valid_standard_aa_only"] = combined["invalid_residue_count"] == 0

    # Save clean full datasets.
    tm.to_csv(outdir / "tmxyl_clean_experimental_labels.csv", index=False)
    topt.to_csv(outdir / "toptxyl_clean_experimental_labels.csv", index=False)
    combined.to_csv(outdir / "tmxyl_toptxyl_combined_experimental_labels.csv", index=False)

    # Deduplicate by exact sequence within target type.
    dedup_records = []

    group_cols = ["sequence_hash", "sequence", "target_type"]

    for keys, g in combined.groupby(group_cols, dropna=False):
        sequence_hash_value, seq, target_type = keys

        temps = g["temperature_c"].dropna().astype(float)

        labels = sorted(set(g["thermal_label_60c"].dropna()))

        if len(temps) == 0:
            temp_mean = np.nan
            temp_median = np.nan
            temp_min = np.nan
            temp_max = np.nan
        else:
            temp_mean = temps.mean()
            temp_median = temps.median()
            temp_min = temps.min()
            temp_max = temps.max()

        label_from_median = binary_label(temp_median)

        conflict = False
        if "thermophilic" in labels and "mesophilic_lower" in labels:
            conflict = True

        first = g.iloc[0]

        dedup_records.append({
            "sequence_hash": sequence_hash_value,
            "sequence": seq,
            "target_type": target_type,
            "source_datasets": ";".join(sorted({str(x) for x in g["source_dataset"].dropna()})),
            "n_measurements": len(g),
            "temperature_mean_c": temp_mean,
            "temperature_median_c": temp_median,
            "temperature_min_c": temp_min,
            "temperature_max_c": temp_max,
            "thermal_label_60c_median": label_from_median,
            "has_conflicting_60c_labels": conflict,
            "sequence_length": len(seq),
            "invalid_residue_count": invalid_residue_count(seq),
            "valid_standard_aa_only": invalid_residue_count(seq) == 0,
            "example_protein_name_or_accession": first.get("protein_name_or_accession", ""),
            "example_organism_or_biology": first.get("organism_or_biology", ""),
            "example_protein_family": first.get("protein_family", ""),
        })

    dedup = pd.DataFrame(dedup_records)

    dedup.to_csv(outdir / "tmxyl_toptxyl_sequence_deduplicated.csv", index=False)

    # Separate deduplicated Tm and Topt.
    dedup_tm = dedup[dedup["target_type"] == "Tm"].copy()
    dedup_topt = dedup[dedup["target_type"] == "Topt"].copy()

    dedup_tm.to_csv(outdir / "tmxyl_sequence_deduplicated.csv", index=False)
    dedup_topt.to_csv(outdir / "toptxyl_sequence_deduplicated.csv", index=False)

    # Summary.
    def label_counts(df, label_col):
        return df[label_col].value_counts(dropna=False).to_dict()

    report = outdir / "TMXYL_TOPTXYL_PREPARATION_REPORT.md"

    with report.open("w") as fh:
        fh.write("# Tmxyl/Toptxyl experimental label preparation report\n\n")

        fh.write("## Raw input datasets\n\n")
        fh.write(f"- Tmxyl raw rows: {len(tm_raw)}\n")
        fh.write(f"- Toptxyl raw rows: {len(topt_raw)}\n\n")

        fh.write("## Standardized combined dataset\n\n")
        fh.write(f"- Combined rows: {len(combined)}\n")
        fh.write(f"- Tm rows: {(combined['target_type'] == 'Tm').sum()}\n")
        fh.write(f"- Topt rows: {(combined['target_type'] == 'Topt').sum()}\n")
        fh.write(f"- Valid standard amino-acid sequences: {combined['valid_standard_aa_only'].sum()} / {len(combined)}\n")
        fh.write(f"- Label counts at 60C threshold: {label_counts(combined, 'thermal_label_60c')}\n\n")

        fh.write("## Exact sequence deduplication\n\n")
        fh.write(f"- Deduplicated sequence-target rows: {len(dedup)}\n")
        fh.write(f"- Deduplicated Tm rows: {len(dedup_tm)}\n")
        fh.write(f"- Deduplicated Topt rows: {len(dedup_topt)}\n")
        fh.write(f"- Tm label counts after median aggregation: {label_counts(dedup_tm, 'thermal_label_60c_median')}\n")
        fh.write(f"- Topt label counts after median aggregation: {label_counts(dedup_topt, 'thermal_label_60c_median')}\n")
        fh.write(f"- Conflicting sequence labels at 60C threshold: {dedup['has_conflicting_60c_labels'].sum()}\n\n")

        fh.write("## Output files\n\n")
        fh.write("- `tmxyl_clean_experimental_labels.csv`\n")
        fh.write("- `toptxyl_clean_experimental_labels.csv`\n")
        fh.write("- `tmxyl_toptxyl_combined_experimental_labels.csv`\n")
        fh.write("- `tmxyl_toptxyl_sequence_deduplicated.csv`\n")
        fh.write("- `tmxyl_sequence_deduplicated.csv`\n")
        fh.write("- `toptxyl_sequence_deduplicated.csv`\n\n")

        fh.write("## Interpretation\n\n")
        fh.write(
            "The prepared dataset provides experimental temperature labels for xylanase ML correction. "
            "The 60C threshold follows the thermophilic/mesophilic split used in the source study. "
            "Exact duplicate sequences were aggregated by median temperature to reduce repeated-measurement bias. "
            "Rows with conflicting labels should be reviewed before final ML training.\n"
        )

    print("[DONE] Prepared Tmxyl/Toptxyl experimental labels.")
    print(f"Report: {report}")
    print(f"Combined: {outdir / 'tmxyl_toptxyl_combined_experimental_labels.csv'}")
    print(f"Deduplicated: {outdir / 'tmxyl_toptxyl_sequence_deduplicated.csv'}")


if __name__ == "__main__":
    main()
