#!/usr/bin/env python3
"""
Match external Tmxyl/Toptxyl experimental temperature labels to the existing xylanase master dataset.

Matching levels:
1. Exact sequence-hash match
2. UniProt/accession-like match extracted from external protein name/accession fields

Outputs:
- external_master_exact_sequence_matches.csv
- external_master_accession_matches.csv
- external_master_all_matches_deduplicated.csv
- external_master_temperature_match_report.md
"""

from pathlib import Path
import re
import hashlib
import pandas as pd
import numpy as np


AA_RE = re.compile(r"[^A-Z]")


def clean_sequence(seq):
    if pd.isna(seq):
        return ""
    seq = str(seq).strip().upper()
    seq = AA_RE.sub("", seq)
    return seq


def sequence_hash(seq):
    return hashlib.md5(str(seq).encode()).hexdigest()


def extract_uniprot_like(text):
    """
    Extract possible UniProt-like accession from text.
    Handles examples like Q59962, P26514, B8XY24, I1RQU5.
    """
    if pd.isna(text):
        return ""

    text = str(text).replace("\xa0", " ").strip()

    patterns = [
        r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b",
        r"\b[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]\b",
        r"\b[A-Z][0-9][A-Z0-9]{3}[0-9]\b",
    ]

    for pat in patterns:
        m = re.search(pat, text)
        if m:
            return m.group(0)

    return ""


def parse_first_number(x):
    if pd.isna(x):
        return np.nan

    s = str(x)

    m = re.search(r"[-+]?\d+(?:\.\d+)?", s)

    if not m:
        return np.nan

    try:
        return float(m.group(0))
    except Exception:
        return np.nan


def find_master_file(root):
    candidates = [
        root / "data/curated/xylanase_master_deduplicated.csv",
        root / "data/curated/xylanase_master_all_curated_with_brenda.csv",
        root / "data/curated/xylanase_master_all_curated.csv",
        root / "data/curated/xylanase_master_all_curated.txt",
        root / "xylanase_master_all_curated.txt",
        root / "xylanase_master_deduplicated.csv",
    ]

    for p in candidates:
        if p.exists():
            return p

    # Fallback search.
    possible = list(root.glob("**/*xylanase*master*curated*.csv")) + list(root.glob("**/*xylanase*master*deduplicated*.csv"))

    if possible:
        return possible[0]

    raise FileNotFoundError("Could not find master xylanase dataset automatically.")


def read_table(path):
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)

    if path.suffix.lower() == ".txt":
        # Try tab first, then comma.
        try:
            return pd.read_csv(path, sep="\t", low_memory=False)
        except Exception:
            return pd.read_csv(path, low_memory=False)

    raise ValueError(f"Unsupported file type: {path}")


def detect_sequence_column(df):
    candidates = [
        "sequence",
        "Sequence",
        "protein_sequence",
        "aa_sequence",
        "amino_acid_sequence",
        "seq",
    ]

    for c in candidates:
        if c in df.columns:
            return c

    # fallback: longest string-like column with many amino acid-looking values
    best_col = None
    best_score = -1

    for c in df.columns:
        vals = df[c].dropna().astype(str).head(100)
        if vals.empty:
            continue

        lengths = vals.str.len()
        score = lengths.median()

        if score > best_score and score > 50:
            best_col = c
            best_score = score

    if best_col is None:
        raise RuntimeError("Could not detect sequence column in master dataset.")

    return best_col


def detect_accession_columns(df):
    cols = []

    for c in df.columns:
        low = c.lower()

        if "uniprot" in low or "accession" in low or low in ["entry", "entry_name", "id"]:
            cols.append(c)

    return cols


def detect_temperature_columns(df):
    cols = []

    keywords = [
        "tm",
        "topt",
        "temp",
        "temperature",
        "thermal",
        "brenda",
        "stability",
    ]

    for c in df.columns:
        low = c.lower()

        if any(k in low for k in keywords):
            cols.append(c)

    return cols


def main():
    root = Path(".").resolve()

    external_file = root / "results/experimental_ml_correction/tmxyl_toptxyl_sequence_deduplicated.csv"

    if not external_file.exists():
        raise FileNotFoundError(external_file)

    master_file = find_master_file(root)

    outdir = root / "results/experimental_ml_correction/master_temperature_matching"
    outdir.mkdir(parents=True, exist_ok=True)

    external = pd.read_csv(external_file, low_memory=False)
    master = read_table(master_file)

    sequence_col = detect_sequence_column(master)
    accession_cols = detect_accession_columns(master)
    temp_cols = detect_temperature_columns(master)

    external["sequence_clean"] = external["sequence"].apply(clean_sequence)
    external["sequence_hash_match"] = external["sequence_clean"].apply(sequence_hash)

    external["external_accession_like"] = external["example_protein_name_or_accession"].apply(extract_uniprot_like)

    master["master_sequence_clean"] = master[sequence_col].apply(clean_sequence)
    master["sequence_hash_match"] = master["master_sequence_clean"].apply(sequence_hash)

    # Build master accession field.
    master_accession_records = []

    for idx, row in master.iterrows():
        vals = []

        for c in accession_cols:
            v = row.get(c, "")
            if pd.isna(v):
                continue

            v = str(v).replace("\xa0", " ").strip()

            extracted = extract_uniprot_like(v)

            if extracted:
                vals.append(extracted)
            elif v:
                vals.append(v)

        vals = sorted(set(vals))

        master_accession_records.append(";".join(vals))

    master["master_accession_like"] = master_accession_records

    # Exact sequence matches.
    seq_matches = external.merge(
        master,
        on="sequence_hash_match",
        how="inner",
        suffixes=("_external", "_master"),
    )

    seq_matches["match_type"] = "exact_sequence"

    # Accession matches.
    external_acc = external[external["external_accession_like"].astype(str).str.len() > 0].copy()
    master_acc = master[master["master_accession_like"].astype(str).str.len() > 0].copy()

    # explode master accession list
    master_acc = master_acc.assign(master_accession_like=master_acc["master_accession_like"].str.split(";")).explode("master_accession_like")
    master_acc["master_accession_like"] = master_acc["master_accession_like"].astype(str).str.strip()

    acc_matches = external_acc.merge(
        master_acc,
        left_on="external_accession_like",
        right_on="master_accession_like",
        how="inner",
        suffixes=("_external", "_master"),
    )

    acc_matches["match_type"] = "accession_like"

    # Save raw matches.
    seq_matches.to_csv(outdir / "external_master_exact_sequence_matches.csv", index=False)
    acc_matches.to_csv(outdir / "external_master_accession_matches.csv", index=False)

    all_matches = pd.concat([seq_matches, acc_matches], ignore_index=True, sort=False)

    if len(all_matches) > 0:
        # Make a conservative deduplicated match key.
        key_cols = [
            "sequence_hash_external" if "sequence_hash_external" in all_matches.columns else None,
            "sequence_hash_match",
            "target_type",
            "temperature_median_c",
        ]

        key_cols = [c for c in key_cols if c in all_matches.columns]

        all_matches = all_matches.drop_duplicates(subset=key_cols + ["match_type"], keep="first")

    all_matches.to_csv(outdir / "external_master_all_matches_deduplicated.csv", index=False)

    # Temperature comparison.
    temp_comparison_rows = []

    if len(all_matches) > 0 and temp_cols:
        for c in temp_cols:
            if c not in all_matches.columns:
                # after suffixing, master columns may retain original names unless conflict
                possible_cols = [x for x in all_matches.columns if x == c or x.endswith("_master") and x.replace("_master", "") == c]
            else:
                possible_cols = [c]

            for pc in possible_cols:
                if pc not in all_matches.columns:
                    continue

                ext_temp = pd.to_numeric(all_matches["temperature_median_c"], errors="coerce")
                master_temp = all_matches[pc].apply(parse_first_number)

                valid = ext_temp.notna() & master_temp.notna()

                if valid.sum() == 0:
                    continue

                diff = (ext_temp[valid] - master_temp[valid]).abs()

                temp_comparison_rows.append({
                    "master_temperature_column": pc,
                    "n_matches_with_both_values": int(valid.sum()),
                    "n_close_within_1c": int((diff <= 1.0).sum()),
                    "n_close_within_2c": int((diff <= 2.0).sum()),
                    "mean_abs_difference_c": float(diff.mean()),
                    "median_abs_difference_c": float(diff.median()),
                    "min_abs_difference_c": float(diff.min()),
                    "max_abs_difference_c": float(diff.max()),
                })

    temp_comparison = pd.DataFrame(temp_comparison_rows)
    temp_comparison.to_csv(outdir / "external_master_temperature_value_comparison.csv", index=False)

    # Compact match table.
    compact_cols = []

    for c in [
        "target_type",
        "temperature_median_c",
        "thermal_label_60c_median",
        "example_protein_name_or_accession",
        "external_accession_like",
        "sequence_length_external",
        "match_type",
        "master_accession_like",
        sequence_col,
    ]:
        if c in all_matches.columns:
            compact_cols.append(c)

    # include likely master identifiers
    for c in accession_cols:
        if c in all_matches.columns and c not in compact_cols:
            compact_cols.append(c)

    compact = all_matches[compact_cols].copy() if len(all_matches) > 0 and compact_cols else pd.DataFrame()
    compact.to_csv(outdir / "external_master_matches_compact.csv", index=False)

    # Report.
    report = outdir / "EXTERNAL_TEMPERATURE_MASTER_MATCH_REPORT.md"

    with report.open("w") as fh:
        fh.write("# External Tmxyl/Toptxyl vs master temperature match report\n\n")

        fh.write("## Files\n\n")
        fh.write(f"- External experimental label file: `{external_file}`\n")
        fh.write(f"- Master dataset used: `{master_file}`\n\n")

        fh.write("## Master dataset detection\n\n")
        fh.write(f"- Master rows: {len(master)}\n")
        fh.write(f"- Detected sequence column: `{sequence_col}`\n")
        fh.write(f"- Detected accession columns: {accession_cols}\n")
        fh.write(f"- Detected possible temperature columns: {temp_cols}\n\n")

        fh.write("## External dataset\n\n")
        fh.write(f"- External deduplicated rows: {len(external)}\n")
        fh.write(f"- Tm external rows: {(external['target_type'] == 'Tm').sum()}\n")
        fh.write(f"- Topt external rows: {(external['target_type'] == 'Topt').sum()}\n")
        fh.write(f"- External rows with accession-like identifier: {(external['external_accession_like'].astype(str).str.len() > 0).sum()}\n\n")

        fh.write("## Match counts\n\n")
        fh.write(f"- Exact sequence matches: {len(seq_matches)}\n")
        fh.write(f"- Accession-like matches: {len(acc_matches)}\n")
        fh.write(f"- Combined deduplicated matches: {len(all_matches)}\n\n")

        if len(all_matches) > 0:
            fh.write("### Matches by target type\n\n")
            fh.write(all_matches["target_type"].value_counts(dropna=False).to_string())
            fh.write("\n\n")

            fh.write("### Matches by match type\n\n")
            fh.write(all_matches["match_type"].value_counts(dropna=False).to_string())
            fh.write("\n\n")

        fh.write("## Temperature-value comparison\n\n")

        if temp_comparison.empty:
            fh.write(
                "No direct numeric temperature comparison could be made. "
                "This may mean that the master dataset does not contain directly comparable numeric Tm/Topt columns, "
                "or that the temperature columns use non-numeric text formats.\n\n"
            )
        else:
            fh.write(temp_comparison.to_string(index=False))
            fh.write("\n\n")

        fh.write("## Interpretation\n\n")
        fh.write(
            "Exact sequence matches are the strongest evidence that an external experimental Tm/Topt label can be attached "
            "to a protein in the master dataset. Accession-like matches are useful but should be checked, especially because "
            "some records may represent engineered variants under the same accession/name. Temperature comparisons should be "
            "treated cautiously if the master dataset contains BRENDA-derived ranges, text descriptions, or heterogeneous "
            "temperature metadata rather than single directly measured Tm/Topt values.\n"
        )

    print("[DONE] External experimental labels matched to master dataset.")
    print(f"Report: {report}")
    print(f"Compact matches: {outdir / 'external_master_matches_compact.csv'}")


if __name__ == "__main__":
    main()
