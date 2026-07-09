#!/usr/bin/env python3
"""
Merge external Tmxyl/Toptxyl labels with previous internal labels after
relabeling the internal dataset using the same 60C threshold.

Important:
- Old internal thermal_class is NOT used.
- Internal labels are recalculated from topt_numeric:
    topt_numeric >= 60 -> thermophilic
    topt_numeric < 60  -> mesophilic_lower
"""

from pathlib import Path
import re
import hashlib
import importlib.util
import numpy as np
import pandas as pd


AA_SET = set("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(seq):
    if pd.isna(seq):
        return ""
    seq = str(seq).strip().upper()
    seq = re.sub(r"[^A-Z]", "", seq)
    return seq


def sequence_hash(seq):
    return hashlib.md5(str(seq).encode()).hexdigest()


def invalid_residue_count(seq):
    return sum(1 for aa in str(seq) if aa not in AA_SET)


def safe_join(values, max_len=1000):
    cleaned = []
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).replace("\xa0", " ").strip()
        if s and s.lower() != "nan":
            cleaned.append(s)
    return ";".join(sorted(set(cleaned)))[:max_len]


def load_classifier_functions(root):
    script_path = root / "scripts/11_experimental_ml_correction/02_train_experimental_sequence_classifier.py"

    if not script_path.exists():
        raise FileNotFoundError(script_path)

    spec = importlib.util.spec_from_file_location("seqclf", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    return mod.extract_features, mod.evaluate_models


def main():
    root = Path(".").resolve()

    external_path = root / "results/experimental_ml_correction/tmxyl_toptxyl_sequence_deduplicated.csv"
    internal_path = root / "results/sequence_thermal_classifier/sequence_thermal_classifier_dataset.csv"
    master_path = root / "data/curated/xylanase_master_deduplicated.csv"

    outdir = root / "results/experimental_ml_correction/combined_60c_aligned_classifier"
    outdir.mkdir(parents=True, exist_ok=True)

    extract_features, evaluate_models = load_classifier_functions(root)

    external = pd.read_csv(external_path, low_memory=False)
    internal = pd.read_csv(internal_path, low_memory=False)
    master = pd.read_csv(master_path, low_memory=False)

    # -------------------------
    # External Tmxyl/Toptxyl
    # -------------------------
    external = external[
        (external["valid_standard_aa_only"] == True)
        & (external["thermal_label_60c_median"].isin(["thermophilic", "mesophilic_lower"]))
        & (external["has_conflicting_60c_labels"] == False)
    ].copy()

    ext_rows = pd.DataFrame()
    ext_rows["sequence"] = external["sequence"].apply(clean_sequence)
    ext_rows["thermal_label"] = external["thermal_label_60c_median"]
    ext_rows["temperature_c"] = pd.to_numeric(external["temperature_median_c"], errors="coerce")
    ext_rows["source_group"] = "external_Tmxyl_Toptxyl"
    ext_rows["target_type"] = external["target_type"].astype(str)
    ext_rows["source_detail"] = external["source_datasets"].astype(str)
    ext_rows["accession_or_name"] = external.get("example_protein_name_or_accession", "").astype(str)

    # -------------------------
    # Internal previous labels
    # Relabel using same 60C rule
    # -------------------------
    internal["topt_numeric"] = pd.to_numeric(internal["topt_numeric"], errors="coerce")

    internal["thermal_label_60c_aligned"] = np.where(
        internal["topt_numeric"] >= 60,
        "thermophilic",
        "mesophilic_lower"
    )

    # Merge internal with master to recover sequence.
    if "sequence" not in internal.columns:
        internal = internal.merge(
            master[["uniprot_accession", "sequence"]],
            on="uniprot_accession",
            how="left"
        )

    int_rows = pd.DataFrame()
    int_rows["sequence"] = internal["sequence"].apply(clean_sequence)
    int_rows["thermal_label"] = internal["thermal_label_60c_aligned"]
    int_rows["temperature_c"] = internal["topt_numeric"]
    int_rows["source_group"] = "internal_previous_relabelled_60c"
    int_rows["target_type"] = "internal_Topt_from_BRENDA"
    int_rows["source_detail"] = "previous_internal_labels_relabelled_from_topt_numeric_using_60C"
    int_rows["accession_or_name"] = internal["uniprot_accession"].astype(str)

    # Save internal relabel review.
    relabel_review_cols = [
        "uniprot_accession",
        "thermal_class",
        "topt_numeric",
        "brenda_temperature_optimum",
        "thermal_label_60c_aligned",
    ]

    internal[[c for c in relabel_review_cols if c in internal.columns]].to_csv(
        outdir / "internal_previous_labels_60c_relabel_review.csv",
        index=False
    )

    # -------------------------
    # Combine and clean
    # -------------------------
    all_rows = pd.concat([ext_rows, int_rows], ignore_index=True)

    all_rows["sequence"] = all_rows["sequence"].apply(clean_sequence)
    all_rows["sequence_length"] = all_rows["sequence"].str.len()
    all_rows["invalid_residue_count"] = all_rows["sequence"].apply(invalid_residue_count)
    all_rows["valid_standard_aa_only"] = all_rows["invalid_residue_count"] == 0
    all_rows["sequence_hash"] = all_rows["sequence"].apply(sequence_hash)

    all_rows = all_rows[
        (all_rows["sequence_length"] > 0)
        & (all_rows["valid_standard_aa_only"] == True)
        & (all_rows["thermal_label"].isin(["thermophilic", "mesophilic_lower"]))
    ].copy()

    all_rows.to_csv(outdir / "combined_60c_aligned_all_rows.csv", index=False)

    # -------------------------
    # Deduplicate by exact sequence
    # Remove conflicting labels
    # -------------------------
    dedup_records = []
    conflict_records = []

    for seq_hash, g in all_rows.groupby("sequence_hash"):
        labels = sorted(set(g["thermal_label"].astype(str)))

        base = {
            "sequence_hash": seq_hash,
            "sequence": g["sequence"].iloc[0],
            "sequence_length": int(g["sequence_length"].iloc[0]),
            "n_source_rows": len(g),
            "source_groups": safe_join(g["source_group"]),
            "target_types": safe_join(g["target_type"]),
            "source_details": safe_join(g["source_detail"], max_len=1000),
            "accession_or_name": safe_join(g["accession_or_name"], max_len=1000),
            "temperature_values_c": ";".join([str(x) for x in g["temperature_c"].dropna().astype(float).tolist()]),
            "labels_seen": ";".join(labels),
        }

        if len(labels) == 1:
            base["thermal_label"] = labels[0]
            dedup_records.append(base)
        else:
            conflict_records.append(base)

    dedup = pd.DataFrame(dedup_records)
    conflicts = pd.DataFrame(conflict_records)

    dedup.to_csv(outdir / "combined_60c_aligned_deduplicated_nonconflict.csv", index=False)
    conflicts.to_csv(outdir / "combined_60c_aligned_conflicting_sequence_labels.csv", index=False)

    # -------------------------
    # Train combined model
    # -------------------------
    feature_rows = [extract_features(seq) for seq in dedup["sequence"]]
    X = pd.DataFrame(feature_rows)
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))

    y = dedup["thermal_label"].map({
        "mesophilic_lower": 0,
        "thermophilic": 1,
    })

    X.to_csv(outdir / "combined_60c_aligned_sequence_features.csv", index=False)

    metrics, best_model_name = evaluate_models(
        X=X,
        y=y,
        outdir=outdir,
        target_name="Combined60CAligned"
    )

    # -------------------------
    # Report
    # -------------------------
    external_hashes = set(ext_rows["sequence"].apply(clean_sequence).apply(sequence_hash))
    internal_hashes = set(int_rows["sequence"].apply(clean_sequence).apply(sequence_hash))

    overlap_hashes = external_hashes & internal_hashes

    report = outdir / "COMBINED_60C_ALIGNED_ML_REPORT.md"

    with report.open("w") as fh:
        fh.write("# Combined 60C-aligned internal + external ML report\n\n")

        fh.write("## Label-definition correction\n\n")
        fh.write(
            "The previous internal thermal_class labels were not merged directly. "
            "They were recalculated from topt_numeric using the same 60C threshold used for the external Tmxyl/Toptxyl labels.\n\n"
        )

        fh.write("Internal relabel rule:\n\n")
        fh.write("- `topt_numeric >= 60` → thermophilic\n")
        fh.write("- `topt_numeric < 60` → mesophilic_lower\n\n")

        fh.write("## Input counts\n\n")
        fh.write(f"- External usable rows: {len(ext_rows)}\n")
        fh.write(f"- Internal previous rows after 60C relabel: {len(int_rows)}\n")
        fh.write(f"- Combined rows before deduplication: {len(all_rows)}\n")
        fh.write(f"- External/internal exact sequence overlaps: {len(overlap_hashes)}\n")
        fh.write(f"- Deduplicated non-conflicting training rows: {len(dedup)}\n")
        fh.write(f"- Conflicting sequence-label rows removed: {len(conflicts)}\n\n")

        fh.write("## Internal 60C relabel counts\n\n")
        fh.write(internal["thermal_label_60c_aligned"].value_counts().to_string())
        fh.write("\n\n")

        fh.write("## Final deduplicated training class counts\n\n")
        fh.write(dedup["thermal_label"].value_counts().to_string())
        fh.write("\n\n")

        fh.write("## Source composition after deduplication\n\n")
        fh.write(dedup["source_groups"].value_counts().head(20).to_string())
        fh.write("\n\n")

        fh.write("## Cross-validation metrics\n\n")
        display_cols = [
            "model",
            "n_samples",
            "class_0_mesophilic_lower",
            "class_1_thermophilic",
            "accuracy_mean",
            "balanced_accuracy_mean",
            "f1_macro_mean",
            "precision_macro_mean",
            "recall_macro_mean",
            "roc_auc_mean",
        ]

        fh.write(metrics[display_cols].to_string(index=False))
        fh.write("\n\n")

        fh.write(f"Best model by balanced accuracy: `{best_model_name}`\n\n")

        fh.write("## Interpretation\n\n")
        fh.write(
            "This combined classifier increases the number of labelled examples while keeping the thermal-class definition consistent. "
            "It should be interpreted as a broad 60C experimental/cross-curated thermal-class classifier. "
            "Because the dataset combines Tm, Topt, and BRENDA-derived Topt labels, the separate Tm-only and Topt-only models should still be reported alongside this combined model.\n"
        )

    print("[DONE] Combined 60C-aligned classifier trained.")
    print(f"Report: {report}")
    print(f"Deduplicated training set: {outdir / 'combined_60c_aligned_deduplicated_nonconflict.csv'}")
    print(f"Conflicts: {outdir / 'combined_60c_aligned_conflicting_sequence_labels.csv'}")


if __name__ == "__main__":
    main()
