#!/usr/bin/env python3
# Purpose: Train sequence-structure ML model.

from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

OUTDIR = Path("results/experimental_ml_correction/sequence_structure_classifier")
OUTDIR.mkdir(parents=True, exist_ok=True)

META_PATH = Path("results/experimental_ml_correction/combined_60c_aligned_classifier/combined_60c_aligned_deduplicated_nonconflict.csv")
SEQ_FEATURE_PATH = Path("results/experimental_ml_correction/combined_60c_aligned_classifier/combined_60c_aligned_sequence_features.csv")
STRUCT_PATH = Path("results/structure_features/combined_structural_features_frozen_thesis.csv")
MAP_PATH = Path("results/experimental_ml_correction/combined_60c_sequence_structure_classifier/sequence_structure_mapping.csv")

meta = pd.read_csv(META_PATH)
seq_features = pd.read_csv(SEQ_FEATURE_PATH)
struct = pd.read_csv(STRUCT_PATH)
mapping = pd.read_csv(MAP_PATH)

if len(meta) != len(seq_features):
    raise ValueError(f"Metadata rows ({len(meta)}) and sequence-feature rows ({len(seq_features)}) do not match.")

seq = pd.concat([meta.reset_index(drop=True), seq_features.reset_index(drop=True)], axis=1)

seq = seq.merge(
    mapping[["sequence_hash", "chosen_structural_accession", "match_method", "has_final_structural_match"]],
    on="sequence_hash",
    how="left"
)

seq = seq[seq["has_final_structural_match"] == True].copy()

struct["uniprot_clean"] = (
    struct["uniprot_accession"]
    .astype(str)
    .str.replace("\xa0", "", regex=False)
    .str.strip()
)

struct_feature_cols = [
    "standard_residue_count",
    "atom_count",
    "hetatm_count",
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

struct_feature_cols = [c for c in struct_feature_cols if c in struct.columns]

struct_agg = (
    struct.groupby("uniprot_clean")
    .agg(
        **{f"mean_{c}": (c, "mean") for c in struct_feature_cols},
        n_structural_records=("structure_id", "count"),
        n_structure_sources=("structure_source_normalized", "nunique"),
    )
    .reset_index()
)

source_summary = (
    struct.assign(
        has_pdb=(struct["structure_source_normalized"].astype(str).str.lower() == "pdb").astype(int)
    )
    .groupby("uniprot_clean")
    .agg(has_pdb_structure=("has_pdb", "max"))
    .reset_index()
)

struct_agg = struct_agg.merge(source_summary, on="uniprot_clean", how="left")

merged = seq.merge(
    struct_agg,
    left_on="chosen_structural_accession",
    right_on="uniprot_clean",
    how="inner"
)

merged = merged[merged["thermal_label"].isin(["mesophilic_lower", "thermophilic"])].copy()
merged["target"] = (merged["thermal_label"] == "thermophilic").astype(int)

non_feature_cols = {
    "target",
    "thermal_label",
    "thermal_label_60c_aligned",
    "labels_seen",
    "sequence",
    "protein_sequence",
    "aa_sequence",
    "source",
    "source_groups",
    "source_details",
    "target_types",
    "accession_or_name",
    "temperature_values_c",
    "sequence_hash",
    "chosen_structural_accession",
    "has_final_structural_match",
    "match_method",
    "uniprot_clean",
}

numeric_cols = [
    c for c in merged.columns
    if c not in non_feature_cols and pd.api.types.is_numeric_dtype(merged[c])
]

X = merged[numeric_cols]
y = merged["target"]

class_counts = y.value_counts().to_dict()
n_splits = min(5, min(class_counts.values()))

models = {
    "logistic_regression": LogisticRegression(max_iter=5000, class_weight="balanced"),
    "random_forest": RandomForestClassifier(
        n_estimators=500,
        random_state=42,
        class_weight="balanced",
        n_jobs=-1,
    ),
}

preprocess = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
])

scoring = {
    "accuracy": "accuracy",
    "balanced_accuracy": "balanced_accuracy",
    "f1_macro": "f1_macro",
    "precision_macro": "precision_macro",
    "recall_macro": "recall_macro",
    "roc_auc": "roc_auc",
}

cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

rows = []

for name, model in models.items():
    pipe = Pipeline([
        ("preprocess", preprocess),
        ("model", model),
    ])

    cvres = cross_validate(pipe, X, y, cv=cv, scoring=scoring, n_jobs=-1)

    row = {
        "model": name,
        "n_samples": len(merged),
        "class_0_mesophilic_lower": int((y == 0).sum()),
        "class_1_thermophilic": int((y == 1).sum()),
        "n_features": len(numeric_cols),
    }

    for metric in scoring:
        vals = cvres[f"test_{metric}"]
        row[f"{metric}_mean"] = float(np.mean(vals))
        row[f"{metric}_sd"] = float(np.std(vals))

    rows.append(row)

metrics = pd.DataFrame(rows)
metrics.to_csv(OUTDIR / "sequence_structure_classifier_cv_metrics.csv", index=False)

best_name = metrics.sort_values("balanced_accuracy_mean", ascending=False).iloc[0]["model"]
best_model = models[best_name]

best_pipe = Pipeline([
    ("preprocess", preprocess),
    ("model", best_model),
])
best_pipe.fit(X, y)

if best_name == "random_forest":
    importance_values = best_pipe.named_steps["model"].feature_importances_
else:
    importance_values = np.abs(best_pipe.named_steps["model"].coef_[0])

importance = (
    pd.DataFrame({
        "feature": numeric_cols,
        "importance": importance_values,
        "model": best_name,
    })
    .sort_values("importance", ascending=False)
)

importance.to_csv(OUTDIR / "sequence_structure_feature_importance.csv", index=False)
merged.to_csv(OUTDIR / "sequence_structure_training_rows.csv", index=False)

report = OUTDIR / "SEQUENCE_STRUCTURE_60C_ML_REPORT.md"

with open(report, "w") as h:
    h.write("# Final maximized sequence + structural 60C classifier report\n\n")
    h.write(f"- Merged labelled sequence+structure rows: {len(merged)}\n")
    h.write(f"- Numeric features used: {len(numeric_cols)}\n")
    h.write(f"- CV folds: {n_splits}\n")
    h.write(f"- Mesophilic_lower: {(y == 0).sum()}\n")
    h.write(f"- Thermophilic: {(y == 1).sum()}\n\n")

    h.write("## Match method counts\n\n")
    h.write(merged["match_method"].value_counts().to_markdown())
    h.write("\n\n")

    h.write("## Cross-validation metrics\n\n")
    h.write(metrics.to_markdown(index=False))
    h.write("\n\n")

    h.write(f"Best model by balanced accuracy: `{best_name}`\n\n")

    h.write("## Top feature importance values\n\n")
    h.write(importance.head(30).to_markdown(index=False))
    h.write("\n")

print("Done")
print("Report:", report)
print("Rows:", len(merged))
print("Class counts:")
print(merged["thermal_label"].value_counts())
print("Best model:", best_name)
