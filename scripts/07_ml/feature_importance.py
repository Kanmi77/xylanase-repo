#!/usr/bin/env python3
# Purpose: Extract feature importance values.

from pathlib import Path
import pandas as pd
import numpy as np

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

OUTDIR = Path("results/experimental_ml_correction/direct217_model_improvement")
OUTDIR.mkdir(parents=True, exist_ok=True)

TRAIN_PATH = Path("results/experimental_ml_correction/sequence_structure_classifier/sequence_structure_training_rows.csv")

df = pd.read_csv(TRAIN_PATH)

df = df[df["thermal_label"].isin(["thermophilic", "mesophilic_lower"])].copy()
df["target"] = (df["thermal_label"] == "thermophilic").astype(int)

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
    "label_source_type",
    "label_confidence",
}

all_numeric = [
    c for c in df.columns
    if c not in non_feature_cols and pd.api.types.is_numeric_dtype(df[c])
]

sequence_cols = [
    c for c in all_numeric
    if (
        c.startswith("aa_frac_")
        or c.endswith("_fraction")
        or c.endswith("_ratio")
        or c == "sequence_length"
    )
]

normalised_struct_cols = [
    c for c in all_numeric
    if c in [
        "mean_hbond_proxy_per_residue",
        "mean_salt_bridge_per_residue",
        "mean_disulfide_per_residue",
        "mean_sasa_per_residue",
        "has_pdb_structure",
    ]
]

feature_cols = list(dict.fromkeys(sequence_cols + normalised_struct_cols))

X = df[feature_cols]
y = df["target"]

model = ExtraTreesClassifier(
    n_estimators=1000,
    max_features="sqrt",
    min_samples_leaf=2,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1,
)

pipe = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
    ("model", model),
])

pipe.fit(X, y)

importances = pipe.named_steps["model"].feature_importances_

importance = (
    pd.DataFrame({
        "feature": feature_cols,
        "importance": importances,
        "feature_group": [
            "normalised_structural" if c in normalised_struct_cols else "sequence"
            for c in feature_cols
        ]
    })
    .sort_values("importance", ascending=False)
)

importance.to_csv(
    OUTDIR / "selected_extratrees_sequence_plus_normalised_structural_feature_importance.csv",
    index=False
)

print("Top 30 selected-model feature importances:")
print(importance.head(30).to_string(index=False))

print("\nFeature-group importance totals:")
print(importance.groupby("feature_group")["importance"].sum().sort_values(ascending=False))

print("\nSaved:")
print(OUTDIR / "selected_extratrees_sequence_plus_normalised_structural_feature_importance.csv")
