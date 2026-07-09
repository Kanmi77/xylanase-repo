#!/usr/bin/env python3
# Purpose: Train sequence-based ML model.

from pathlib import Path
import json
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.metrics import confusion_matrix, classification_report


AA = list("ACDEFGHIKLMNPQRSTVWY")

GROUPS = {
    "hydrophobic_fraction": set("AVILMFWY"),
    "aromatic_fraction": set("FWY"),
    "polar_fraction": set("STNQC"),
    "charged_fraction": set("KRHDE"),
    "acidic_fraction": set("DE"),
    "basic_fraction": set("KRH"),
    "tiny_fraction": set("AGSC"),
    "small_fraction": set("AGSTVCPDN"),
    "proline_fraction": set("P"),
    "glycine_fraction": set("G"),
    "cysteine_fraction": set("C"),
}


def extract_features(seq: str) -> dict:
    seq = str(seq).strip().upper()
    n = len(seq)

    feats = {
        "sequence_length": n,
    }

    if n == 0:
        for aa in AA:
            feats[f"aa_frac_{aa}"] = np.nan
        for group in GROUPS:
            feats[group] = np.nan
        return feats

    for aa in AA:
        feats[f"aa_frac_{aa}"] = seq.count(aa) / n

    for group_name, residues in GROUPS.items():
        feats[group_name] = sum(seq.count(r) for r in residues) / n

    # Simple ratios often relevant to stability.
    feats["charged_to_polar_ratio"] = (
        feats["charged_fraction"] / feats["polar_fraction"]
        if feats["polar_fraction"] > 0 else np.nan
    )

    feats["hydrophobic_to_polar_ratio"] = (
        feats["hydrophobic_fraction"] / feats["polar_fraction"]
        if feats["polar_fraction"] > 0 else np.nan
    )

    feats["acidic_to_basic_ratio"] = (
        feats["acidic_fraction"] / feats["basic_fraction"]
        if feats["basic_fraction"] > 0 else np.nan
    )

    return feats


def prepare_dataset(path: Path, target_name: str) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    df = pd.read_csv(path, low_memory=False)

    # Keep clean rows.
    df = df[
        (df["valid_standard_aa_only"] == True)
        & (df["thermal_label_60c_median"].isin(["thermophilic", "mesophilic_lower"]))
        & (df["has_conflicting_60c_labels"] == False)
    ].copy()

    df = df[df["sequence"].notna()].copy()
    df = df[df["sequence"].astype(str).str.len() > 0].copy()

    feature_rows = [extract_features(seq) for seq in df["sequence"]]
    X = pd.DataFrame(feature_rows)

    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))

    y = df["thermal_label_60c_median"].map({
        "mesophilic_lower": 0,
        "thermophilic": 1,
    })

    df["binary_label"] = y
    df["target_dataset"] = target_name

    return X, y, df


def evaluate_models(X: pd.DataFrame, y: pd.Series, outdir: Path, target_name: str):
    min_class_count = y.value_counts().min()
    n_splits = min(5, int(min_class_count))

    if n_splits < 2:
        raise RuntimeError(f"Not enough samples for CV in {target_name}")

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    models = {
        "dummy_majority": DummyClassifier(strategy="most_frequent"),
        "logistic_regression": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
                random_state=42,
            )),
        ]),
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
            min_samples_leaf=2,
        ),
    }

    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "f1_macro": "f1_macro",
        "precision_macro": "precision_macro",
        "recall_macro": "recall_macro",
        "roc_auc": "roc_auc",
    }

    rows = []

    for model_name, model in models.items():
        scores = cross_validate(
            model,
            X,
            y,
            cv=cv,
            scoring=scoring,
            return_train_score=False,
            n_jobs=-1,
            error_score=np.nan,
        )

        row = {
            "target_dataset": target_name,
            "model": model_name,
            "n_samples": len(y),
            "n_features": X.shape[1],
            "n_splits": n_splits,
            "class_0_mesophilic_lower": int((y == 0).sum()),
            "class_1_thermophilic": int((y == 1).sum()),
        }

        for metric in scoring:
            vals = scores[f"test_{metric}"]
            row[f"{metric}_mean"] = np.nanmean(vals)
            row[f"{metric}_sd"] = np.nanstd(vals)

        rows.append(row)

    metrics = pd.DataFrame(rows)
    metrics.to_csv(outdir / f"{target_name.lower()}_experimental_classifier_cv_metrics.csv", index=False)

    # Pick best model by balanced accuracy.
    best_row = metrics.sort_values("balanced_accuracy_mean", ascending=False).iloc[0]
    best_model_name = best_row["model"]
    best_model = models[best_model_name]

    best_model.fit(X, y)

    # Apparent training confusion matrix only for diagnostics, not final performance.
    y_pred = best_model.predict(X)

    cm = confusion_matrix(y, y_pred)
    report_dict = classification_report(
        y,
        y_pred,
        target_names=["mesophilic_lower", "thermophilic"],
        output_dict=True,
        zero_division=0,
    )

    with open(outdir / f"{target_name.lower()}_best_model_training_classification_report.json", "w") as fh:
        json.dump(report_dict, fh, indent=2)

    pd.DataFrame(
        cm,
        index=["true_mesophilic_lower", "true_thermophilic"],
        columns=["pred_mesophilic_lower", "pred_thermophilic"],
    ).to_csv(outdir / f"{target_name.lower()}_best_model_training_confusion_matrix.csv")

    # Feature importance.
    importance_rows = []

    if best_model_name == "random_forest":
        importances = best_model.feature_importances_
        for f, imp in zip(X.columns, importances):
            importance_rows.append({
                "feature": f,
                "importance": imp,
                "model": best_model_name,
            })

    elif best_model_name == "logistic_regression":
        clf = best_model.named_steps["clf"]
        coefs = clf.coef_[0]
        for f, coef in zip(X.columns, coefs):
            importance_rows.append({
                "feature": f,
                "importance": abs(coef),
                "coefficient": coef,
                "model": best_model_name,
            })

    if importance_rows:
        imp_df = pd.DataFrame(importance_rows).sort_values("importance", ascending=False)
        imp_df.to_csv(outdir / f"{target_name.lower()}_best_model_feature_importance.csv", index=False)

    return metrics, best_model_name


def main():
    root = Path(".").resolve()

    input_files = {
        "Tm": root / "results/experimental_ml_correction/tmxyl_sequence_deduplicated.csv",
        "Topt": root / "results/experimental_ml_correction/toptxyl_sequence_deduplicated.csv",
    }

    outdir = root / "results/experimental_ml_correction/sequence_classifier"
    outdir.mkdir(parents=True, exist_ok=True)

    all_metrics = []
    dataset_summaries = []

    for target_name, path in input_files.items():
        if not path.exists():
            raise FileNotFoundError(path)

        X, y, clean_df = prepare_dataset(path, target_name)

        clean_df.to_csv(outdir / f"{target_name.lower()}_classifier_input_records.csv", index=False)
        X.to_csv(outdir / f"{target_name.lower()}_sequence_features.csv", index=False)

        metrics, best_model_name = evaluate_models(X, y, outdir, target_name)

        all_metrics.append(metrics)

        dataset_summaries.append({
            "target_dataset": target_name,
            "input_file": str(path),
            "n_records_after_filtering": len(clean_df),
            "n_features": X.shape[1],
            "mesophilic_lower": int((y == 0).sum()),
            "thermophilic": int((y == 1).sum()),
            "best_model_by_balanced_accuracy": best_model_name,
        })

    all_metrics_df = pd.concat(all_metrics, ignore_index=True)
    all_metrics_df.to_csv(outdir / "experimental_sequence_classifier_all_cv_metrics.csv", index=False)

    summary_df = pd.DataFrame(dataset_summaries)
    summary_df.to_csv(outdir / "experimental_sequence_classifier_dataset_summary.csv", index=False)

    report = outdir / "EXPERIMENTAL_SEQUENCE_CLASSIFIER_REPORT.md"

    with report.open("w") as fh:
        fh.write("# Experimental-label sequence classifier report\n\n")
        fh.write("This ML correction uses experimental Tm and Topt labels, not FoldX-derived labels.\n\n")

        fh.write("## Dataset summary\n\n")
        fh.write(summary_df.to_string(index=False))
        fh.write("\n\n")

        fh.write("## Cross-validation metrics\n\n")
        display_cols = [
            "target_dataset",
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

        fh.write(all_metrics_df[display_cols].to_string(index=False))
        fh.write("\n\n")

        fh.write("## Interpretation\n\n")
        fh.write(
            "This resolves the circularity problem because the target variable is derived from experimental "
            "Tm/Topt measurements rather than from FoldX energy. The model should still be interpreted cautiously "
            "because the external dataset contains heterogeneous literature-derived records and may include engineered "
            "mutant variants as well as natural enzymes. Therefore, this classifier is best reported as a non-circular "
            "experimental-label ML correction layer, not as definitive proof of thermostability.\n"
        )

    print("[DONE] Experimental sequence classifiers trained.")
    print(f"Report: {report}")
    print(f"All metrics: {outdir / 'experimental_sequence_classifier_all_cv_metrics.csv'}")


if __name__ == "__main__":
    main()
