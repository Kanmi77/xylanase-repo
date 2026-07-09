#!/usr/bin/env python3

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd
import yaml

from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.dummy import DummyRegressor, DummyClassifier
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def safe_numeric_frame(df):
    numeric = df.select_dtypes(include=["number"]).copy()

    drop_columns = [
        column for column in numeric.columns
        if numeric[column].isna().all()
    ]

    if drop_columns:
        numeric = numeric.drop(columns=drop_columns)

    return numeric


def regression_metrics(y_true, y_pred):
    return {
        "mae": mean_absolute_error(y_true, y_pred),
        "rmse": mean_squared_error(y_true, y_pred) ** 0.5,
        "r2": r2_score(y_true, y_pred) if len(y_true) > 1 else np.nan,
    }


def run_regression_task(task_name, dataset, target_column, output_dir):
    dataset = dataset.copy()
    dataset[target_column] = pd.to_numeric(dataset[target_column], errors="coerce")
    dataset = dataset.dropna(subset=[target_column])

    feature_data = safe_numeric_frame(dataset.drop(columns=[target_column], errors="ignore"))

    if target_column in feature_data.columns:
        feature_data = feature_data.drop(columns=[target_column])

    valid_rows = feature_data.notna().any(axis=1)
    feature_data = feature_data[valid_rows]
    target = dataset.loc[valid_rows, target_column]

    if len(feature_data) < 6 or feature_data.shape[1] == 0:
        return {
            "task": task_name,
            "target": target_column,
            "model": "skipped",
            "n_samples": len(feature_data),
            "n_features": feature_data.shape[1],
            "status": "skipped_insufficient_data",
        }, pd.DataFrame(), pd.DataFrame()

    test_size = 0.33 if len(feature_data) >= 9 else 0.40

    x_train, x_test, y_train, y_test = train_test_split(
        feature_data,
        target,
        test_size=test_size,
        random_state=42,
    )

    models = {
        "dummy_mean": DummyRegressor(strategy="mean"),
        "random_forest": RandomForestRegressor(
            n_estimators=200,
            random_state=42,
            min_samples_leaf=1,
        ),
    }

    rows = []
    predictions = []
    importance_rows = []

    for model_name, model in models.items():
        pipe = Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", model),
            ]
        )

        pipe.fit(x_train, y_train)
        y_pred = pipe.predict(x_test)

        metrics = regression_metrics(y_test, y_pred)

        rows.append(
            {
                "task": task_name,
                "target": target_column,
                "model": model_name,
                "n_samples": len(feature_data),
                "n_features": feature_data.shape[1],
                "test_samples": len(x_test),
                "mae": metrics["mae"],
                "rmse": metrics["rmse"],
                "r2": metrics["r2"],
                "status": "ready",
            }
        )

        pred_df = pd.DataFrame(
            {
                "task": task_name,
                "model": model_name,
                "observed": y_test.values,
                "predicted": y_pred,
            }
        )
        predictions.append(pred_df)

        fitted_model = pipe.named_steps["model"]

        if hasattr(fitted_model, "feature_importances_"):
            importance = pd.DataFrame(
                {
                    "task": task_name,
                    "model": model_name,
                    "feature": feature_data.columns,
                    "importance": fitted_model.feature_importances_,
                }
            )
            importance_rows.append(importance)

    summary = pd.DataFrame(rows)
    prediction_table = pd.concat(predictions, ignore_index=True)
    importance_table = (
        pd.concat(importance_rows, ignore_index=True)
        if importance_rows
        else pd.DataFrame()
    )

    return summary, prediction_table, importance_table


def run_classification_task(task_name, dataset, target_column, output_dir):
    dataset = dataset.copy()
    dataset[target_column] = pd.to_numeric(dataset[target_column], errors="coerce")
    dataset = dataset.dropna(subset=[target_column])

    lower = -0.5
    upper = 0.5

    dataset["foldx_class"] = np.where(
        dataset[target_column] <= lower,
        "stabilising",
        np.where(dataset[target_column] >= upper, "destabilising", "neutral"),
    )

    classified = dataset[dataset["foldx_class"] != "neutral"].copy()

    if classified["foldx_class"].nunique() < 2 or len(classified) < 8:
        return {
            "task": task_name,
            "target": "foldx_class",
            "model": "skipped",
            "n_samples": len(classified),
            "n_features": 0,
            "status": "skipped_insufficient_classes_or_samples",
        }, pd.DataFrame()

    feature_data = safe_numeric_frame(
        classified.drop(columns=[target_column, "foldx_class"], errors="ignore")
    )

    valid_rows = feature_data.notna().any(axis=1)
    feature_data = feature_data[valid_rows]
    target = classified.loc[valid_rows, "foldx_class"]

    if target.nunique() < 2:
        return {
            "task": task_name,
            "target": "foldx_class",
            "model": "skipped",
            "n_samples": len(feature_data),
            "n_features": feature_data.shape[1],
            "status": "skipped_single_class_after_filtering",
        }, pd.DataFrame()

    x_train, x_test, y_train, y_test = train_test_split(
        feature_data,
        target,
        test_size=0.33,
        random_state=42,
        stratify=target,
    )

    models = {
        "dummy_most_frequent": DummyClassifier(strategy="most_frequent"),
        "random_forest": RandomForestClassifier(
            n_estimators=200,
            random_state=42,
            class_weight="balanced",
        ),
    }

    rows = []

    for model_name, model in models.items():
        pipe = Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", model),
            ]
        )

        pipe.fit(x_train, y_train)
        y_pred = pipe.predict(x_test)

        rows.append(
            {
                "task": task_name,
                "target": "foldx_class",
                "model": model_name,
                "n_samples": len(feature_data),
                "n_features": feature_data.shape[1],
                "test_samples": len(x_test),
                "accuracy": accuracy_score(y_test, y_pred),
                "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
                "macro_f1": f1_score(y_test, y_pred, average="macro"),
                "status": "ready",
            }
        )

    return pd.DataFrame(rows), pd.DataFrame()


def write_report(output_dir, summary, dataset_shapes):
    report = output_dir / "ml_report.md"

    lines = [
        "# Machine-learning workflow report",
        "",
        "This report was generated by the Snakemake ML branch.",
        "",
        "## Input dataset sizes",
        "",
    ]

    for name, shape in dataset_shapes.items():
        lines.append(f"- {name}: {shape[0]} rows × {shape[1]} columns")

    lines.extend(
        [
            "",
            "## Model summary",
            "",
            summary.to_markdown(index=False),
            "",
            "## Interpretation note",
            "",
            "When the workflow is run in test mode, ML results only confirm that the branch executes correctly. "
            "Final thesis interpretation requires the full FoldX and docking datasets.",
            "",
        ]
    )

    report.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Run ML proof-of-workflow models for FoldX and docking outputs."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--foldx-wt", required=True)
    parser.add_argument("--foldx-mutations", required=True)
    parser.add_argument("--docking", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    read_yaml(args.config)

    output_path = Path(args.output)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    foldx_wt = pd.read_csv(args.foldx_wt)
    mutations = pd.read_csv(args.foldx_mutations)
    docking = pd.read_csv(args.docking)

    docking_ready = docking[docking["status"] == "ready"].copy()

    docking_summary = (
        docking_ready
        .groupby(["uniprot_accession", "pdb_id", "foldx_mutation_code"], as_index=False)
        .agg(
            mean_vina_affinity=("vina_affinity", "mean"),
            best_vina_affinity=("vina_affinity", "min"),
            docking_jobs=("vina_affinity", "count"),
        )
    )

    mutation_dataset = mutations.merge(
        docking_summary,
        on=["uniprot_accession", "pdb_id", "foldx_mutation_code"],
        how="left",
    )

    wt_summary, wt_predictions, wt_importance = run_regression_task(
        "wt_foldx_energy_regression",
        foldx_wt,
        "foldx_energy_per_residue",
        output_dir,
    )

    ddg_summary, ddg_predictions, ddg_importance = run_regression_task(
        "mutation_ddg_regression",
        mutation_dataset,
        "foldx_ddg",
        output_dir,
    )

    docking_reg_summary, docking_predictions, docking_importance = run_regression_task(
        "docking_affinity_regression",
        mutation_dataset.dropna(subset=["mean_vina_affinity"]),
        "mean_vina_affinity",
        output_dir,
    )

    class_summary, class_predictions = run_classification_task(
        "foldx_mutation_classification",
        mutation_dataset,
        "foldx_ddg",
        output_dir,
    )

    summary_parts = []

    for part in [wt_summary, ddg_summary, docking_reg_summary, class_summary]:
        if isinstance(part, dict):
            summary_parts.append(pd.DataFrame([part]))
        else:
            summary_parts.append(part)

    summary = pd.concat(summary_parts, ignore_index=True, sort=False)
    summary.to_csv(output_path, index=False)

    prediction_tables = [
        table for table in [wt_predictions, ddg_predictions, docking_predictions, class_predictions]
        if isinstance(table, pd.DataFrame) and not table.empty
    ]

    if prediction_tables:
        pd.concat(prediction_tables, ignore_index=True).to_csv(
            output_dir / "ml_predictions.csv",
            index=False,
        )

    importance_tables = [
        table for table in [wt_importance, ddg_importance, docking_importance]
        if isinstance(table, pd.DataFrame) and not table.empty
    ]

    if importance_tables:
        pd.concat(importance_tables, ignore_index=True).to_csv(
            output_dir / "ml_feature_importance.csv",
            index=False,
        )

    mutation_dataset.to_csv(output_dir / "ml_ready_mutation_dataset.csv", index=False)

    dataset_shapes = {
        "foldx_wt": foldx_wt.shape,
        "foldx_mutations": mutations.shape,
        "docking": docking.shape,
        "ml_ready_mutation_dataset": mutation_dataset.shape,
    }

    write_report(output_dir, summary, dataset_shapes)

    print(f"Saved ML summary: {output_path}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
