#!/usr/bin/env python3
"""
Train wild-type experimental stability ML models.

Target:
    Experimental Tm/Topt median temperature in Celsius.

FoldX is not used as target.
This script focuses on the clean WT experimental sequence panel.

Feature sets:
    1. sequence_only
    2. sequence_plus_gh
    3. sequence_plus_metadata
    4. sequence_plus_metadata_assay
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


warnings.filterwarnings("ignore")


BASE_DIR = Path.home() / "xylanase-thesis"

INPUT_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "datasets"
    / "wt_experimental_stability_panel_enriched.csv"
)

OUTPUT_DIR = BASE_DIR / "results" / "stability_ml_goal" / "models"

REPORT_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "reports"
    / "wt_experimental_stability_ml_report.md"
)

METRICS_FILE = OUTPUT_DIR / "wt_experimental_stability_cv_metrics.csv"
SUMMARY_FILE = OUTPUT_DIR / "wt_experimental_stability_cv_summary.csv"
PREDICTIONS_FILE = OUTPUT_DIR / "wt_experimental_stability_predictions.csv"
IMPORTANCE_FILE = OUTPUT_DIR / "wt_experimental_stability_feature_importance.csv"


try:
    from xgboost import XGBRegressor

    HAS_XGBOOST = True
except ImportError:
    HAS_XGBOOST = False


def make_one_hot_encoder():
    """Create one-hot encoder compatible with older/newer sklearn."""
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def get_sequence_columns(data):
    """Select sequence descriptor columns."""
    columns = []

    direct_columns = [
        "seq_length_calculated",
        "molecular_weight",
        "isoelectric_point",
        "aromaticity",
        "instability_index",
        "gravy",
        "helix_fraction_protparam",
        "turn_fraction_protparam",
        "sheet_fraction_protparam",
        "hydrophobic_fraction",
        "aromatic_fraction",
        "polar_fraction",
        "charged_fraction",
        "acidic_fraction",
        "basic_fraction",
        "tiny_fraction",
        "small_fraction",
        "proline_fraction",
        "glycine_fraction",
        "cysteine_fraction",
        "charged_to_polar_ratio",
        "hydrophobic_to_polar_ratio",
        "acidic_to_basic_ratio",
    ]

    columns.extend([col for col in direct_columns if col in data.columns])
    columns.extend([col for col in data.columns if col.startswith("aa_frac_")])

    return list(dict.fromkeys(columns))


def make_preprocessor(data, feature_columns):
    """Create preprocessing pipeline."""
    selected = data[feature_columns].copy()

    numeric_columns = selected.select_dtypes(
        include=["number", "bool"]
    ).columns.tolist()

    categorical_columns = [
        column for column in feature_columns
        if column not in numeric_columns
    ]

    transformers = []

    if numeric_columns:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                    ]
                ),
                numeric_columns,
            )
        )

    if categorical_columns:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", make_one_hot_encoder()),
                    ]
                ),
                categorical_columns,
            )
        )

    return ColumnTransformer(transformers=transformers)


def get_models():
    """Define models."""
    models = {
        "dummy_mean": DummyRegressor(strategy="mean"),
        "ridge": Ridge(alpha=1.0),
        "extra_trees": ExtraTreesRegressor(
            n_estimators=600,
            random_state=42,
            min_samples_leaf=2,
            max_features="sqrt",
            n_jobs=-1,
        ),
        "hist_gradient_boosting": HistGradientBoostingRegressor(
            random_state=42,
            max_iter=300,
            learning_rate=0.04,
            l2_regularization=0.1,
        ),
    }

    if HAS_XGBOOST:
        models["xgboost"] = XGBRegressor(
            n_estimators=500,
            learning_rate=0.04,
            max_depth=3,
            subsample=0.85,
            colsample_bytree=0.85,
            objective="reg:squarederror",
            random_state=42,
            n_jobs=-1,
        )

    return models


def rmse(y_true, y_pred):
    """Calculate RMSE."""
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


def calculate_metrics(y_true, y_pred):
    """Calculate regression metrics."""
    y_true_array = np.asarray(y_true, dtype=float)
    y_pred_array = np.asarray(y_pred, dtype=float)

    return {
        "mae": mean_absolute_error(y_true_array, y_pred_array),
        "rmse": rmse(y_true_array, y_pred_array),
        "r2": r2_score(y_true_array, y_pred_array),
        "spearman": pd.Series(y_true_array).corr(
            pd.Series(y_pred_array),
            method="spearman",
        ),
    }


def run_cv(data, feature_set_name, feature_columns):
    """Run grouped 5-fold CV."""
    y = data["target_temperature_c"].astype(float)
    groups = data["sequence_hash"].astype(str)

    splitter = GroupKFold(n_splits=5)

    metric_rows = []
    prediction_rows = []

    for model_name, model in get_models().items():
        for fold, (train_idx, test_idx) in enumerate(
            splitter.split(data, y, groups),
            start=1,
        ):
            train = data.iloc[train_idx].copy()
            test = data.iloc[test_idx].copy()

            x_train = train[feature_columns]
            x_test = test[feature_columns]

            y_train = y.iloc[train_idx]
            y_test = y.iloc[test_idx]

            pipeline = Pipeline(
                steps=[
                    ("preprocess", make_preprocessor(data, feature_columns)),
                    ("model", model),
                ]
            )

            pipeline.fit(x_train, y_train)
            pred = pipeline.predict(x_test)

            metrics = calculate_metrics(y_test, pred)

            row = {
                "feature_set": feature_set_name,
                "model": model_name,
                "fold": fold,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
            }
            row.update(metrics)
            metric_rows.append(row)

            for row_index, true_value, predicted_value in zip(
                test.index,
                y_test,
                pred,
            ):
                prediction_rows.append(
                    {
                        "row_index": row_index,
                        "feature_set": feature_set_name,
                        "model": model_name,
                        "fold": fold,
                        "sequence_hash": test.loc[row_index, "sequence_hash"],
                        "target_type": test.loc[row_index, "target_type"],
                        "thermal_label": test.loc[row_index, "thermal_label"],
                        "true_temperature_c": true_value,
                        "predicted_temperature_c": predicted_value,
                        "residual_c": true_value - predicted_value,
                    }
                )

    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def summarise_metrics(metrics):
    """Summarise metrics across folds."""
    summary = (
        metrics.groupby(["feature_set", "model"])[
            ["mae", "rmse", "r2", "spearman"]
        ]
        .agg(["mean", "std"])
        .reset_index()
    )

    summary.columns = [
        "_".join(col).strip("_")
        for col in summary.columns.to_flat_index()
    ]

    return summary


def get_feature_names(pipeline):
    """Get transformed feature names."""
    preprocessor = pipeline.named_steps["preprocess"]

    try:
        return preprocessor.get_feature_names_out()
    except Exception:
        names = []

        for _, _, columns in preprocessor.transformers_:
            if isinstance(columns, list):
                names.extend(columns)

        return names


def train_importance_model(data, feature_columns):
    """Train ExtraTrees on the richest feature set and save importance."""
    y = data["target_temperature_c"].astype(float)

    model = ExtraTreesRegressor(
        n_estimators=600,
        random_state=42,
        min_samples_leaf=2,
        max_features="sqrt",
        n_jobs=-1,
    )

    pipeline = Pipeline(
        steps=[
            ("preprocess", make_preprocessor(data, feature_columns)),
            ("model", model),
        ]
    )

    pipeline.fit(data[feature_columns], y)

    feature_names = get_feature_names(pipeline)
    importances = pipeline.named_steps["model"].feature_importances_

    importance = pd.DataFrame(
        {
            "feature": feature_names,
            "importance": importances,
        }
    ).sort_values("importance", ascending=False)

    importance.to_csv(IMPORTANCE_FILE, index=False)

    return importance


def write_report(data, feature_sets, summary, importance):
    """Write markdown report."""
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    best_rows = (
        summary.sort_values("rmse_mean", ascending=True)
        .groupby("feature_set")
        .head(1)
        .copy()
    )

    lines = [
        "# Wild-Type Experimental Stability ML Report",
        "",
        "## Dataset",
        "",
        f"Rows used: {len(data)}",
        f"Unique sequence hashes: {data['sequence_hash'].nunique()}",
        "",
        "## Target",
        "",
        "Target: experimental Tm/Topt median temperature in Celsius.",
        "",
        data["target_temperature_c"].describe().to_string(),
        "",
        "## Feature sets",
        "",
    ]

    for name, columns in feature_sets.items():
        lines.append(f"- {name}: {len(columns)} columns")

    lines.extend(
        [
            "",
            "## Best model per feature set",
            "",
            best_rows[
                [
                    "feature_set",
                    "model",
                    "mae_mean",
                    "rmse_mean",
                    "r2_mean",
                    "spearman_mean",
                ]
            ].to_string(index=False),
            "",
            "## Full CV summary",
            "",
            summary.to_string(index=False),
            "",
            "## Top ExtraTrees features",
            "",
            importance.head(30).to_string(index=False),
            "",
            "## Interpretation note",
            "",
            "Organism type was included only in sensitivity feature sets because many rows remained unknown after metadata enrichment.",
            "FoldX WT stability was not used as target. It should be evaluated separately as a computational predictor on the matched subset.",
            "",
        ]
    )

    REPORT_FILE.write_text("\n".join(lines) + "\n")


def main():
    """Run WT stability ML."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(INPUT_FILE)

    data = data[
        data["target_temperature_c"].notna()
        & (data["target_temperature_c"] >= 20)
        & (data["target_temperature_c"] <= 120)
        & data["thermal_label"].isin(["thermophilic", "mesophilic_lower"])
    ].copy()

    data["gh_family_enriched"] = data["gh_family_enriched"].fillna("unknown")
    data["organism_type_enriched"] = data[
        "organism_type_enriched"
    ].fillna("unknown")
    data["metadata_quality"] = data["metadata_quality"].fillna("unknown")

    sequence_columns = get_sequence_columns(data)

    feature_sets = {
        "sequence_only": sequence_columns,
        "sequence_plus_gh": sequence_columns + ["gh_family_enriched"],
        "sequence_plus_metadata": sequence_columns
        + [
            "gh_family_enriched",
            "organism_type_enriched",
            "metadata_quality",
        ],
        "sequence_plus_metadata_assay": sequence_columns
        + [
            "gh_family_enriched",
            "organism_type_enriched",
            "metadata_quality",
            "target_type",
        ],
    }

    print(f"Rows used: {len(data)}")
    print(f"Sequence features: {len(sequence_columns)}")
    print(f"XGBoost available: {HAS_XGBOOST}")

    all_metrics = []
    all_predictions = []

    for feature_set_name, feature_columns in feature_sets.items():
        print(f"\nRunning: {feature_set_name}")
        print(f"Columns: {len(feature_columns)}")

        metrics, predictions = run_cv(
            data,
            feature_set_name,
            feature_columns,
        )

        all_metrics.append(metrics)
        all_predictions.append(predictions)

    metrics = pd.concat(all_metrics, ignore_index=True)
    predictions = pd.concat(all_predictions, ignore_index=True)

    summary = summarise_metrics(metrics)

    metrics.to_csv(METRICS_FILE, index=False)
    summary.to_csv(SUMMARY_FILE, index=False)
    predictions.to_csv(PREDICTIONS_FILE, index=False)

    importance = train_importance_model(
        data,
        feature_sets["sequence_plus_metadata_assay"],
    )

    write_report(data, feature_sets, summary, importance)

    print("\nWrote:")
    print(METRICS_FILE)
    print(SUMMARY_FILE)
    print(PREDICTIONS_FILE)
    print(IMPORTANCE_FILE)
    print(REPORT_FILE)

    print("\nCV summary:")
    print(
        summary[
            [
                "feature_set",
                "model",
                "mae_mean",
                "rmse_mean",
                "r2_mean",
                "spearman_mean",
            ]
        ]
        .sort_values("rmse_mean")
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
