#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
import pandas as pd


def minmax_score_lower_is_better(values):
    values = pd.to_numeric(values, errors="coerce")
    minimum = values.min()
    maximum = values.max()

    if pd.isna(minimum) or pd.isna(maximum):
        return pd.Series(np.nan, index=values.index)

    if maximum == minimum:
        return pd.Series(1.0, index=values.index)

    return (maximum - values) / (maximum - minimum)


def classify_foldx_ddg(ddg):
    if pd.isna(ddg):
        return "unknown"

    if ddg <= -0.5:
        return "stabilising"

    if ddg >= 0.5:
        return "destabilising"

    return "neutral"


def main():
    parser = argparse.ArgumentParser(
        description="Integrate FoldX, docking, and ML branch outputs."
    )
    parser.add_argument("--foldx-mutations", required=True)
    parser.add_argument("--docking", required=True)
    parser.add_argument("--ml-summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mutations = pd.read_csv(args.foldx_mutations)
    docking = pd.read_csv(args.docking)
    ml_summary = pd.read_csv(args.ml_summary)

    docking_ready = docking[docking["status"] == "ready"].copy()
    docking_ready["vina_affinity"] = pd.to_numeric(
        docking_ready["vina_affinity"],
        errors="coerce",
    )

    docking_summary = (
        docking_ready
        .groupby(
            ["uniprot_accession", "pdb_id", "foldx_mutation_code"],
            as_index=False,
        )
        .agg(
            mean_vina_affinity=("vina_affinity", "mean"),
            best_vina_affinity=("vina_affinity", "min"),
            docking_jobs=("vina_affinity", "count"),
            ligand_count=("ligand", "nunique"),
        )
    )

    integrated = mutations.merge(
        docking_summary,
        on=["uniprot_accession", "pdb_id", "foldx_mutation_code"],
        how="left",
    )

    integrated["foldx_ddg"] = pd.to_numeric(
        integrated["foldx_ddg"],
        errors="coerce",
    )

    integrated["foldx_class"] = integrated["foldx_ddg"].apply(classify_foldx_ddg)

    integrated["foldx_score"] = minmax_score_lower_is_better(
        integrated["foldx_ddg"]
    )

    integrated["docking_score"] = minmax_score_lower_is_better(
        integrated["mean_vina_affinity"]
    )

    integrated["best_docking_score"] = minmax_score_lower_is_better(
        integrated["best_vina_affinity"]
    )

    integrated["docking_complete"] = integrated["ligand_count"].fillna(0) >= 2
    integrated["docking_complete_score"] = integrated["docking_complete"].astype(float)

    integrated["integrated_score"] = (
        0.50 * integrated["foldx_score"].fillna(0)
        + 0.30 * integrated["docking_score"].fillna(0)
        + 0.10 * integrated["best_docking_score"].fillna(0)
        + 0.10 * integrated["docking_complete_score"].fillna(0)
    )

    integrated["candidate_category"] = np.where(
        (integrated["foldx_class"] == "stabilising")
        & (integrated["docking_complete"]),
        "priority_stabilising_candidate",
        np.where(
            (integrated["foldx_class"] == "neutral")
            & (integrated["docking_complete"]),
            "neutral_with_docking_support",
            "lower_priority_or_destabilising",
        ),
    )

    ready_ml_tasks = ml_summary[ml_summary["status"] == "ready"]["task"].nunique()
    skipped_ml_tasks = ml_summary[ml_summary["status"].str.startswith("skipped", na=False)]["task"].nunique()

    integrated["ml_ready_tasks"] = ready_ml_tasks
    integrated["ml_skipped_tasks"] = skipped_ml_tasks
    integrated["workflow_mode"] = "test_mode"

    sort_columns = [
        "integrated_score",
        "foldx_score",
        "docking_score",
        "best_docking_score",
    ]

    integrated = integrated.sort_values(sort_columns, ascending=False).reset_index(drop=True)
    integrated.insert(0, "rank", range(1, len(integrated) + 1))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    integrated.to_csv(output_path, index=False)

    print(f"Saved integration table: {output_path}")
    print(f"Rows: {len(integrated)}")
    print(integrated.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
