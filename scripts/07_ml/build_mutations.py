#!/usr/bin/env python3
# Purpose: Build mutation prioritisation panel.

from pathlib import Path
import re

import numpy as np
import pandas as pd


BASE_DIR = Path.home() / "xylanase-thesis"

HOMOLOGY_FOLDX_FILE = (
    BASE_DIR
    / "results"
    / "foldx_clean"
    / "tier2_ddg_ranked_industrial_annotated.csv"
)

PDB_FOLDX_FILE = (
    BASE_DIR
    / "results"
    / "optionC_original_only"
    / "pdb60_foldx_mutations"
    / "parsed_ddg_and_docking"
    / "pdb60_foldx_ddg_parsed_all.csv"
)

DOCKING_FILE = (
    BASE_DIR
    / "results"
    / "optionC_original_only"
    / "pdb60_foldx_mutations"
    / "parsed_ddg_and_docking"
    / "paired_mutant_wt_docking_delta"
    / "pdb60_paired_mutant_wt_docking_delta.csv"
)

INTEGRATED_FILE = (
    BASE_DIR
    / "results"
    / "integration"
    / "final_integrated_candidate_ranking.csv"
)

EVIDENCE_FILE = (
    BASE_DIR
    / "results"
    / "integration"
    / "stratified_pdb_vs_homology_new_ml_pdb_foldx_docking"
    / "final_evidence_tiered_summary_annotated_metadata_filled.csv"
)

OUTPUT_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "datasets"
    / "mutation_prioritisation_panel.csv"
)

TOP15_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "candidates"
    / "top15_mutation_candidates_goal.csv"
)

REPORT_FILE = (
    BASE_DIR
    / "results"
    / "stability_ml_goal"
    / "reports"
    / "mutation_prioritisation_panel_report.md"
)


def read_table(path):
    """Read a CSV/TSV table."""
    if not path.exists():
        return pd.DataFrame()

    try:
        return pd.read_csv(path)
    except Exception:
        return pd.read_csv(path, sep="\t")


def clean_text(value):
    """Clean general text."""
    if pd.isna(value):
        return ""
    return str(value).strip()


def clean_protein(value):
    """Clean protein/accession identifier."""
    text = clean_text(value).replace("\xa0", "")
    text = re.sub(r"\s+", "", text)
    return text.upper()


def clean_mutation(value):
    """Clean FoldX mutation code."""
    text = clean_text(value).replace(";", "")
    text = re.sub(r"\s+", "", text)
    return text.upper()


def make_candidate_id(protein, mutation):
    """Create stable candidate key."""
    protein = clean_protein(protein)
    mutation = clean_mutation(mutation)

    if not protein or not mutation:
        return ""

    return f"{protein}_{mutation}"


def first_existing(columns, candidates):
    """Return the first available column name."""
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def normalise_lower_is_better(values):
    """Normalise a numeric series where lower values are better."""
    numeric = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=numeric.index, dtype=float)

    valid = numeric.dropna()

    if valid.empty:
        return result

    minimum = valid.min()
    maximum = valid.max()

    if maximum == minimum:
        result.loc[valid.index] = 0.5
        return result

    result.loc[valid.index] = (maximum - valid) / (maximum - minimum)

    return result


def normalise_higher_is_better(values):
    """Normalise a numeric series where higher values are better."""
    numeric = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=numeric.index, dtype=float)

    valid = numeric.dropna()

    if valid.empty:
        return result

    minimum = valid.min()
    maximum = valid.max()

    if maximum == minimum:
        result.loc[valid.index] = 0.5
        return result

    result.loc[valid.index] = (valid - minimum) / (maximum - minimum)

    return result


def standardise_homology_foldx():
    """Load and standardise homology FoldX mutation results."""
    data = read_table(HOMOLOGY_FOLDX_FILE)

    if data.empty:
        return pd.DataFrame()

    out = pd.DataFrame()
    out["protein"] = data["protein"].apply(clean_protein)
    out["uniprot_accession"] = data["uniprot_accession"].apply(clean_protein)
    out["mutation"] = data["mutation"].apply(clean_mutation)
    out["candidate_id"] = [
        make_candidate_id(protein, mutation)
        for protein, mutation in zip(out["protein"], out["mutation"])
    ]

    out["foldx_ddg"] = pd.to_numeric(data["ddg"], errors="coerce")
    out["foldx_effect"] = data["effect"].astype(str)
    out["mutant_model"] = data["mutant_model"].astype(str)
    out["foldx_structure_branch"] = "homology_model"
    out["foldx_source_file"] = str(HOMOLOGY_FOLDX_FILE.relative_to(BASE_DIR))

    optional_columns = [
        "organism",
        "organism_type",
        "gh_family",
        "foldx_energy_per_residue",
        "sasa_per_res",
        "hbond_per_res",
        "disulfide_per_res",
        "brenda_temperature_optimum",
        "brenda_temperature_range",
        "brenda_temperature_stability",
        "brenda_ph_optimum",
        "brenda_ph_range",
    ]

    for column in optional_columns:
        if column in data.columns:
            out[column] = data[column]

    return out[out["candidate_id"].ne("")].copy()


def standardise_pdb_foldx():
    """Load and standardise PDB FoldX mutation results, if available."""
    data = read_table(PDB_FOLDX_FILE)

    if data.empty:
        return pd.DataFrame()

    protein_col = first_existing(
        data.columns,
        [
            "uniprot_accession",
            "protein",
            "candidate_accession",
            "std_uniprot",
        ],
    )

    mutation_col = first_existing(
        data.columns,
        [
            "mutation",
            "foldx_mutation_code",
            "std_mutation_x",
            "std_mutation_y",
        ],
    )

    ddg_col = first_existing(
        data.columns,
        [
            "foldx_ddg",
            "ddg",
            "delta_delta_g",
            "total_ddg",
        ],
    )

    if protein_col is None or mutation_col is None or ddg_col is None:
        return pd.DataFrame()

    out = pd.DataFrame()
    out["protein"] = data[protein_col].apply(clean_protein)
    out["uniprot_accession"] = data[protein_col].apply(clean_protein)
    out["mutation"] = data[mutation_col].apply(clean_mutation)
    out["candidate_id"] = [
        make_candidate_id(protein, mutation)
        for protein, mutation in zip(out["protein"], out["mutation"])
    ]

    out["foldx_ddg"] = pd.to_numeric(data[ddg_col], errors="coerce")
    out["foldx_effect"] = np.where(
        out["foldx_ddg"] < 0,
        "stabilizing",
        "destabilizing",
    )

    out["mutant_model"] = (
        data["mutant_job_id"].astype(str)
        if "mutant_job_id" in data.columns
        else ""
    )

    out["foldx_structure_branch"] = "pdb_structure"
    out["foldx_source_file"] = str(PDB_FOLDX_FILE.relative_to(BASE_DIR))

    optional_columns = [
        "organism",
        "organism_type",
        "gh_family",
        "pdb_id",
        "structure_id",
        "foldx_energy_per_residue",
        "sasa_per_res",
        "hbond_per_res",
        "disulfide_per_res",
    ]

    for column in optional_columns:
        if column in data.columns:
            out[column] = data[column]

    return out[out["candidate_id"].ne("")].copy()


def build_foldx_base():
    foldx_rows = pd.concat(
        [
            standardise_homology_foldx(),
            standardise_pdb_foldx(),
        ],
        ignore_index=True,
        sort=False,
    )

    foldx_rows = foldx_rows[
        foldx_rows["candidate_id"].ne("")
        & foldx_rows["foldx_ddg"].notna()
    ].copy()

    if foldx_rows.empty:
        raise SystemExit("No usable FoldX mutation rows found.")

    best_index = foldx_rows.groupby("candidate_id")["foldx_ddg"].idxmin()
    best = foldx_rows.loc[best_index].copy()

    source_summary = (
        foldx_rows.groupby("candidate_id")
        .agg(
            n_foldx_sources=("foldx_source_file", "count"),
            foldx_ddg_min=("foldx_ddg", "min"),
            foldx_ddg_mean=("foldx_ddg", "mean"),
            foldx_ddg_max=("foldx_ddg", "max"),
            foldx_branches=(
                "foldx_structure_branch",
                lambda x: ";".join(sorted(set(x.dropna().astype(str)))),
            ),
            foldx_source_files=(
                "foldx_source_file",
                lambda x: ";".join(sorted(set(x.dropna().astype(str)))),
            ),
        )
        .reset_index()
    )

    best = best.drop(
        columns=[
            "n_foldx_sources",
            "foldx_ddg_min",
            "foldx_ddg_mean",
            "foldx_ddg_max",
            "foldx_branches",
            "foldx_source_files",
        ],
        errors="ignore",
    )

    base = best.merge(source_summary, on="candidate_id", how="left")

    return base


def build_docking_summary():
    data = read_table(DOCKING_FILE)

    if data.empty:
        return pd.DataFrame(columns=["candidate_id"])

    data["protein"] = data["uniprot_accession"].apply(clean_protein)
    data["mutation"] = data["foldx_mutation_code"].apply(clean_mutation)
    data["candidate_id"] = [
        make_candidate_id(protein, mutation)
        for protein, mutation in zip(data["protein"], data["mutation"])
    ]

    data["delta_affinity_mut_minus_wt"] = pd.to_numeric(
        data["delta_affinity_mut_minus_wt"],
        errors="coerce",
    )

    data["mutant_affinity_kcal_mol"] = pd.to_numeric(
        data["mutant_affinity_kcal_mol"],
        errors="coerce",
    )

    data["wt_affinity_kcal_mol"] = pd.to_numeric(
        data["wt_affinity_kcal_mol"],
        errors="coerce",
    )

    data["docking_improved_vs_wt"] = (
        data["docking_improved_vs_wt"]
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )

    data["docking_retained_vs_wt"] = (
        data["delta_affinity_mut_minus_wt"] <= 0.25
    )

    summary = (
        data.groupby("candidate_id")
        .agg(
            n_docking_rows=("candidate_id", "count"),
            docking_ligands=(
                "ligand_name",
                lambda x: ";".join(sorted(set(x.dropna().astype(str)))),
            ),
            docking_delta_mean=(
                "delta_affinity_mut_minus_wt",
                "mean",
            ),
            docking_delta_best=(
                "delta_affinity_mut_minus_wt",
                "min",
            ),
            docking_delta_worst=(
                "delta_affinity_mut_minus_wt",
                "max",
            ),
            mutant_affinity_mean=("mutant_affinity_kcal_mol", "mean"),
            wt_affinity_mean=("wt_affinity_kcal_mol", "mean"),
            docking_improved_fraction=("docking_improved_vs_wt", "mean"),
            docking_retained_fraction=("docking_retained_vs_wt", "mean"),
        )
        .reset_index()
    )

    return summary


def build_integrated_summary():
    """Load final integrated mutation ranking."""
    data = read_table(INTEGRATED_FILE)

    if data.empty:
        return pd.DataFrame(columns=["candidate_id"])

    data["protein_clean"] = data["protein"].apply(clean_protein)
    data["mutation_clean"] = data["mutation"].apply(clean_mutation)
    data["candidate_id"] = [
        make_candidate_id(protein, mutation)
        for protein, mutation in zip(data["protein_clean"], data["mutation_clean"])
    ]

    keep_columns = [
        "candidate_id",
        "final_rank",
        "final_integrated_score",
        "candidate_class",
        "docking_retention_fraction",
        "mean_delta_binding",
        "md_373_delta_rmsd",
        "md_373_delta_rmsf",
        "md_373_delta_rg",
        "md_373_delta_hbond",
        "predicted_foldx_energy_per_residue",
        "foldx_component",
        "docking_component",
        "md_component",
        "ml_component",
    ]

    rename = {
        column: f"integrated_{column}"
        for column in keep_columns
        if column != "candidate_id"
    }

    return data[keep_columns].rename(columns=rename)


def build_evidence_summary():
    """Load evidence-tiered PDB/homology summary."""
    data = read_table(EVIDENCE_FILE)

    if data.empty:
        return pd.DataFrame(columns=["candidate_id"])

    protein = pd.Series("", index=data.index, dtype=object)
    mutation = pd.Series("", index=data.index, dtype=object)

    for column in ["protein", "candidate_accession", "uniprot_accession", "std_uniprot"]:
        if column in data.columns:
            mask = protein.eq("") & data[column].notna()
            protein.loc[mask] = data.loc[mask, column].apply(clean_protein)

    for column in ["mutation", "foldx_mutation_code", "std_mutation_x", "std_mutation_y"]:
        if column in data.columns:
            mask = mutation.eq("") & data[column].notna()
            mutation.loc[mask] = data.loc[mask, column].apply(clean_mutation)

    data["candidate_id_rebuilt"] = [
        make_candidate_id(protein_value, mutation_value)
        for protein_value, mutation_value in zip(protein, mutation)
    ]

    if "candidate_id" in data.columns:
        data["candidate_id_final"] = data["candidate_id"].fillna("").astype(str)
        data.loc[
            data["candidate_id_final"].isin(["", "nan", "None"]),
            "candidate_id_final",
        ] = data.loc[
            data["candidate_id_final"].isin(["", "nan", "None"]),
            "candidate_id_rebuilt",
        ]
    else:
        data["candidate_id_final"] = data["candidate_id_rebuilt"]

    data["candidate_id_final"] = data["candidate_id_final"].astype(str)

    keep_columns = [
        "candidate_id_final",
        "candidate_tier",
        "integrated_rank_score",
        "evidence_tier",
        "evidence_tier_order",
        "source_branch",
        "evidence_count",
        "score_foldx",
        "score_docking",
        "score_ml",
        "branch_integrated_score",
        "mean_delta_affinity",
        "delta_binding_mut_minus_wt_mean",
        "delta_binding_mut_minus_wt_best",
        "delta_binding_mut_minus_wt_worst",
        "both_ligands_improved",
        "any_ligand_improved",
        "docking_improved_or_retained_all_ligands",
        "docking_improved_or_retained_any_ligand",
        "foldx_predicted_stabilising",
        "candidate_organism",
        "candidate_organism_type",
        "candidate_gh_family",
        "candidate_enzyme_label",
        "candidate_enzyme_display",
    ]

    keep_columns = [column for column in keep_columns if column in data.columns]

    evidence = data[keep_columns].copy()
    evidence = evidence.rename(columns={"candidate_id_final": "candidate_id"})

    rename = {
        column: f"evidence_{column}"
        for column in evidence.columns
        if column != "candidate_id"
    }

    evidence = evidence.rename(columns=rename)

    if "evidence_evidence_tier_order" in evidence.columns:
        evidence["evidence_evidence_tier_order"] = pd.to_numeric(
            evidence["evidence_evidence_tier_order"],
            errors="coerce",
        )

        evidence = evidence.sort_values(
            [
                "evidence_evidence_tier_order",
                "evidence_integrated_rank_score",
            ],
            ascending=[True, False],
        )
    elif "evidence_integrated_rank_score" in evidence.columns:
        evidence = evidence.sort_values(
            "evidence_integrated_rank_score",
            ascending=False,
        )

    evidence = evidence.drop_duplicates("candidate_id", keep="first")

    return evidence[evidence["candidate_id"].ne("")].copy()


def add_interpretation_columns(panel):
    """Add flags, classes, and ranking score."""
    panel["foldx_predicted_stabilising"] = panel["foldx_ddg"] < 0
    panel["strong_foldx_stabilising"] = panel["foldx_ddg"] <= -1.0

    panel["has_docking_proxy"] = (
        panel.get("n_docking_rows", pd.Series(0, index=panel.index))
        .fillna(0)
        .gt(0)
    )

    if "evidence_mean_delta_affinity" in panel.columns:
        panel["has_docking_proxy"] = panel["has_docking_proxy"] | panel[
            "evidence_mean_delta_affinity"
        ].notna()

    if "integrated_mean_delta_binding" in panel.columns:
        panel["has_docking_proxy"] = panel["has_docking_proxy"] | panel[
            "integrated_mean_delta_binding"
        ].notna()

    panel["has_integrated_ranking"] = panel[
        "integrated_final_integrated_score"
    ].notna()

    md_columns = [
        "integrated_md_373_delta_rmsd",
        "integrated_md_373_delta_rmsf",
        "integrated_md_373_delta_rg",
        "integrated_md_373_delta_hbond",
    ]

    panel["has_md_validation"] = panel[
        [column for column in md_columns if column in panel.columns]
    ].notna().any(axis=1)

    docking_retained = pd.Series(False, index=panel.index)

    if "docking_retained_fraction" in panel.columns:
        docking_retained = docking_retained | (
            pd.to_numeric(panel["docking_retained_fraction"], errors="coerce")
            .fillna(0)
            .gt(0)
        )

    if "integrated_docking_retention_fraction" in panel.columns:
        docking_retained = docking_retained | (
            pd.to_numeric(
                panel["integrated_docking_retention_fraction"],
                errors="coerce",
            )
            .fillna(0)
            .gt(0)
        )

    if "evidence_docking_improved_or_retained_any_ligand" in panel.columns:
        docking_retained = docking_retained | (
            panel["evidence_docking_improved_or_retained_any_ligand"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes"])
        )

    panel["docking_retained_or_improved"] = docking_retained

    panel["activity_proxy_status"] = "not_tested"
    panel.loc[
        panel["has_docking_proxy"] & panel["docking_retained_or_improved"],
        "activity_proxy_status",
    ] = "retained_or_improved"
    panel.loc[
        panel["has_docking_proxy"] & ~panel["docking_retained_or_improved"],
        "activity_proxy_status",
    ] = "possible_binding_penalty"

    panel["recommendation_class"] = "not_prioritised"

    panel.loc[
        panel["foldx_predicted_stabilising"],
        "recommendation_class",
    ] = "foldx_stabilising_only"

    panel.loc[
        panel["strong_foldx_stabilising"] & ~panel["has_docking_proxy"],
        "recommendation_class",
    ] = "strong_foldx_only_no_docking"

    panel.loc[
        panel["strong_foldx_stabilising"]
        & panel["docking_retained_or_improved"],
        "recommendation_class",
    ] = "foldx_and_docking_supported"

    panel.loc[
        panel["has_integrated_ranking"],
        "recommendation_class",
    ] = panel.loc[
        panel["has_integrated_ranking"],
        "integrated_candidate_class",
    ].fillna("integrated_ranked_candidate")

    panel["foldx_component_goal"] = normalise_lower_is_better(panel["foldx_ddg"])

    docking_delta = pd.Series(np.nan, index=panel.index, dtype=float)

    for column in [
        "docking_delta_mean",
        "evidence_mean_delta_affinity",
        "evidence_delta_binding_mut_minus_wt_mean",
        "integrated_mean_delta_binding",
    ]:
        if column in panel.columns:
            docking_delta = docking_delta.fillna(
                pd.to_numeric(panel[column], errors="coerce")
            )

    panel["docking_delta_for_score"] = docking_delta
    panel["docking_component_goal"] = normalise_lower_is_better(docking_delta)

    if "integrated_docking_component" in panel.columns:
        panel["docking_component_goal"] = panel[
            "docking_component_goal"
        ].fillna(
            pd.to_numeric(panel["integrated_docking_component"], errors="coerce")
        )

    panel["docking_component_goal"] = panel["docking_component_goal"].fillna(0.0)

    panel["md_component_goal"] = 0.0
    if "integrated_md_component" in panel.columns:
        panel["md_component_goal"] = pd.to_numeric(
            panel["integrated_md_component"],
            errors="coerce",
        ).fillna(0.0)

    panel["computed_priority_score"] = (
        0.75 * panel["foldx_component_goal"].fillna(0.0)
        + 0.25 * panel["docking_component_goal"].fillna(0.0)
    )

    panel["goal_priority_score"] = panel["computed_priority_score"]

    panel.loc[
        panel["has_integrated_ranking"],
        "goal_priority_score",
    ] = panel.loc[
        panel["has_integrated_ranking"],
        "integrated_final_integrated_score",
    ]

    panel = panel.sort_values(
        [
            "goal_priority_score",
            "foldx_ddg",
        ],
        ascending=[False, True],
    ).reset_index(drop=True)

    panel.insert(0, "goal_rank", range(1, len(panel) + 1))

    return panel


def write_report(panel):
    """Write markdown summary report."""
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    top_columns = [
        "goal_rank",
        "protein",
        "mutation",
        "organism_type",
        "gh_family",
        "foldx_ddg",
        "activity_proxy_status",
        "has_docking_proxy",
        "has_md_validation",
        "goal_priority_score",
        "recommendation_class",
    ]

    top_columns = [column for column in top_columns if column in panel.columns]

    lines = [
        "# Mutation Prioritisation Panel",
        "",
        "This panel ranks mutations using computational evidence only.",
        "",
        "FoldX ΔΔG is treated as the thermal-stability predictor.",
        "Docking change is treated as a substrate-binding/activity-related proxy.",
        "MD metrics are treated as validation evidence where available.",
        "",
        f"Rows: {len(panel)}",
        f"Unique proteins: {panel['protein'].nunique()}",
        f"FoldX-predicted stabilising rows: {panel['foldx_predicted_stabilising'].sum()}",
        f"Strong FoldX stabilising rows, ΔΔG <= -1 kcal/mol: {panel['strong_foldx_stabilising'].sum()}",
        f"Rows with docking proxy: {panel['has_docking_proxy'].sum()}",
        f"Rows with integrated ranking: {panel['has_integrated_ranking'].sum()}",
        f"Rows with MD validation metrics: {panel['has_md_validation'].sum()}",
        "",
        "## Recommendation classes",
        "",
        panel["recommendation_class"].value_counts(dropna=False).to_string(),
        "",
        "## Activity-proxy status",
        "",
        panel["activity_proxy_status"].value_counts(dropna=False).to_string(),
        "",
        "## Organism type",
        "",
        panel["organism_type"].value_counts(dropna=False).to_string()
        if "organism_type" in panel.columns
        else "Not available.",
        "",
        "## GH family",
        "",
        panel["gh_family"].value_counts(dropna=False).to_string()
        if "gh_family" in panel.columns
        else "Not available.",
        "",
        "## Top 20 mutation candidates",
        "",
        panel[top_columns].head(20).to_string(index=False),
        "",
        "## Thesis wording",
        "",
        "These candidates should be described as mutations predicted to improve thermal stability while retaining or improving substrate-binding compatibility. They should not be described as experimentally proven thermostable or experimentally more active.",
        "",
    ]

    REPORT_FILE.write_text("\n".join(lines) + "\n")


def main():
    """Build mutation prioritisation panel."""
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    TOP15_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    foldx_base = build_foldx_base()
    docking = build_docking_summary()
    integrated = build_integrated_summary()
    evidence = build_evidence_summary()

    panel = foldx_base.merge(docking, on="candidate_id", how="left")
    panel = panel.merge(integrated, on="candidate_id", how="left")
    panel = panel.merge(evidence, on="candidate_id", how="left")

    for column in ["organism", "organism_type", "gh_family"]:
        evidence_column = f"evidence_candidate_{column}"

        if evidence_column in panel.columns:
            panel[column] = panel[column].fillna(panel[evidence_column])

    panel = add_interpretation_columns(panel)

    panel.to_csv(OUTPUT_FILE, index=False)
    panel.head(15).to_csv(TOP15_FILE, index=False)

    write_report(panel)

    print(f"Wrote: {OUTPUT_FILE}")
    print(f"Wrote: {TOP15_FILE}")
    print(f"Wrote: {REPORT_FILE}")
    print()
    print(f"Rows: {len(panel)}")
    print(f"FoldX stabilising: {panel['foldx_predicted_stabilising'].sum()}")
    print(f"Strong FoldX stabilising: {panel['strong_foldx_stabilising'].sum()}")
    print(f"With docking proxy: {panel['has_docking_proxy'].sum()}")
    print(f"With integrated ranking: {panel['has_integrated_ranking'].sum()}")
    print(f"With MD validation: {panel['has_md_validation'].sum()}")
    print()
    print("Top 15:")
    print(
        panel[
            [
                "goal_rank",
                "protein",
                "mutation",
                "organism_type",
                "gh_family",
                "foldx_ddg",
                "activity_proxy_status",
                "has_md_validation",
                "goal_priority_score",
                "recommendation_class",
            ]
        ]
        .head(15)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
