#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "xylanase-thesis"

FOLDX = BASE / "results/foldx_clean/tier2_ddg_ranked_annotated.csv"
DOCK = BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_wt_mutant_docking_comparison.csv"
MD = BASE / "md_tier2_wt_mutant_compact/analysis/compact_md_wt_mutant_comparison.csv"
ML = BASE / "results/ml/structural_stability_ml_predictions.csv"

OUT_DIR = BASE / "results/integration"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FULL = OUT_DIR / "integrated_thermostability_catalytic_retention_scores.csv"
OUT_RANKED = OUT_DIR / "final_integrated_candidate_ranking.csv"
OUT_GROUP = OUT_DIR / "integrated_group_summary.csv"


def require_file(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path}")


def minmax_score(series, higher_is_better=True):
    x = pd.to_numeric(series, errors="coerce")
    lo = x.min()
    hi = x.max()

    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(0.5, index=series.index)

    score = (x - lo) / (hi - lo)

    if not higher_is_better:
        score = 1 - score

    return score.clip(0, 1)


def load_foldx():
    require_file(FOLDX)
    df = pd.read_csv(FOLDX)

    # Keep best mutation per protein according to FoldX ΔΔG
    df = df.sort_values("ddg").groupby("protein", as_index=False).first()

    keep = [
        "protein", "mutation", "ddg", "organism", "organism_type", "gh_family",
        "foldx_energy_per_residue", "sasa_per_res", "hbond_per_res",
        "disulfide_per_res"
    ]

    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    df = df.rename(columns={"ddg": "foldx_ddg"})

    # More negative ΔΔG is better, so -ΔΔG is scored higher.
    df["foldx_stability_signal"] = -pd.to_numeric(df["foldx_ddg"], errors="coerce")

    return df


def load_docking():
    require_file(DOCK)
    df = pd.read_csv(DOCK)

    # Aggregate across xylobiose and xylotriose
    agg = (
        df.groupby(["protein", "mutation"], as_index=False)
        .agg(
            docking_ligands=("ligand", "nunique"),
            mean_wt_binding=("wt_binding_energy", "mean"),
            mean_mut_binding=("mut_binding_energy", "mean"),
            mean_delta_binding=("delta_binding_mut_minus_wt", "mean"),
            best_delta_binding=("delta_binding_mut_minus_wt", "min"),
            worst_delta_binding=("delta_binding_mut_minus_wt", "max"),
            retained_or_improved_ligands=(
                "functional_integrity",
                lambda x: ((x == "retained_binding") | (x == "improved_binding") | (x == "retained_or_improved")).sum()
            ),
            improved_ligands=(
                "functional_integrity",
                lambda x: (x == "improved_binding").sum()
            ),
            weakened_ligands=(
                "functional_integrity",
                lambda x: (x == "weakened_binding").sum()
            ),
        )
    )

    agg["docking_retention_fraction"] = agg["retained_or_improved_ligands"] / agg["docking_ligands"]

    # More negative delta means mutant binds better. Retention fraction is also important.
    agg["docking_binding_signal"] = -pd.to_numeric(agg["mean_delta_binding"], errors="coerce")

    return agg


def load_md():
    require_file(MD)
    df = pd.read_csv(MD)

    # Create stability signals.
    # Lower RMSD/RMSF/Rg change is better. Higher H-bond change is better.
    df["md_rmsd_signal"] = -pd.to_numeric(df["delta_mut_minus_wt_rmsd_late_mean_nm"], errors="coerce")
    df["md_rmsf_signal"] = -pd.to_numeric(df["delta_mut_minus_wt_rmsf_late_mean_nm"], errors="coerce")
    df["md_rg_signal"] = -pd.to_numeric(df["delta_mut_minus_wt_rg_late_mean_nm"], errors="coerce")
    df["md_hbond_signal"] = pd.to_numeric(df["delta_mut_minus_wt_hbond_late_mean"], errors="coerce")

    # Emphasize 373 K because it is the high-temperature stress condition.
    df["temperature_weight"] = np.where(df["temperature_K"] == 373, 2.0, 1.0)

    md_metric_cols = ["md_rmsd_signal", "md_rmsf_signal", "md_rg_signal", "md_hbond_signal"]

    rows = []

    for (protein, mutation), sub in df.groupby(["protein", "mutation"]):
        out = {
            "protein": protein,
            "mutation": mutation,
            "md_temperatures": ",".join(str(x) for x in sorted(sub["temperature_K"].unique())),
            "md_n_temperatures": sub["temperature_K"].nunique(),
        }

        for col in md_metric_cols:
            values = pd.to_numeric(sub[col], errors="coerce")
            weights = pd.to_numeric(sub["temperature_weight"], errors="coerce")
            ok = values.notna() & weights.notna()

            if ok.any():
                out[f"{col}_weighted"] = np.average(values[ok], weights=weights[ok])
            else:
                out[f"{col}_weighted"] = np.nan

        # Store explicit 373K values if available
        t373 = sub[sub["temperature_K"] == 373]
        if not t373.empty:
            r = t373.iloc[0]
            out["md_373_delta_rmsd"] = r.get("delta_mut_minus_wt_rmsd_late_mean_nm", np.nan)
            out["md_373_delta_rmsf"] = r.get("delta_mut_minus_wt_rmsf_late_mean_nm", np.nan)
            out["md_373_delta_rg"] = r.get("delta_mut_minus_wt_rg_late_mean_nm", np.nan)
            out["md_373_delta_hbond"] = r.get("delta_mut_minus_wt_hbond_late_mean", np.nan)
            out["md_373_rmsd_interpretation"] = r.get("rmsd_interpretation", "")
            out["md_373_rmsf_interpretation"] = r.get("rmsf_interpretation", "")
            out["md_373_hbond_interpretation"] = r.get("hbond_interpretation", "")

        rows.append(out)

    return pd.DataFrame(rows)


def load_ml():
    require_file(ML)
    df = pd.read_csv(ML)

    keep = [
        "uniprot_accession", "foldx_energy_per_residue",
        "predicted_foldx_energy_per_residue", "prediction_error"
    ]

    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    df = df.rename(columns={"uniprot_accession": "protein"})

    # Lower predicted FoldX energy per residue is better.
    if "predicted_foldx_energy_per_residue" in df.columns:
        df["ml_stability_signal"] = -pd.to_numeric(df["predicted_foldx_energy_per_residue"], errors="coerce")
    else:
        df["ml_stability_signal"] = np.nan

    return df


def classify_candidate(row):
    score = row.get("final_integrated_score", np.nan)
    foldx = row.get("foldx_ddg", np.nan)
    dock_frac = row.get("docking_retention_fraction", np.nan)
    md_rmsd_373 = row.get("md_373_delta_rmsd", np.nan)
    md_rmsf_373 = row.get("md_373_delta_rmsf", np.nan)

    if pd.isna(score):
        return "insufficient_data"

    if (
        foldx < 0
        and dock_frac >= 0.5
        and pd.notna(md_rmsd_373)
        and pd.notna(md_rmsf_373)
        and md_rmsd_373 < 0
        and md_rmsf_373 < 0
    ):
        return "high_confidence_stabilizing_and_functionally_retained"

    if foldx < 0 and dock_frac >= 0.5 and pd.notna(md_rmsd_373) and md_rmsd_373 < 0:
        return "moderate_confidence_stabilizing"

    if foldx < 0 and dock_frac >= 0.5:
        return "foldx_and_docking_supported"

    if foldx < 0:
        return "foldx_only_supported"

    return "weak_or_conflicting_support"


def main():
    foldx = load_foldx()
    dock = load_docking()
    md = load_md()
    ml = load_ml()

    df = foldx.merge(dock, on=["protein", "mutation"], how="left")
    df = df.merge(md, on=["protein", "mutation"], how="left")
    df = df.merge(ml, on="protein", how="left", suffixes=("", "_ml"))

    # Component scores
    df["score_foldx"] = minmax_score(df["foldx_stability_signal"], higher_is_better=True)
    df["score_docking_energy"] = minmax_score(df["docking_binding_signal"], higher_is_better=True)
    df["score_docking_retention"] = pd.to_numeric(df["docking_retention_fraction"], errors="coerce").fillna(0)

    df["score_md_rmsd"] = minmax_score(df["md_rmsd_signal_weighted"], higher_is_better=True)
    df["score_md_rmsf"] = minmax_score(df["md_rmsf_signal_weighted"], higher_is_better=True)
    df["score_md_rg"] = minmax_score(df["md_rg_signal_weighted"], higher_is_better=True)
    df["score_md_hbond"] = minmax_score(df["md_hbond_signal_weighted"], higher_is_better=True)

    df["score_ml"] = minmax_score(df["ml_stability_signal"], higher_is_better=True)

    # Main component blocks
    df["foldx_component"] = df["score_foldx"]
    df["docking_component"] = (
        0.6 * df["score_docking_retention"] +
        0.4 * df["score_docking_energy"].fillna(0)
    )

    df["md_component"] = (
        0.35 * df["score_md_rmsd"].fillna(0) +
        0.30 * df["score_md_rmsf"].fillna(0) +
        0.15 * df["score_md_rg"].fillna(0) +
        0.20 * df["score_md_hbond"].fillna(0)
    )

    df["ml_component"] = df["score_ml"].fillna(0)

    # Final score:
    # FoldX = thermostability prediction
    # Docking = catalytic retention proxy
    # MD = dynamic structural validation
    # ML = supporting structural-stability prediction
    df["final_integrated_score"] = (
        0.35 * df["foldx_component"] +
        0.25 * df["docking_component"] +
        0.30 * df["md_component"] +
        0.10 * df["ml_component"]
    )

    df["candidate_class"] = df.apply(classify_candidate, axis=1)

    df = df.sort_values("final_integrated_score", ascending=False).reset_index(drop=True)
    df["final_rank"] = range(1, len(df) + 1)

    df.to_csv(OUT_FULL, index=False)

    # Keep concise ranked output for thesis table
    ranked_cols = [
        "final_rank", "protein", "mutation", "organism", "organism_type", "gh_family",
        "final_integrated_score", "candidate_class",
        "foldx_ddg",
        "docking_retention_fraction", "mean_delta_binding",
        "md_373_delta_rmsd", "md_373_delta_rmsf", "md_373_delta_rg", "md_373_delta_hbond",
        "predicted_foldx_energy_per_residue",
        "foldx_component", "docking_component", "md_component", "ml_component"
    ]

    ranked_cols = [c for c in ranked_cols if c in df.columns]
    ranked = df[ranked_cols].copy()
    ranked.to_csv(OUT_RANKED, index=False)

    group = (
        df.groupby(["organism_type", "gh_family"], dropna=False)
        .agg(
            n_candidates=("protein", "nunique"),
            mean_final_score=("final_integrated_score", "mean"),
            best_score=("final_integrated_score", "max"),
            mean_foldx_ddg=("foldx_ddg", "mean"),
            mean_docking_retention=("docking_retention_fraction", "mean"),
            mean_md_component=("md_component", "mean"),
        )
        .reset_index()
        .sort_values("best_score", ascending=False)
    )

    group.to_csv(OUT_GROUP, index=False)

    print("Saved:")
    print(OUT_FULL)
    print(OUT_RANKED)
    print(OUT_GROUP)

    print("\nTop 20 integrated candidates:")
    print(ranked.head(20).to_string(index=False))

    print("\nGroup summary:")
    print(group.to_string(index=False))


if __name__ == "__main__":
    main()
