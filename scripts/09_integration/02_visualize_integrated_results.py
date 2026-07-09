#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

BASE = Path.home() / "xylanase-thesis"

INTEGRATION = BASE / "results/integration"
FIGDIR = BASE / "results/figures/integration"
FIGDIR.mkdir(parents=True, exist_ok=True)

RANKING = INTEGRATION / "final_integrated_candidate_ranking.csv"
FULL = INTEGRATION / "integrated_thermostability_catalytic_retention_scores.csv"
GROUP = INTEGRATION / "integrated_group_summary.csv"

OUT_TOP15 = FIGDIR / "top15_integrated_candidates.png"
OUT_COMPONENTS = FIGDIR / "top10_score_components.png"
OUT_GROUP = FIGDIR / "best_score_by_group.png"
OUT_MD = FIGDIR / "md_373_delta_rmsd_vs_rmsf.png"
OUT_FOLDX_MD = FIGDIR / "foldx_ddg_vs_md_373_delta_rmsd.png"


def require(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")


def candidate_label(row):
    protein = str(row.get("protein", ""))
    mutation = str(row.get("mutation", ""))

    if mutation and mutation.lower() != "nan":
        return f"{protein}_{mutation}"

    return protein


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)


def plot_top15(df):
    d = df.copy()
    d["candidate"] = d.apply(candidate_label, axis=1)
    d = d.sort_values("final_integrated_score", ascending=False).head(15)
    d = d.sort_values("final_integrated_score", ascending=True)

    plt.figure(figsize=(9, 7))
    plt.barh(d["candidate"], d["final_integrated_score"])
    plt.xlabel("Final integrated score")
    plt.ylabel("Candidate")
    plt.title("Top 15 integrated thermostability and catalytic-retention candidates")
    savefig(OUT_TOP15)


def plot_components(df):
    d = df.copy()
    d["candidate"] = d.apply(candidate_label, axis=1)
    d = d.sort_values("final_integrated_score", ascending=False).head(10)

    component_cols = [
        "foldx_component",
        "md_component",
        "docking_component",
        "ml_component",
    ]

    component_cols = [c for c in component_cols if c in d.columns]

    if not component_cols:
        print("No component columns found; skipping component plot.")
        return

    plot_df = d.set_index("candidate")[component_cols]

    plt.figure(figsize=(11, 6))
    bottom = np.zeros(len(plot_df))

    for col in component_cols:
        values = plot_df[col].fillna(0).values
        plt.bar(plot_df.index, values, bottom=bottom, label=col.replace("_component", ""))
        bottom += values

    plt.ylabel("Component score")
    plt.xlabel("Candidate")
    plt.title("Score-component profile of top 10 integrated candidates")
    plt.xticks(rotation=45, ha="right")
    plt.legend()
    savefig(OUT_COMPONENTS)


def plot_group_summary(df, group):
    if group.exists():
        g = pd.read_csv(group)
    else:
        g = (
            df.groupby(["organism_type", "gh_family"], dropna=False)
            .agg(best_score=("final_integrated_score", "max"))
            .reset_index()
        )

    if "best_score" not in g.columns:
        if "final_integrated_score" in g.columns:
            g["best_score"] = g["final_integrated_score"]
        else:
            print("No best_score column found; skipping group plot.")
            return

    g["group"] = g["organism_type"].astype(str) + "_" + g["gh_family"].astype(str)
    g = g.sort_values("best_score", ascending=True)

    plt.figure(figsize=(8, 5))
    plt.barh(g["group"], g["best_score"])
    plt.xlabel("Best final integrated score")
    plt.ylabel("Group")
    plt.title("Best integrated candidate score by organism type and GH family")
    savefig(OUT_GROUP)


def plot_md_373(df):
    needed = [
        "md_373_delta_rmsd",
        "md_373_delta_rmsf",
        "protein",
        "mutation",
    ]

    if not all(c in df.columns for c in needed):
        print("MD 373 K columns not found; skipping MD scatter.")
        return

    d = df.dropna(subset=["md_373_delta_rmsd", "md_373_delta_rmsf"]).copy()

    if d.empty:
        print("No MD 373 K values found; skipping MD scatter.")
        return

    d["candidate"] = d.apply(candidate_label, axis=1)

    plt.figure(figsize=(8, 6))
    plt.scatter(d["md_373_delta_rmsd"], d["md_373_delta_rmsf"])

    for _, r in d.iterrows():
        plt.text(
            r["md_373_delta_rmsd"],
            r["md_373_delta_rmsf"],
            r["candidate"],
            fontsize=8
        )

    plt.axvline(0, linestyle="--", linewidth=1)
    plt.axhline(0, linestyle="--", linewidth=1)

    plt.xlabel("ΔRMSD at 373 K, mutant − WT (nm)")
    plt.ylabel("ΔRMSF at 373 K, mutant − WT (nm)")
    plt.title("MD high-temperature structural response at 373 K")
    savefig(OUT_MD)


def plot_foldx_vs_md(df):
    needed = [
        "foldx_ddg",
        "md_373_delta_rmsd",
        "protein",
        "mutation",
    ]

    if not all(c in df.columns for c in needed):
        print("FoldX/MD columns not found; skipping FoldX vs MD plot.")
        return

    d = df.dropna(subset=["foldx_ddg", "md_373_delta_rmsd"]).copy()

    if d.empty:
        print("No FoldX/MD paired values found; skipping FoldX vs MD plot.")
        return

    d["candidate"] = d.apply(candidate_label, axis=1)

    plt.figure(figsize=(8, 6))
    plt.scatter(d["foldx_ddg"], d["md_373_delta_rmsd"])

    for _, r in d.iterrows():
        plt.text(
            r["foldx_ddg"],
            r["md_373_delta_rmsd"],
            r["candidate"],
            fontsize=8
        )

    plt.axvline(0, linestyle="--", linewidth=1)
    plt.axhline(0, linestyle="--", linewidth=1)

    plt.xlabel("FoldX ΔΔG (kcal/mol)")
    plt.ylabel("ΔRMSD at 373 K, mutant − WT (nm)")
    plt.title("FoldX mutation stability vs MD high-temperature RMSD response")
    savefig(OUT_FOLDX_MD)


def main():
    require(RANKING)

    df = pd.read_csv(RANKING)

    print("Loaded:", RANKING)
    print("Rows:", len(df))
    print("Columns:", ", ".join(df.columns))

    plot_top15(df)
    plot_components(df)
    plot_group_summary(df, GROUP)
    plot_md_373(df)
    plot_foldx_vs_md(df)

    print("\nDone. Figures saved in:")
    print(FIGDIR)


if __name__ == "__main__":
    main()
