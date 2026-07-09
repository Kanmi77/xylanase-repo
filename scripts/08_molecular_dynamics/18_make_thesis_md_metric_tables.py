#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import pandas as pd


def mean_sd(mean, sd, digits=3):
    if pd.isna(mean):
        return ""
    if pd.isna(sd):
        return f"{mean:.{digits}f}"
    return f"{mean:.{digits}f} ± {sd:.{digits}f}"


def mutation_from_system_id(system_id, state):
    if str(state).upper() != "MUT":
        return "—"
    m = re.search(r"_MUT_([^_]+)_", str(system_id))
    return m.group(1) if m else "MUT"


def make_phase2_table(phase2_wide, outdir):
    df = pd.read_csv(phase2_wide)

    rows = []
    for _, r in df.iterrows():
        rows.append({
            "Protein": r.get("protein", ""),
            "System ID": r.get("system_id", ""),
            "Variant": r.get("state", ""),
            "Mutation": mutation_from_system_id(r.get("system_id", ""), r.get("state", "")),
            "Temperature": str(r.get("temperature", "")).replace("K", " K"),
            "RMSD mean ± SD (nm)": mean_sd(pd.to_numeric(r.get("rmsd_mean_all"), errors="coerce"),
                                           pd.to_numeric(r.get("rmsd_sd_all"), errors="coerce")),
            "RMSF mean ± SD (nm)": mean_sd(pd.to_numeric(r.get("rmsf_mean_all"), errors="coerce"),
                                           pd.to_numeric(r.get("rmsf_sd_all"), errors="coerce")),
            "Rg mean ± SD (nm)": mean_sd(pd.to_numeric(r.get("rg_mean_all"), errors="coerce"),
                                         pd.to_numeric(r.get("rg_sd_all"), errors="coerce")),
            "H-bonds mean ± SD": mean_sd(pd.to_numeric(r.get("hbond_mean_all"), errors="coerce"),
                                         pd.to_numeric(r.get("hbond_sd_all"), errors="coerce"), digits=1),
            "Late-50 RMSD mean (nm)": round(pd.to_numeric(r.get("rmsd_late50_mean"), errors="coerce"), 3),
            "Late-50 Rg mean (nm)": round(pd.to_numeric(r.get("rg_late50_mean"), errors="coerce"), 3),
            "RMSF max residue/index": r.get("rmsf_max_residue_or_index", ""),
            "RMSF max value (nm)": round(pd.to_numeric(r.get("rmsf_max_rmsf_value"), errors="coerce"), 3),
            "Interpretation": ""
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(["Protein", "Temperature", "Variant"])

    out_csv = outdir / "table_phase2_wt_mutant_md_metrics_compact.csv"
    out_md = outdir / "table_phase2_wt_mutant_md_metrics_compact.md"

    out.to_csv(out_csv, index=False)
    out.to_markdown(out_md, index=False)

    return out_csv, out_md


def make_phase1_table(phase1_long, outdir):
    df = pd.read_csv(phase1_long)

    keep_metrics = ["rmsd", "rmsf", "rg", "hbond", "sasa"]
    df = df[df["metric"].isin(keep_metrics)].copy()

    index_cols = ["protein", "system_id", "state", "temperature"]

    wide = df.pivot_table(
        index=index_cols,
        columns="metric",
        values=["mean_all", "sd_all", "late50_mean", "max_rmsf_residue_or_index", "max_rmsf_value"],
        aggfunc="first"
    )

    wide.columns = [f"{metric}_{value}" for value, metric in wide.columns]
    wide = wide.reset_index()

    rows = []
    for _, r in wide.iterrows():
        rows.append({
            "Enzyme / PDB ID": r.get("protein", ""),
            "System ID": r.get("system_id", ""),
            "State": r.get("state", ""),
            "Temperature": str(r.get("temperature", "")).replace("K", " K"),
            "RMSD mean ± SD (nm)": mean_sd(pd.to_numeric(r.get("rmsd_mean_all"), errors="coerce"),
                                           pd.to_numeric(r.get("rmsd_sd_all"), errors="coerce")),
            "RMSF mean ± SD (nm)": mean_sd(pd.to_numeric(r.get("rmsf_mean_all"), errors="coerce"),
                                           pd.to_numeric(r.get("rmsf_sd_all"), errors="coerce")),
            "Rg mean ± SD (nm)": mean_sd(pd.to_numeric(r.get("rg_mean_all"), errors="coerce"),
                                         pd.to_numeric(r.get("rg_sd_all"), errors="coerce")),
            "H-bonds mean ± SD": mean_sd(pd.to_numeric(r.get("hbond_mean_all"), errors="coerce"),
                                         pd.to_numeric(r.get("hbond_sd_all"), errors="coerce"), digits=1),
            "SASA mean ± SD (nm²)": mean_sd(pd.to_numeric(r.get("sasa_mean_all"), errors="coerce"),
                                            pd.to_numeric(r.get("sasa_sd_all"), errors="coerce")),
            "Late-50 RMSD mean (nm)": round(pd.to_numeric(r.get("rmsd_late50_mean"), errors="coerce"), 3),
            "Late-50 Rg mean (nm)": round(pd.to_numeric(r.get("rg_late50_mean"), errors="coerce"), 3),
            "RMSF max residue/index": r.get("rmsf_max_rmsf_residue_or_index", ""),
            "RMSF max value (nm)": round(pd.to_numeric(r.get("rmsf_max_max_rmsf_value"), errors="coerce"), 3),
            "Interpretation": ""
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(["Enzyme / PDB ID", "Temperature"])

    out_csv = outdir / "table_phase1_pdb_md_metrics_compact.csv"
    out_md = outdir / "table_phase1_pdb_md_metrics_compact.md"

    out.to_csv(out_csv, index=False)
    out.to_markdown(out_md, index=False)

    return out_csv, out_md


def make_phase2_delta_table(phase2_wide, outdir):
    df = pd.read_csv(phase2_wide)

    metrics = {
        "rmsd_mean_all": "ΔRMSD MUT−WT (nm)",
        "rmsf_mean_all": "ΔRMSF MUT−WT (nm)",
        "rg_mean_all": "ΔRg MUT−WT (nm)",
        "hbond_mean_all": "ΔH-bonds MUT−WT"
    }

    needed = ["protein", "temperature", "state", "system_id"] + list(metrics.keys())
    df = df[[c for c in needed if c in df.columns]].copy()

    wt = df[df["state"].astype(str).str.upper() == "WT"].copy()
    mut = df[df["state"].astype(str).str.upper() == "MUT"].copy()

    merged = mut.merge(
        wt,
        on=["protein", "temperature"],
        suffixes=("_mut", "_wt"),
        how="inner"
    )

    rows = []
    for _, r in merged.iterrows():
        row = {
            "Protein": r["protein"],
            "Mutation": mutation_from_system_id(r.get("system_id_mut", ""), "MUT"),
            "Temperature": str(r["temperature"]).replace("K", " K"),
        }

        for col, label in metrics.items():
            mut_v = pd.to_numeric(r.get(f"{col}_mut"), errors="coerce")
            wt_v = pd.to_numeric(r.get(f"{col}_wt"), errors="coerce")
            row[label] = round(mut_v - wt_v, 4) if pd.notna(mut_v) and pd.notna(wt_v) else ""

        row["Interpretation"] = ""
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(["Protein", "Temperature"])

    out_csv = outdir / "table_phase2_mutant_minus_wt_md_delta.csv"
    out_md = outdir / "table_phase2_mutant_minus_wt_md_delta.md"

    out.to_csv(out_csv, index=False)
    out.to_markdown(out_md, index=False)

    return out_csv, out_md


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase1-long", required=True)
    parser.add_argument("--phase2-wide", required=True)
    parser.add_argument("--outdir", default="results/md_tables")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    p1_csv, p1_md = make_phase1_table(args.phase1_long, outdir)
    p2_csv, p2_md = make_phase2_table(args.phase2_wide, outdir)
    d_csv, d_md = make_phase2_delta_table(args.phase2_wide, outdir)

    print("Saved:")
    print(p1_csv)
    print(p1_md)
    print(p2_csv)
    print(p2_md)
    print(d_csv)
    print(d_md)


if __name__ == "__main__":
    main()
