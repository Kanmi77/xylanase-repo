#!/usr/bin/env python3
# Purpose: Summarise DSSP outputs.

from pathlib import Path
import subprocess
import shutil
import re
import pandas as pd
import numpy as np

BASE = Path.home() / "xylanase-thesis"

ROOT = BASE / "md_tier2_wt_mutant_compact"
SYSTEMS = ROOT / "systems"
ANALYSIS = ROOT / "analysis"
ANALYSIS.mkdir(parents=True, exist_ok=True)

COMPLETION_TABLE = ROOT / "compact_md_deep_completion_summary.csv"
METRICS = ANALYSIS / "compact_md_metrics_summary.csv"

OUT_24 = ANALYSIS / "md_24_simulation_mean_metrics_with_dssp.csv"
OUT_PAIR = ANALYSIS / "md_wt_mut_333K_373K_table_with_dssp.csv"
OUT_DSSP_LONG = ANALYSIS / "md_dssp_secondary_structure_long.csv"

GMX = "gmx"


def run_cmd(cmd, cwd=None, input_text=None, timeout=1800):
    try:
        p = subprocess.run(
            cmd,
            cwd=cwd,
            input=input_text,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "TIMEOUT"


def parse_system_name(name):
    protein = name.split("_")[0]

    state = "MUT" if "_MUT_" in name else "WT"

    mutation = ""
    if state == "MUT":
        m = re.search(r"_MUT_([^_]+)_", name)
        if m:
            mutation = m.group(1)

    organism_type = ""
    if "bacterial" in name:
        organism_type = "bacterial"
    elif "fungal" in name:
        organism_type = "fungal"

    gh_family = ""
    if "GH10" in name:
        gh_family = "GH10"
    elif "GH11" in name:
        gh_family = "GH11"

    label = name.replace(f"{protein}_", "")

    return protein, state, mutation, organism_type, gh_family, label


def read_xvg(path):
    rows = []

    if not path.exists():
        return pd.DataFrame()

    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    rows.append([float(x) for x in parts])
                except Exception:
                    continue

    return pd.DataFrame(rows)


def summarize_xvg(path, value_col=1):
    df = read_xvg(path)

    if df.empty or value_col not in df.columns:
        return {
            "mean": np.nan,
            "sd": np.nan,
            "late_mean": np.nan,
            "late_sd": np.nan,
            "final": np.nan,
            "min": np.nan,
            "max": np.nan,
        }

    time = df[0].astype(float)
    values = df[value_col].astype(float)

    cutoff = time.max() * 0.5
    late = values[time >= cutoff]

    return {
        "mean": values.mean(),
        "sd": values.std(ddof=1),
        "late_mean": late.mean() if len(late) else np.nan,
        "late_sd": late.std(ddof=1) if len(late) > 1 else np.nan,
        "final": values.iloc[-1],
        "min": values.min(),
        "max": values.max(),
    }


def dssp_available():
    exe = shutil.which("mkdssp") or shutil.which("dssp")
    return exe

def run_dssp_for_system(sysdir):
    """
    Robust DSSP runner.

    Tries:
    1. md.xtc with Protein selection
    2. md.xtc with numeric group 1
    3. PBC/fit protein trajectory, then DSSP with Protein selection
    4. PBC/fit protein trajectory, then DSSP with numeric group 1
    """
    exe = dssp_available()
    if exe is None:
        return False, "DSSP executable not found"

    ss_xpm = sysdir / "ss.xpm"
    scount_xvg = sysdir / "scount.xvg"

    if ss_xpm.exists() and scount_xvg.exists() and scount_xvg.stat().st_size > 0:
        return True, "DSSP already exists"

    env_prefix = f"export DSSP={exe}; export GMX_DSSP={exe}; "

    attempts = []

    # Attempt 1: direct, Protein group
    attempts.append((
        "direct Protein selection",
        env_prefix
        + f'printf "Protein\\n" | {GMX} do_dssp '
        + "-s md.tpr "
        + "-f md.xtc "
        + "-o ss.xpm "
        + "-sc scount.xvg"
    ))

    # Attempt 2: direct, group 1
    attempts.append((
        "direct numeric group 1",
        env_prefix
        + f'printf "1\\n" | {GMX} do_dssp '
        + "-s md.tpr "
        + "-f md.xtc "
        + "-o ss.xpm "
        + "-sc scount.xvg"
    ))

    # Attempt 3: create fitted protein trajectory
    protein_fit = sysdir / "protein_fit.xtc"

    if not protein_fit.exists():
        trjconv_cmd = (
            f'printf "Protein\\nProtein\\n" | {GMX} trjconv '
            + "-s md.tpr "
            + "-f md.xtc "
            + "-o protein_fit.xtc "
            + "-fit rot+trans "
            + "-pbc mol "
            + "-center"
        )
        code, out, err = run_cmd(trjconv_cmd, cwd=sysdir, timeout=2400)

    if protein_fit.exists():
        attempts.append((
            "protein_fit.xtc Protein selection",
            env_prefix
            + f'printf "Protein\\n" | {GMX} do_dssp '
            + "-s md.tpr "
            + "-f protein_fit.xtc "
            + "-o ss.xpm "
            + "-sc scount.xvg"
        ))

        attempts.append((
            "protein_fit.xtc numeric group 1",
            env_prefix
            + f'printf "1\\n" | {GMX} do_dssp '
            + "-s md.tpr "
            + "-f protein_fit.xtc "
            + "-o ss.xpm "
            + "-sc scount.xvg"
        ))

    messages = []

    for label, cmd in attempts:
        # Remove failed partial files before each attempt
        for f in [ss_xpm, scount_xvg]:
            if f.exists() and f.stat().st_size == 0:
                f.unlink()

        code, out, err = run_cmd(cmd, cwd=sysdir, timeout=2400)
        combined = (out + "\n" + err)[-3000:]

        messages.append(f"ATTEMPT: {label}\nRETURN_CODE: {code}\n{combined}")

        if (
            code == 0
            and scount_xvg.exists()
            and scount_xvg.stat().st_size > 0
        ):
            return True, f"DSSP completed using: {label}"

    return False, "\n\n".join(messages)[-5000:]
def parse_scount_xvg(path):
    """
    scount.xvg columns vary, but commonly:
    time, structure, coil, B-sheet, B-bridge, bend, turn, A-helix, 3-helix, 5-helix

    We parse legends from @ sN legend lines if present.
    """
    if not path.exists():
        return {}, pd.DataFrame()

    legends = {}
    data = []

    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith("@") and "legend" in line:
                m = re.search(r's(\d+)\s+legend\s+"([^"]+)"', line)
                if m:
                    legends[int(m.group(1)) + 1] = m.group(2)  # +1 because col0 is time
            elif not line.startswith(("#", "@")):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        data.append([float(x) for x in parts])
                    except Exception:
                        pass

    if not data:
        return {}, pd.DataFrame()

    df = pd.DataFrame(data)

    # Rename columns
    rename = {0: "time_ps_or_ns"}
    for col in df.columns:
        if col == 0:
            continue
        rename[col] = legends.get(col, f"ss_col_{col}")

    df = df.rename(columns=rename)

    # Compute means for available DSSP categories
    result = {}

    numeric_cols = [c for c in df.columns if c != "time_ps_or_ns"]

    for col in numeric_cols:
        values = pd.to_numeric(df[col], errors="coerce")
        result[f"dssp_mean_{clean_col(col)}"] = values.mean()
        result[f"dssp_final_{clean_col(col)}"] = values.iloc[-1]

    # Create broad secondary structure summaries.
    helix_cols = [c for c in numeric_cols if "helix" in str(c).lower()]
    sheet_cols = [
        c for c in numeric_cols
        if "sheet" in str(c).lower() or "b-bridge" in str(c).lower() or "bridge" in str(c).lower()
    ]
    turn_bend_cols = [
        c for c in numeric_cols
        if "turn" in str(c).lower() or "bend" in str(c).lower()
    ]
    coil_cols = [c for c in numeric_cols if "coil" in str(c).lower()]

    def sum_mean(cols):
        if not cols:
            return np.nan
        return df[cols].apply(pd.to_numeric, errors="coerce").sum(axis=1).mean()

    def sum_final(cols):
        if not cols:
            return np.nan
        return df[cols].apply(pd.to_numeric, errors="coerce").sum(axis=1).iloc[-1]

    result["dssp_helix_mean_count"] = sum_mean(helix_cols)
    result["dssp_sheet_mean_count"] = sum_mean(sheet_cols)
    result["dssp_turn_bend_mean_count"] = sum_mean(turn_bend_cols)
    result["dssp_coil_mean_count"] = sum_mean(coil_cols)

    result["dssp_helix_final_count"] = sum_final(helix_cols)
    result["dssp_sheet_final_count"] = sum_final(sheet_cols)
    result["dssp_turn_bend_final_count"] = sum_final(turn_bend_cols)
    result["dssp_coil_final_count"] = sum_final(coil_cols)

    result["dssp_total_mean_count"] = sum_mean(numeric_cols)
    result["dssp_total_final_count"] = sum_final(numeric_cols)

    if pd.notna(result["dssp_total_mean_count"]) and result["dssp_total_mean_count"] != 0:
        result["dssp_helix_mean_fraction"] = result["dssp_helix_mean_count"] / result["dssp_total_mean_count"]
        result["dssp_sheet_mean_fraction"] = result["dssp_sheet_mean_count"] / result["dssp_total_mean_count"]
        result["dssp_turn_bend_mean_fraction"] = result["dssp_turn_bend_mean_count"] / result["dssp_total_mean_count"]
        result["dssp_coil_mean_fraction"] = result["dssp_coil_mean_count"] / result["dssp_total_mean_count"]
    else:
        result["dssp_helix_mean_fraction"] = np.nan
        result["dssp_sheet_mean_fraction"] = np.nan
        result["dssp_turn_bend_mean_fraction"] = np.nan
        result["dssp_coil_mean_fraction"] = np.nan

    return result, df


def clean_col(x):
    return (
        str(x)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
    )


def compute_secondary_retention(pair_df):
    """
    Adds WT vs MUT DSSP retention values:
    delta_mut_minus_wt_helix_fraction etc.
    """
    rows = []

    for (protein, temperature), sub in pair_df.groupby(["protein", "temperature_K"]):
        wt = sub[sub["state"] == "WT"]
        mut = sub[sub["state"] == "MUT"]

        if wt.empty or mut.empty:
            continue

        wt = wt.iloc[0]
        mut = mut.iloc[0]

        out = {
            "protein": protein,
            "temperature_K": temperature,
            "mutation": mut.get("mutation", ""),
            "organism_type": mut.get("organism_type", wt.get("organism_type", "")),
            "gh_family": mut.get("gh_family", wt.get("gh_family", "")),
        }

        metrics = [
            "rmsd_mean_nm",
            "rmsd_late_mean_nm",
            "rmsf_mean_nm",
            "rg_mean_nm",
            "rg_late_mean_nm",
            "hbond_mean",
            "hbond_late_mean",
            "dssp_helix_mean_fraction",
            "dssp_sheet_mean_fraction",
            "dssp_turn_bend_mean_fraction",
            "dssp_coil_mean_fraction",
            "dssp_helix_mean_count",
            "dssp_sheet_mean_count",
            "dssp_turn_bend_mean_count",
            "dssp_coil_mean_count",
        ]

        for m in metrics:
            out[f"wt_{m}"] = wt.get(m, np.nan)
            out[f"mut_{m}"] = mut.get(m, np.nan)
            out[f"delta_mut_minus_wt_{m}"] = mut.get(m, np.nan) - wt.get(m, np.nan)

        rows.append(out)

    return pd.DataFrame(rows)


def main():
    if not COMPLETION_TABLE.exists():
        raise FileNotFoundError(f"Missing deep completion table: {COMPLETION_TABLE}")

    completion = pd.read_csv(COMPLETION_TABLE)
    complete = completion[completion["completion_status"] == "DEEP_TRUE_COMPLETE"].copy()

    if complete.empty:
        raise RuntimeError("No DEEP_TRUE_COMPLETE systems found.")

    print(f"Deep-complete simulations found: {len(complete)}")

    all_rows = []
    dssp_long_rows = []

    for _, r in complete.iterrows():
        sysdir = Path(r["path"])
        system = r["system"]
        temperature_K = int(str(r["temperature"]).replace("K", ""))

        protein, state, mutation, organism_type, gh_family, label = parse_system_name(system)

        print(f"Processing: {system} {temperature_K}K")

        # Core metrics from existing xvg files
        rmsd = summarize_xvg(sysdir / "rmsd_protein.xvg")
        rmsf = summarize_xvg(sysdir / "rmsf_protein.xvg")
        rg = summarize_xvg(sysdir / "gyrate_protein.xvg")
        hb = summarize_xvg(sysdir / "hbnum_protein.xvg")

        row = {
            "system": system,
            "protein": protein,
            "state": state,
            "mutation": mutation,
            "organism_type": organism_type,
            "gh_family": gh_family,
            "label": label,
            "temperature_K": temperature_K,
            "path": str(sysdir),

            "rmsd_mean_nm": rmsd["mean"],
            "rmsd_sd_nm": rmsd["sd"],
            "rmsd_late_mean_nm": rmsd["late_mean"],
            "rmsd_late_sd_nm": rmsd["late_sd"],
            "rmsd_final_nm": rmsd["final"],

            "rmsf_mean_nm": rmsf["mean"],
            "rmsf_sd_nm": rmsf["sd"],
            "rmsf_late_mean_nm": rmsf["late_mean"],
            "rmsf_late_sd_nm": rmsf["late_sd"],
            "rmsf_final_nm": rmsf["final"],

            "rg_mean_nm": rg["mean"],
            "rg_sd_nm": rg["sd"],
            "rg_late_mean_nm": rg["late_mean"],
            "rg_late_sd_nm": rg["late_sd"],
            "rg_final_nm": rg["final"],

            "hbond_mean": hb["mean"],
            "hbond_sd": hb["sd"],
            "hbond_late_mean": hb["late_mean"],
            "hbond_late_sd": hb["late_sd"],
            "hbond_final": hb["final"],
        }

        # DSSP
        ok, msg = run_dssp_for_system(sysdir)
        row["dssp_available"] = ok
        row["dssp_message"] = msg

        if ok:
            dssp_stats, dssp_df = parse_scount_xvg(sysdir / "scount.xvg")
            row.update(dssp_stats)

            if not dssp_df.empty:
                dssp_df = dssp_df.copy()
                dssp_df["system"] = system
                dssp_df["protein"] = protein
                dssp_df["state"] = state
                dssp_df["mutation"] = mutation
                dssp_df["organism_type"] = organism_type
                dssp_df["gh_family"] = gh_family
                dssp_df["temperature_K"] = temperature_K
                dssp_long_rows.append(dssp_df)

        all_rows.append(row)

    df = pd.DataFrame(all_rows)
    df = df.sort_values(["protein", "state", "temperature_K"]).reset_index(drop=True)
    df.to_csv(OUT_24, index=False)

    pair = compute_secondary_retention(df)
    pair = pair.sort_values(["protein", "temperature_K"]).reset_index(drop=True)
    pair.to_csv(OUT_PAIR, index=False)

    if dssp_long_rows:
        long = pd.concat(dssp_long_rows, ignore_index=True)
        long.to_csv(OUT_DSSP_LONG, index=False)
    else:
        pd.DataFrame().to_csv(OUT_DSSP_LONG, index=False)

    print("\nSaved:")
    print(OUT_24)
    print(OUT_PAIR)
    print(OUT_DSSP_LONG)

    print("\n24-simulation table preview:")
    show_cols = [
        "protein", "state", "mutation", "organism_type", "gh_family", "temperature_K",
        "rmsd_mean_nm", "rmsf_mean_nm", "rg_mean_nm", "hbond_mean",
        "dssp_helix_mean_fraction", "dssp_sheet_mean_fraction",
        "dssp_coil_mean_fraction",
        "dssp_available"
    ]
    show_cols = [c for c in show_cols if c in df.columns]
    print(df[show_cols].to_string(index=False))

    print("\nWT/MUT table preview:")
    show_pair_cols = [
        "protein", "mutation", "temperature_K",
        "wt_rmsd_mean_nm", "mut_rmsd_mean_nm", "delta_mut_minus_wt_rmsd_mean_nm",
        "wt_rmsf_mean_nm", "mut_rmsf_mean_nm", "delta_mut_minus_wt_rmsf_mean_nm",
        "wt_rg_mean_nm", "mut_rg_mean_nm", "delta_mut_minus_wt_rg_mean_nm",
        "wt_hbond_mean", "mut_hbond_mean", "delta_mut_minus_wt_hbond_mean",
        "wt_dssp_helix_mean_fraction", "mut_dssp_helix_mean_fraction", "delta_mut_minus_wt_dssp_helix_mean_fraction",
        "wt_dssp_sheet_mean_fraction", "mut_dssp_sheet_mean_fraction", "delta_mut_minus_wt_dssp_sheet_mean_fraction",
    ]
    show_pair_cols = [c for c in show_pair_cols if c in pair.columns]
    print(pair[show_pair_cols].to_string(index=False))


if __name__ == "__main__":
    main()
