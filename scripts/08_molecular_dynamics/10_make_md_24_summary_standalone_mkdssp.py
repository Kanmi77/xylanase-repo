#!/usr/bin/env python3

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

OUT_24 = ANALYSIS / "md_24_simulation_mean_metrics_with_standalone_dssp.csv"
OUT_PAIR = ANALYSIS / "md_wt_mut_333K_373K_table_with_standalone_dssp.csv"
OUT_DSSP_LONG = ANALYSIS / "md_standalone_dssp_framewise_long.csv"

GMX = "gmx"
MKDSSP = shutil.which("mkdssp") or shutil.which("dssp")

# Use 11 frames per simulation: 0, 1000, ..., 10000 ps.
# Increase density later if needed, but this is good for thesis-level retention.
DSSP_FRAME_TIMES_PS = list(range(0, 10001, 1000))


def run_cmd(cmd, cwd=None, input_text=None, timeout=1200):
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

def extract_frame_pdb(sysdir, frame_time_ps, out_pdb):
    """
    Extract one protein-only PDB frame at requested time.

    Important:
    out_pdb is inside standalone_dssp/, so we must give gmx trjconv
    the full output path, not only out_pdb.name.
    """
    if out_pdb.exists() and out_pdb.stat().st_size > 0:
        return True, "frame exists"

    cmd = (
        f'printf "Protein\\n" | {GMX} trjconv '
        f'-s md.tpr '
        f'-f md.xtc '
        f'-o "{out_pdb}" '
        f'-dump {frame_time_ps}'
    )

    code, out, err = run_cmd(cmd, cwd=sysdir, timeout=600)
    combined = out + "\n" + err

    if out_pdb.exists() and out_pdb.stat().st_size > 0:
        return True, "frame extracted"

    return False, combined[-1500:]

def make_dssp_compatible_pdb(pdb):
    """
    Create a clean DSSP-compatible PDB.

    Fixes:
    - Adds HEADER and TITLE
    - Keeps only heavy protein ATOM records
    - Forces chain ID to A
    - Removes GROMACS REMARK/MODEL/CRYST1 records
    - Omits TER to avoid mkdssp warning
    """
    if not pdb.exists() or pdb.stat().st_size == 0:
        return False, "PDB frame missing or empty"

    fixed = pdb.with_name(pdb.stem + "_dssp_clean.pdb")

    if fixed.exists() and fixed.stat().st_size > 0:
        return True, str(fixed)

    out_lines = [
        "HEADER    GROMACS MD FRAME CLEANED FOR DSSP\n",
        f"TITLE     {pdb.name}\n",
    ]

    atom_count = 0
    serial = 1

    with open(pdb, errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) >= 78 else ""

            # Remove hydrogens
            if atom_name.startswith("H") or element == "H":
                continue

            line = line.rstrip("\n")
            if len(line) < 80:
                line = line.ljust(80)

            # Force ATOM record, new serial number, and chain ID A
            new_line = (
                "ATOM  "
                + f"{serial:5d}"
                + line[11:21]
                + "A"
                + line[22:80]
            )

            out_lines.append(new_line + "\n")
            atom_count += 1
            serial += 1

    out_lines.append("END\n")

    if atom_count == 0:
        return False, "No ATOM records found after cleaning"

    fixed.write_text("".join(out_lines))

    return True, str(fixed)

def run_mkdssp_on_pdb(pdb, dssp_out):
    """
    Run mkdssp on a cleaned DSSP-compatible PDB.
    """
    if dssp_out.exists() and dssp_out.stat().st_size > 0:
        return True, "DSSP exists"

    ok, fixed_msg = make_dssp_compatible_pdb(pdb)

    if not ok:
        return False, fixed_msg

    fixed_pdb = Path(fixed_msg)

    attempts = [
        f'{MKDSSP} "{fixed_pdb}" "{dssp_out}"',
        f'{MKDSSP} "{fixed_pdb}" > "{dssp_out}"',
    ]

    messages = []

    for cmd in attempts:
        code, out, err = run_cmd(cmd, cwd=pdb.parent, timeout=600)

        messages.append(
            f"CMD: {cmd}\nCODE: {code}\nSTDOUT:\n{out[-1000:]}\nSTDERR:\n{err[-1000:]}"
        )

        if dssp_out.exists() and dssp_out.stat().st_size > 0:
            text = dssp_out.read_text(errors="ignore")
            if "RESIDUE" in text and "STRUCTURE" in text:
                return True, f"DSSP completed with: {cmd}"

    return False, "\n\n".join(messages)[-3000:]

def parse_dssp_file(dssp_file):
    """
    Parse classic DSSP output robustly.

    Supports DSSP output where the residue table begins after a line containing:
    '#  RESIDUE AA STRUCTURE'

    Secondary structure:
    H/G/I = helix
    E/B   = sheet/bridge
    T/S   = turn/bend
    blank/other = coil
    """
    if not dssp_file.exists() or dssp_file.stat().st_size == 0:
        return None

    text = dssp_file.read_text(errors="ignore")

    started = False

    counts = {
        "helix": 0,
        "sheet": 0,
        "turn_bend": 0,
        "coil": 0,
        "total": 0,
    }

    for line in text.splitlines():
        if "RESIDUE" in line and "AA" in line and "STRUCTURE" in line:
            started = True
            continue

        if not started:
            continue

        if len(line) < 20:
            continue

        # Classic DSSP residue rows usually contain:
        # columns 6-10 = DSSP residue index
        # column 14 = amino acid
        # column 17 = secondary structure code
        residue_index = line[0:5].strip()
        aa = line[13:14].strip() if len(line) >= 14 else ""
        ss = line[16:17].strip() if len(line) >= 17 else ""

        # Skip non-residue or chain-break rows
        if not residue_index or not residue_index.replace("-", "").isdigit():
            continue

        if aa in {"!", "*"}:
            continue

        counts["total"] += 1

        if ss in {"H", "G", "I"}:
            counts["helix"] += 1
        elif ss in {"E", "B"}:
            counts["sheet"] += 1
        elif ss in {"T", "S"}:
            counts["turn_bend"] += 1
        else:
            counts["coil"] += 1

    if counts["total"] == 0:
        return None

    return {
        "dssp_helix_count": counts["helix"],
        "dssp_sheet_count": counts["sheet"],
        "dssp_turn_bend_count": counts["turn_bend"],
        "dssp_coil_count": counts["coil"],
        "dssp_total_count": counts["total"],
        "dssp_helix_fraction": counts["helix"] / counts["total"],
        "dssp_sheet_fraction": counts["sheet"] / counts["total"],
        "dssp_turn_bend_fraction": counts["turn_bend"] / counts["total"],
        "dssp_coil_fraction": counts["coil"] / counts["total"],
    }
def summarize_dssp_for_system(sysdir):
    if MKDSSP is None:
        return False, "mkdssp/dssp executable not found", {}, pd.DataFrame()

    dssp_dir = sysdir / "standalone_dssp"
    dssp_dir.mkdir(exist_ok=True)

    frame_rows = []
    messages = []

    for t in DSSP_FRAME_TIMES_PS:
        pdb = dssp_dir / f"frame_{t:05d}ps.pdb"
        dssp_out = dssp_dir / f"frame_{t:05d}ps.dssp"

        ok_frame, msg_frame = extract_frame_pdb(sysdir, t, pdb)
        if not ok_frame:
            messages.append(f"{t} ps frame failed: {msg_frame}")
            continue

        ok_dssp, msg_dssp = run_mkdssp_on_pdb(pdb, dssp_out)
        if not ok_dssp:
            messages.append(f"{t} ps DSSP failed: {msg_dssp}")
            continue

        parsed = parse_dssp_file(dssp_out)
        if parsed is None:
            messages.append(f"{t} ps DSSP parsed zero residues")
            continue

        parsed["time_ps"] = t
        frame_rows.append(parsed)

    if not frame_rows:
        return False, "; ".join(messages)[-2000:], {}, pd.DataFrame()

    frame_df = pd.DataFrame(frame_rows).sort_values("time_ps")

    stats = {}

    for col in [
        "dssp_helix_count",
        "dssp_sheet_count",
        "dssp_turn_bend_count",
        "dssp_coil_count",
        "dssp_total_count",
        "dssp_helix_fraction",
        "dssp_sheet_fraction",
        "dssp_turn_bend_fraction",
        "dssp_coil_fraction",
    ]:
        stats[f"{col}_mean"] = frame_df[col].mean()
        stats[f"{col}_sd"] = frame_df[col].std(ddof=1)
        stats[f"{col}_final"] = frame_df[col].iloc[-1]

    # Secondary-structure retention proxy:
    # total ordered fraction = helix + sheet
    frame_df["dssp_ordered_fraction"] = (
        frame_df["dssp_helix_fraction"] + frame_df["dssp_sheet_fraction"]
    )

    stats["dssp_ordered_fraction_mean"] = frame_df["dssp_ordered_fraction"].mean()
    stats["dssp_ordered_fraction_sd"] = frame_df["dssp_ordered_fraction"].std(ddof=1)
    stats["dssp_ordered_fraction_final"] = frame_df["dssp_ordered_fraction"].iloc[-1]

    return True, f"DSSP completed on {len(frame_df)} frames", stats, frame_df


def compute_wt_mut_table(df):
    rows = []

    for (protein, temperature), sub in df.groupby(["protein", "temperature_K"]):
        wt = sub[sub["state"] == "WT"]
        mut = sub[sub["state"] == "MUT"]

        if wt.empty or mut.empty:
            continue

        wt = wt.iloc[0]
        mut = mut.iloc[0]

        out = {
            "protein": protein,
            "mutation": mut.get("mutation", ""),
            "organism_type": mut.get("organism_type", wt.get("organism_type", "")),
            "gh_family": mut.get("gh_family", wt.get("gh_family", "")),
            "temperature_K": temperature,
        }

        metrics = [
            "rmsd_mean_nm",
            "rmsf_mean_nm",
            "rg_mean_nm",
            "hbond_mean",

            "rmsd_late_mean_nm",
            "rmsf_late_mean_nm",
            "rg_late_mean_nm",
            "hbond_late_mean",

            "dssp_helix_fraction_mean",
            "dssp_sheet_fraction_mean",
            "dssp_turn_bend_fraction_mean",
            "dssp_coil_fraction_mean",
            "dssp_ordered_fraction_mean",
            "dssp_helix_count_mean",
            "dssp_sheet_count_mean",
            "dssp_total_count_mean",
        ]

        for m in metrics:
            out[f"wt_{m}"] = wt.get(m, np.nan)
            out[f"mut_{m}"] = mut.get(m, np.nan)
            out[f"delta_mut_minus_wt_{m}"] = mut.get(m, np.nan) - wt.get(m, np.nan)

        rows.append(out)

    return pd.DataFrame(rows).sort_values(["protein", "temperature_K"])


def main():
    if MKDSSP is None:
        print("WARNING: mkdssp/dssp not found. DSSP will be unavailable.")
    else:
        print("Using DSSP executable:", MKDSSP)

    if not COMPLETION_TABLE.exists():
        raise FileNotFoundError(f"Missing deep completion table: {COMPLETION_TABLE}")

    completion = pd.read_csv(COMPLETION_TABLE)
    complete = completion[completion["completion_status"] == "DEEP_TRUE_COMPLETE"].copy()

    print("Deep-complete simulations:", len(complete))

    rows = []
    long_rows = []

    for _, r in complete.iterrows():
        sysdir = Path(r["path"])
        system = r["system"]
        temperature_K = int(str(r["temperature"]).replace("K", ""))

        protein, state, mutation, organism_type, gh_family, label = parse_system_name(system)

        print(f"Processing {system} {temperature_K}K")

        rmsd = summarize_xvg(sysdir / "rmsd_protein.xvg")
        rmsf = summarize_xvg(sysdir / "rmsf_protein.xvg")
        rg = summarize_xvg(sysdir / "gyrate_protein.xvg")
        hb = summarize_xvg(sysdir / "hbnum_protein.xvg")

        ok_dssp, msg_dssp, dssp_stats, dssp_frame_df = summarize_dssp_for_system(sysdir)

        row = {
            "system": system,
            "protein": protein,
            "state": state,
            "mutation": mutation,
            "organism_type": organism_type,
            "gh_family": gh_family,
            "temperature_K": temperature_K,
            "path": str(sysdir),

            "rmsd_mean_nm": rmsd["mean"],
            "rmsd_sd_nm": rmsd["sd"],
            "rmsd_late_mean_nm": rmsd["late_mean"],
            "rmsd_final_nm": rmsd["final"],

            "rmsf_mean_nm": rmsf["mean"],
            "rmsf_sd_nm": rmsf["sd"],
            "rmsf_late_mean_nm": rmsf["late_mean"],
            "rmsf_final_nm": rmsf["final"],

            "rg_mean_nm": rg["mean"],
            "rg_sd_nm": rg["sd"],
            "rg_late_mean_nm": rg["late_mean"],
            "rg_final_nm": rg["final"],

            "hbond_mean": hb["mean"],
            "hbond_sd": hb["sd"],
            "hbond_late_mean": hb["late_mean"],
            "hbond_final": hb["final"],

            "dssp_available": ok_dssp,
            "dssp_message": msg_dssp,
        }

        row.update(dssp_stats)
        rows.append(row)

        if ok_dssp and not dssp_frame_df.empty:
            tmp = dssp_frame_df.copy()
            tmp["system"] = system
            tmp["protein"] = protein
            tmp["state"] = state
            tmp["mutation"] = mutation
            tmp["organism_type"] = organism_type
            tmp["gh_family"] = gh_family
            tmp["temperature_K"] = temperature_K
            long_rows.append(tmp)

    df = pd.DataFrame(rows).sort_values(["protein", "state", "temperature_K"])
    df.to_csv(OUT_24, index=False)

    pair = compute_wt_mut_table(df)
    pair.to_csv(OUT_PAIR, index=False)

    if long_rows:
        long = pd.concat(long_rows, ignore_index=True)
        long.to_csv(OUT_DSSP_LONG, index=False)
    else:
        pd.DataFrame().to_csv(OUT_DSSP_LONG, index=False)

    print("\nSaved:")
    print(OUT_24)
    print(OUT_PAIR)
    print(OUT_DSSP_LONG)

    print("\nDSSP availability:")
    print(df["dssp_available"].value_counts(dropna=False).to_string())

    print("\n24 simulation summary:")
    show = [
        "protein", "state", "mutation", "organism_type", "gh_family", "temperature_K",
        "rmsd_mean_nm", "rmsf_mean_nm", "rg_mean_nm", "hbond_mean",
        "dssp_helix_fraction_mean",
        "dssp_sheet_fraction_mean",
        "dssp_ordered_fraction_mean",
        "dssp_coil_fraction_mean",
        "dssp_available",
    ]
    show = [c for c in show if c in df.columns]
    print(df[show].to_string(index=False))

    print("\nWT/MUT 333K/373K table:")
    show_pair = [
        "protein", "mutation", "organism_type", "gh_family", "temperature_K",
        "wt_rmsd_mean_nm", "mut_rmsd_mean_nm", "delta_mut_minus_wt_rmsd_mean_nm",
        "wt_rmsf_mean_nm", "mut_rmsf_mean_nm", "delta_mut_minus_wt_rmsf_mean_nm",
        "wt_rg_mean_nm", "mut_rg_mean_nm", "delta_mut_minus_wt_rg_mean_nm",
        "wt_hbond_mean", "mut_hbond_mean", "delta_mut_minus_wt_hbond_mean",
        "wt_dssp_helix_fraction_mean", "mut_dssp_helix_fraction_mean",
        "delta_mut_minus_wt_dssp_helix_fraction_mean",
        "wt_dssp_sheet_fraction_mean", "mut_dssp_sheet_fraction_mean",
        "delta_mut_minus_wt_dssp_sheet_fraction_mean",
        "wt_dssp_ordered_fraction_mean", "mut_dssp_ordered_fraction_mean",
        "delta_mut_minus_wt_dssp_ordered_fraction_mean",
    ]
    show_pair = [c for c in show_pair if c in pair.columns]
    print(pair[show_pair].to_string(index=False))


if __name__ == "__main__":
    main()
