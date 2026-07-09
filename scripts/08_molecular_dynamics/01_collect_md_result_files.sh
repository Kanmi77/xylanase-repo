#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"
cd "$BASE"

OUT="results/molecular_dynamics_results_for_git"
PDB_OUT="$OUT/phase1_pdb_referenced_baseline"
MOD_OUT="$OUT/phase2_modeller_wt_mutant_validation"
MAN_OUT="$OUT/manifests"

rm -rf "$OUT"
mkdir -p "$PDB_OUT" "$MOD_OUT" "$MAN_OUT"

echo "Extracting MD result data for Git..."
echo "Output directory: $OUT"

copy_file_if_exists () {
    local src="$1"
    local destroot="$2"
    if [ -f "$src" ]; then
        mkdir -p "$destroot/$(dirname "$src")"
        cp -p "$src" "$destroot/$src"
        echo "COPIED: $src"
    fi
}

copy_find_results () {
    local search_root="$1"
    local destroot="$2"
    shift 2

    if [ -d "$search_root" ]; then
        find "$search_root" -type f "$@" -print0 | while IFS= read -r -d '' f; do
            mkdir -p "$destroot/$(dirname "$f")"
            cp -p "$f" "$destroot/$f"
        done
    fi
}

###############################################################################
# PHASE 1: PDB-REFERENCED BASELINE MD
###############################################################################

echo
echo "=== PHASE 1: PDB-referenced baseline MD ==="

PDB_IDS=("1O8S" "1T6G" "1VBU" "2C79" "3NIY" "7K4X" "1VBR" "7K4P" "8RD5")

# Main PDB/reference MD summary datasets
PDB_SUMMARY_FILES=(
    "results/reports/previous_pdb_md_presence_summary.csv"
    "results/reports/previous_pdb_md_summary_by_pdb.csv"
    "results/reports/md_previous_pdb_vs_current_long_metrics.csv"
    "results/reports/md_previous_pdb_vs_current_long_metrics.csv"
    "results/reports/md_previous_pdb_vs_current_summary_wide.csv"
    "results/reports/md_previous_pdb_vs_current_summary_wide.csv"
    "results/reports/md_previous_pdb_vs_current_group_comparison.csv"
    "results/reports/md_previous_pdb_vs_current_group_comparison.csv"
    "results/md_stage1_metrics_summary.csv"
    "md/summary/stage1_md_metrics.csv"
    "md_10ns/plots/md_summary_metrics.csv"
)

for f in "${PDB_SUMMARY_FILES[@]}"; do
    copy_file_if_exists "$f" "$PDB_OUT"
done

# PDB baseline plots
if [ -d "md_10ns/plots" ]; then
    for pdb in "${PDB_IDS[@]}"; do
        find md_10ns/plots -type f \
            \( -name "${pdb}*rmsd*.png" \
            -o -name "${pdb}*rmsf*.png" \
            -o -name "${pdb}*rg*.png" \
            -o -name "${pdb}*gyrate*.png" \
            -o -name "${pdb}*hbond*.png" \
            -o -name "${pdb}*sasa*.png" \) \
            -print0 | while IFS= read -r -d '' f; do
                mkdir -p "$PDB_OUT/$(dirname "$f")"
                cp -p "$f" "$PDB_OUT/$f"
            done
    done
fi

# PDB per-run analysis files from possible PDB MD folders
PDB_SEARCH_ROOTS=(
    "md_10ns/systems"
    "md/systems"
    "md/_archive"
)

for root in "${PDB_SEARCH_ROOTS[@]}"; do
    [ -d "$root" ] || continue

    for pdb in "${PDB_IDS[@]}"; do
        find "$root" -type f -path "*${pdb}*" \
            \( -name "*rmsd*.xvg" \
            -o -name "*rmsf*.xvg" \
            -o -name "*rg*.xvg" \
            -o -name "*gyrate*.xvg" \
            -o -name "*hbond*.xvg" \
            -o -name "*hbnum*.xvg" \
            -o -name "*sasa*.xvg" \
            -o -name "ss_counts.xvg" \
            -o -name "ss.xpm" \
            -o -name "*.png" \
            -o -name "*.csv" \
            -o -name "*.log" \) \
            -print0 | while IFS= read -r -d '' f; do
                mkdir -p "$PDB_OUT/$(dirname "$f")"
                cp -p "$f" "$PDB_OUT/$f"
            done
    done
done

# Manifest of raw/heavy PDB MD files, not copied to Git folder
{
    echo "Raw/heavy PDB MD files present in VM but not copied into Git result folder:"
    for root in "${PDB_SEARCH_ROOTS[@]}"; do
        [ -d "$root" ] || continue
        find "$root" -type f \
            \( -name "*.xtc" \
            -o -name "*.trr" \
            -o -name "*.tpr" \
            -o -name "*.edr" \
            -o -name "*.cpt" \
            -o -name "*.gro" \) \
            -printf "%p\t%s bytes\n" 2>/dev/null || true
    done
} > "$MAN_OUT/phase1_pdb_raw_heavy_files_manifest.tsv"

###############################################################################
# PHASE 2: MODELLER / HOMOLOGY WT-MUTANT MD
###############################################################################

echo
echo "=== PHASE 2: MODELLER / homology WT-mutant MD ==="

# Main MODELLER WT-mutant MD summary datasets
MOD_SUMMARY_FILES=(
    "md_tier2_wt_mutant_compact/compact_md_completion_summary.csv"
    "md_tier2_wt_mutant_compact/compact_md_deep_completion_summary.csv"
    "md_tier2_wt_mutant_compact/analysis/compact_md_metrics_summary.csv"
    "md_tier2_wt_mutant_compact/analysis/compact_md_performance_summary.csv"
    "md_tier2_wt_mutant_compact/analysis/compact_md_rmsf_summary.csv"
    "md_tier2_wt_mutant_compact/analysis/compact_md_temperature_comparison.csv"
    "md_tier2_wt_mutant_compact/analysis/compact_md_wt_mutant_comparison.csv"
    "md_tier2_wt_mutant_compact/analysis/md_24_simulation_mean_metrics_with_dssp.csv"
    "md_tier2_wt_mutant_compact/analysis/md_24_simulation_mean_metrics_with_standalone_dssp.csv"
    "results/reports/md_current_mutant_minus_wt_differences.csv"
)

for f in "${MOD_SUMMARY_FILES[@]}"; do
    copy_file_if_exists "$f" "$MOD_OUT"
done

# Copy every CSV/TSV/PNG generated inside the compact MD analysis folder
copy_find_results "md_tier2_wt_mutant_compact/analysis" "$MOD_OUT" \
    \( -name "*.csv" -o -name "*.tsv" -o -name "*.png" -o -name "*.txt" \)

# Copy per-run analysis outputs for WT-mutant systems
copy_find_results "md_tier2_wt_mutant_compact/systems" "$MOD_OUT" \
    \( -name "*rmsd*.xvg" \
    -o -name "*rmsf*.xvg" \
    -o -name "*rg*.xvg" \
    -o -name "*gyrate*.xvg" \
    -o -name "*hbond*.xvg" \
    -o -name "*hbnum*.xvg" \
    -o -name "*sasa*.xvg" \
    -o -name "ss_counts.xvg" \
    -o -name "ss.xpm" \
    -o -name "*.dssp" \
    -o -name "*dssp*.csv" \
    -o -name "*.png" \
    -o -name "*.log" \)

# Copy FoldX-vs-MD correlation result figures/data if present
copy_find_results "results/correlation" "$MOD_OUT" \
    \( -name "13_foldx_ddg_vs_delta*.png" \
    -o -name "13_foldx_ddg_vs_delta*.csv" \
    -o -name "*md*mut_minus_wt*.csv" \
    -o -name "*md*mutant_minus_wt*.csv" \)

copy_find_results "results/figures/integration" "$MOD_OUT" \
    \( -name "*md*.png" -o -name "*MD*.png" -o -name "*.csv" \)

# Manifest of raw/heavy MODELLER MD files, not copied to Git folder
{
    echo "Raw/heavy MODELLER WT-mutant MD files present in VM but not copied into Git result folder:"
    find md_tier2_wt_mutant_compact -type f \
        \( -name "*.xtc" \
        -o -name "*.trr" \
        -o -name "*.tpr" \
        -o -name "*.edr" \
        -o -name "*.cpt" \
        -o -name "*.gro" \) \
        -printf "%p\t%s bytes\n" 2>/dev/null || true
} > "$MAN_OUT/phase2_modeller_raw_heavy_files_manifest.tsv"

###############################################################################
# GENERAL MANIFESTS AND README
###############################################################################

echo
echo "=== Writing manifests ==="

find "$PDB_OUT" -type f | sort > "$MAN_OUT/phase1_pdb_files_copied.txt"
find "$MOD_OUT" -type f | sort > "$MAN_OUT/phase2_modeller_files_copied.txt"
find "$OUT" -type f | sort > "$MAN_OUT/all_md_files_copied.txt"

{
    echo "# Molecular dynamics result extraction for Git"
    echo
    echo "This folder separates molecular dynamics result data according to the thesis methodology."
    echo
    echo "## Phase 1: PDB-referenced baseline MD"
    echo
    echo "Folder: phase1_pdb_referenced_baseline"
    echo
    echo "Includes result summaries, per-run analysis files, plots, and logs for experimentally resolved PDB reference xylanases."
    echo
    echo "Reference systems:"
    echo "1O8S, 1T6G, 1VBU, 2C79, 3NIY, 7K4X, 1VBR, 7K4P, and 8RD5."
    echo
    echo "Expected temperatures:"
    echo "333 K, 353 K, and 373 K."
    echo
    echo "## Phase 2: MODELLER / homology WT-mutant validation"
    echo
    echo "Folder: phase2_modeller_wt_mutant_validation"
    echo
    echo "Includes result summaries, per-run analysis files, DSSP outputs, plots, and logs for WT-mutant comparative MD systems."
    echo
    echo "Selected systems:"
    echo "D5UGW9, Q60041, B0ZSE5, Q9HGX1, Q9HEN6, and E3WF08."
    echo
    echo "Expected temperatures:"
    echo "333 K and 373 K."
    echo
    echo "## Excluded from Git result folder"
    echo
    echo "Large raw trajectory and binary simulation files were not copied into this folder."
    echo "Their paths and sizes are listed in:"
    echo
    echo "- manifests/phase1_pdb_raw_heavy_files_manifest.tsv"
    echo "- manifests/phase2_modeller_raw_heavy_files_manifest.tsv"
} > "$OUT/README.md"

echo
echo "=== Extraction summary ==="
echo "PDB copied files:"
wc -l "$MAN_OUT/phase1_pdb_files_copied.txt"
echo "MODELLER copied files:"
wc -l "$MAN_OUT/phase2_modeller_files_copied.txt"
echo "All copied files:"
wc -l "$MAN_OUT/all_md_files_copied.txt"
echo
du -sh "$OUT"

echo
echo "Done. Extracted MD results are in:"
echo "$OUT"
