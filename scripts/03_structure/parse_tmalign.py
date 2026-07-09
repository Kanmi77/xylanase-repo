#!/usr/bin/env python3
# Purpose: Parse best TM-align matches.

from pathlib import Path
import pandas as pd

BASE = Path.home() / "xylanase-thesis"

IN = BASE / "results/structures/tmalign_results_GH10_GH11_only.csv"
OUT = BASE / "results/structures/tmalign_best_reference_per_model.csv"
SUMMARY = BASE / "results/reports/tmalign_best_reference_per_model_summary.txt"

df = pd.read_csv(IN)

# Standardize likely column names
cols = list(df.columns)

# Identify model column
model_col = None
for c in ["model_uniprot", "uniprot_accession", "model_accession", "model_id", "query_id"]:
    if c in cols:
        model_col = c
        break

if model_col is None:
    raise ValueError(f"Could not identify model accession column. Columns are: {cols}")

# Make best TM-score column if not already present
if "tm_score_best" not in df.columns:
    tm_cols = [c for c in df.columns if "tm_score" in c.lower()]
    if len(tm_cols) >= 2:
        df["tm_score_best"] = df[tm_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
    elif len(tm_cols) == 1:
        df["tm_score_best"] = pd.to_numeric(df[tm_cols[0]], errors="coerce")
    else:
        raise ValueError(f"No TM-score column found. Columns are: {cols}")

# Numeric cleanup
df["tm_score_best"] = pd.to_numeric(df["tm_score_best"], errors="coerce")

if "rmsd" in df.columns:
    df["rmsd"] = pd.to_numeric(df["rmsd"], errors="coerce")

# Sort: highest TM-score first, then lowest RMSD if present
sort_cols = [model_col, "tm_score_best"]
ascending = [True, False]

if "rmsd" in df.columns:
    sort_cols.append("rmsd")
    ascending.append(True)

best = (
    df.sort_values(sort_cols, ascending=ascending)
      .groupby(model_col, as_index=False)
      .head(1)
      .reset_index(drop=True)
)

OUT.parent.mkdir(parents=True, exist_ok=True)
SUMMARY.parent.mkdir(parents=True, exist_ok=True)

best.to_csv(OUT, index=False)

summary_lines = []
summary_lines.append("TM-align Best Reference Per Model Summary")
summary_lines.append("=" * 50)
summary_lines.append(f"Input: {IN}")
summary_lines.append(f"Output: {OUT}")
summary_lines.append(f"Input rows: {len(df)}")
summary_lines.append(f"Unique models: {df[model_col].nunique()}")
summary_lines.append(f"Best-reference rows: {len(best)}")
summary_lines.append(f"Model column used: {model_col}")
summary_lines.append("")
summary_lines.append("Columns:")
summary_lines.extend(best.columns.tolist())
summary_lines.append("")
summary_lines.append("TM-score summary:")
summary_lines.append(str(best["tm_score_best"].describe()))

if "rmsd" in best.columns:
    summary_lines.append("")
    summary_lines.append("RMSD summary:")
    summary_lines.append(str(best["rmsd"].describe()))

SUMMARY.write_text("\n".join(summary_lines))

print(f"Saved: {OUT}")
print(f"Saved: {SUMMARY}")
print(f"Input rows: {len(df)}")
print(f"Best-reference rows: {len(best)}")
print(f"Model column used: {model_col}")
print(best.head())
