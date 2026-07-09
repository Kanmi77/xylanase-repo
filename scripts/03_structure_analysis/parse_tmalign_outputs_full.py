#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

indir = Path("results/structures/tmalign_outputs_full")
outdir = Path("results/structures/tmalign_numeric_summary")
outdir.mkdir(parents=True, exist_ok=True)

rows = []

for f in sorted(indir.glob("*.out")):
    txt = f.read_text(errors="ignore")

    stem = f.stem
    parts = stem.split("_")

    query_accession = parts[0] if parts else None
    reference_id = "_".join(parts[1:]).replace("_Repair", "") if len(parts) > 1 else None

    # TM-align outputs often contain two TM-scores.
    # We collect both and keep the first as the query-normalised score where available.
    tm_scores = re.findall(r"TM-score=\s*([0-9.]+)", txt)
    tm_score_1 = float(tm_scores[0]) if len(tm_scores) >= 1 else None
    tm_score_2 = float(tm_scores[1]) if len(tm_scores) >= 2 else None

    # Example:
    # Aligned length=  300, RMSD=   1.23, Seq_ID=n_identical/n_aligned= 0.456
    m = re.search(
        r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+),\s*Seq_ID=n_identical/n_aligned=\s*([0-9.]+)",
        txt
    )

    aligned_length = int(m.group(1)) if m else None
    rmsd = float(m.group(2)) if m else None
    seq_id = float(m.group(3)) if m else None

    rows.append({
        "file": str(f),
        "query_accession": query_accession,
        "reference_id": reference_id,
        "tm_score_1": tm_score_1,
        "tm_score_2": tm_score_2,
        "tm_score_best_of_two": max([x for x in [tm_score_1, tm_score_2] if x is not None], default=None),
        "aligned_length": aligned_length,
        "rmsd": rmsd,
        "seq_id_aligned": seq_id,
    })

df = pd.DataFrame(rows)

full_out = outdir / "tmalign_all_pairwise_summary.csv"
df.to_csv(full_out, index=False)

# Best reference per query accession using tm_score_best_of_two
valid = df.dropna(subset=["tm_score_best_of_two"]).copy()
best = (
    valid.sort_values(
        ["query_accession", "tm_score_best_of_two", "seq_id_aligned", "aligned_length"],
        ascending=[True, False, False, False]
    )
    .groupby("query_accession", as_index=False)
    .head(1)
)

best_out = outdir / "tmalign_best_reference_per_query.csv"
best.to_csv(best_out, index=False)

# General summary
summary = {
    "total_output_files": len(df),
    "parsed_tm_score_rows": int(df["tm_score_best_of_two"].notna().sum()),
    "unique_query_accessions": int(df["query_accession"].nunique()),
    "unique_reference_ids": int(df["reference_id"].nunique()),
    "best_reference_rows": len(best),
}

numeric_cols = ["tm_score_best_of_two", "rmsd", "aligned_length", "seq_id_aligned"]
desc = df[numeric_cols].describe().T

summary_report = outdir / "TMALIGN_NUMERIC_SUMMARY_REPORT.md"
with open(summary_report, "w") as h:
    h.write("# TM-align numeric summary report\n\n")
    for k, v in summary.items():
        h.write(f"- {k}: {v}\n")

    h.write("\n## Overall numeric summary\n\n")
    h.write(desc.to_markdown())
    h.write("\n\n## Best-reference numeric summary\n\n")
    if len(best):
        h.write(best[numeric_cols].describe().T.to_markdown())
    else:
        h.write("No valid best-reference rows were parsed.\n")

    h.write("\n\n## Top 20 best model-reference matches\n\n")
    if len(best):
        cols = [
            "query_accession",
            "reference_id",
            "tm_score_best_of_two",
            "rmsd",
            "aligned_length",
            "seq_id_aligned",
        ]
        h.write(best.sort_values("tm_score_best_of_two", ascending=False)[cols].head(20).to_markdown(index=False))
    else:
        h.write("No valid best-reference rows were parsed.\n")

print("Done")
print("Full pairwise summary:", full_out)
print("Best reference summary:", best_out)
print("Report:", summary_report)
