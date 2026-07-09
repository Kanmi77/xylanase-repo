#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

OUTDIR = Path("results/experimental_ml_correction/combined_60c_sequence_structure_classifier")
OUTDIR.mkdir(parents=True, exist_ok=True)

PRIMARY_MAPPING = OUTDIR / "sequence_structure_mapping_initial.csv"
STRUCTURE_AVAILABILITY = OUTDIR / "missing_labelled_structure_availability.csv"

final = pd.read_csv(PRIMARY_MAPPING)
avail = pd.read_csv(STRUCTURE_AVAILABILITY)

# Keep original final matches
final = final.copy()

# Prepare extra matches from rows that were already in the structural table but missed
recoverable = avail[
    avail["status"].eq("already_in_structural_table_but_mapping_missed")
].copy()

recoverable["already_in_structural_table"] = recoverable["already_in_structural_table"].fillna("")
recoverable["chosen_recovered_accession"] = recoverable["already_in_structural_table"].apply(
    lambda x: str(x).split(";")[0] if str(x).strip() else ""
)

# Merge recovery information back onto sequence-structure mapping table by sequence_hash if available;
# otherwise use accession_or_name + thermal_label + sequence_length as fallback.
if "sequence_hash" in recoverable.columns and "sequence_hash" in final.columns:
    final = final.merge(
        recoverable[["sequence_hash", "chosen_recovered_accession"]],
        on="sequence_hash",
        how="left"
    )
else:
    final = final.merge(
        recoverable[["accession_or_name", "thermal_label", "sequence_length", "chosen_recovered_accession"]],
        on=["accession_or_name", "thermal_label", "sequence_length"],
        how="left"
    )

final["chosen_recovered_accession"] = final["chosen_recovered_accession"].fillna("")

# Fill missing final matches with recovered accessions
mask = (~final["has_final_structural_match"].astype(bool)) & (final["chosen_recovered_accession"] != "")

final.loc[mask, "chosen_structural_accession"] = final.loc[mask, "chosen_recovered_accession"]
final.loc[mask, "match_method"] = "availability_structural_table"
final.loc[mask, "has_final_structural_match"] = True

final = final.drop(columns=["chosen_recovered_accession"], errors="ignore")

matched = final[final["has_final_structural_match"].astype(bool)].copy()
missing = final[~final["has_final_structural_match"].astype(bool)].copy()

print("Original final matched rows:", pd.read_csv(PRIMARY_MAPPING)["has_final_structural_match"].sum())
print("Recovered from structure-availability table:", int(mask.sum()))
print("New final matched rows:", len(matched))
print("New unique structural accessions:", matched["chosen_structural_accession"].nunique())

print("\nMatch method counts:")
print(matched["match_method"].value_counts())

print("\nMatched class counts:")
print(matched["thermal_label"].value_counts())

print("\nStill missing class counts:")
print(missing["thermal_label"].value_counts())

final.to_csv(OUTDIR / "sequence_structure_mapping.csv", index=False)
matched.to_csv(OUTDIR / "sequence_structure_matched_rows.csv", index=False)
missing.to_csv(OUTDIR / "sequence_structure_unmatched_rows.csv", index=False)

print("\nSaved:")
print(OUTDIR / "sequence_structure_mapping.csv")
print(OUTDIR / "sequence_structure_matched_rows.csv")
print(OUTDIR / "sequence_structure_unmatched_rows.csv")
