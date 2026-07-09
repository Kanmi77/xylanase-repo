#!/usr/bin/env python3
# Purpose: Build the final master table.

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "xylanase-thesis"

MASTER = BASE / "data/curated/xylanase_master_deduplicated.csv"
REFSEQ = BASE / "data/curated/refseq_inventory.csv"
STRUCT_MANIFEST = BASE / "data/curated/combined_structure_manifest.csv"
MODELLER = BASE / "data/curated/modeller_model_manifest.csv"
FOLDX = BASE / "results/foldx/foldx_normalized.csv"
PROTPARAM = BASE / "results/sequence_features/protparam_features.csv"
SIGNAL = BASE / "results/sequence_features/signal_peptide_proxy.csv"
ML = BASE / "results/ml/structural_stability_ml_predictions.csv"

OUT = BASE / "data/curated/xylanase_thesis_master_final.csv"
SUMMARY = BASE / "results/reports/xylanase_thesis_master_summary.txt"

def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path, dtype=str).fillna("")

def nonempty(s):
    return s.astype(str).str.strip().ne("")

def as_numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def yesno(series):
    return np.where(series, True, False)

# ---------------------------------------------------------------------
# 1. Load protein-level master
# ---------------------------------------------------------------------
master = read_csv(MASTER)

# Ensure one row per UniProt accession
master = master.drop_duplicates(subset=["uniprot_accession"]).copy()

# Sequence length from sequence column
master["sequence_length"] = master["sequence"].astype(str).str.len()

# BRENDA flags
temp_cols = [
    "brenda_temperature_optimum",
    "brenda_temperature_range",
    "brenda_temperature_stability",
]
ph_cols = [
    "brenda_ph_optimum",
    "brenda_ph_range",
]
all_brenda_cols = temp_cols + ph_cols

for c in all_brenda_cols:
    if c not in master.columns:
        master[c] = ""

master["has_brenda_temp"] = yesno(master[temp_cols].apply(lambda r: any(str(x).strip() for x in r), axis=1))
master["has_brenda_ph"] = yesno(master[ph_cols].apply(lambda r: any(str(x).strip() for x in r), axis=1))
master["has_any_brenda"] = yesno(master[all_brenda_cols].apply(lambda r: any(str(x).strip() for x in r), axis=1))

# ---------------------------------------------------------------------
# 2. RefSeq support
# ---------------------------------------------------------------------
refseq = read_csv(REFSEQ)
refseq_summary = (
    refseq.groupby("uniprot_accession", as_index=False)
    .agg(
        refseq_acc_all=("refseq_acc", lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))),
        refseq_count=("refseq_acc", lambda x: len(set(str(v).strip() for v in x if str(v).strip())))
    )
)

df = master.merge(refseq_summary, on="uniprot_accession", how="left")
df["refseq_acc_all"] = df["refseq_acc_all"].fillna("")
df["refseq_count"] = pd.to_numeric(df["refseq_count"], errors="coerce").fillna(0).astype(int)
df["has_refseq_link"] = df["refseq_count"] > 0

# ---------------------------------------------------------------------
# 3. Structure manifest summary
# ---------------------------------------------------------------------
struct = read_csv(STRUCT_MANIFEST)

# Normalize available columns
if "structure_exists" in struct.columns:
    struct["structure_exists_bool"] = struct["structure_exists"].astype(str).str.lower().isin(["true", "1", "yes"])
else:
    struct["structure_exists_bool"] = True

struct_summary = (
    struct.groupby("uniprot_accession", as_index=False)
    .agg(
        structure_count=("structure_id", lambda x: len(set(str(v).strip() for v in x if str(v).strip()))),
        structure_sources=("structure_source", lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))),
        structure_ids=("structure_id", lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))),
        representative_structure_id=("structure_id", lambda x: next((str(v).strip() for v in x if str(v).strip()), "")),
        representative_structure_path=("structure_path", lambda x: next((str(v).strip() for v in x if str(v).strip()), "")),
        any_structure_exists=("structure_exists_bool", "any"),
    )
)

struct_summary["has_experimental_structure"] = struct_summary["structure_sources"].str.contains("pdb|experimental", case=False, regex=True)
struct_summary["has_modelled_structure"] = struct_summary["structure_sources"].str.contains("modeller|model", case=False, regex=True)
struct_summary["has_any_structure"] = struct_summary["any_structure_exists"] & (struct_summary["structure_count"] > 0)

df = df.merge(struct_summary, on="uniprot_accession", how="left")

for c in ["structure_sources", "structure_ids", "representative_structure_id", "representative_structure_path"]:
    df[c] = df[c].fillna("")

for c in ["structure_count"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

for c in ["has_experimental_structure", "has_modelled_structure", "has_any_structure"]:
    df[c] = df[c].fillna(False).astype(bool)

# ---------------------------------------------------------------------
# 4. MODELLER summary
# ---------------------------------------------------------------------
modeller = read_csv(MODELLER)

modeller_summary = (
    modeller.groupby("uniprot_accession", as_index=False)
    .agg(
        modeller_model_count=("model_path", lambda x: len(set(str(v).strip() for v in x if str(v).strip()))),
        modeller_template_codes=("template_code", lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))),
        modeller_best_template_identity=("template_identity", lambda x: pd.to_numeric(x, errors="coerce").max()),
        modeller_best_template_coverage=("template_coverage", lambda x: pd.to_numeric(x, errors="coerce").max()),
        modeller_run_status=("run_status", lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))),
    )
)

df = df.merge(modeller_summary, on="uniprot_accession", how="left")
df["modeller_model_count"] = pd.to_numeric(df["modeller_model_count"], errors="coerce").fillna(0).astype(int)
for c in ["modeller_template_codes", "modeller_run_status"]:
    df[c] = df[c].fillna("")

# ---------------------------------------------------------------------
# 5. FoldX and structural features
# ---------------------------------------------------------------------
foldx = read_csv(FOLDX)

foldx = as_numeric(foldx, [
    "foldx_wt_total_energy",
    "foldx_energy_per_residue",
    "chain_length",
    "disulfide_count",
    "hbond_proxy_count",
    "sasa_total",
    "hbond_per_res",
    "disulfide_per_res",
    "sasa_per_res",
])

# Keep best available FoldX row per UniProt by lowest foldx_energy_per_residue
foldx_best = (
    foldx.sort_values(["uniprot_accession", "foldx_energy_per_residue"], ascending=[True, True])
    .groupby("uniprot_accession", as_index=False)
    .head(1)
    .copy()
)

foldx_keep = [
    "uniprot_accession",
    "structure_source",
    "structure_id_x",
    "foldx_input_pdb",
    "foldx_ready",
    "foldx_wt_total_energy",
    "chain_length",
    "disulfide_count",
    "hbond_proxy_count",
    "sasa_total",
    "hbond_per_res",
    "disulfide_per_res",
    "sasa_per_res",
    "foldx_energy_per_residue",
]
foldx_keep = [c for c in foldx_keep if c in foldx_best.columns]

foldx_best = foldx_best[foldx_keep].rename(columns={
    "structure_source": "foldx_structure_source",
    "structure_id_x": "foldx_structure_id",
})

df = df.merge(foldx_best, on="uniprot_accession", how="left")
df["has_foldx_result"] = df["foldx_energy_per_residue"].notna()

# ---------------------------------------------------------------------
# 6. ProtParam features
# ---------------------------------------------------------------------
prot = read_csv(PROTPARAM)
prot = as_numeric(prot, [
    "length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "isoelectric_point",
    "gravy",
    "helix_fraction",
    "turn_fraction",
    "sheet_fraction",
])

prot_keep = [
    "uniprot_accession",
    "thermal_label",
    "length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "isoelectric_point",
    "gravy",
    "helix_fraction",
    "turn_fraction",
    "sheet_fraction",
]
prot_keep = [c for c in prot_keep if c in prot.columns]

df = df.merge(prot[prot_keep].drop_duplicates("uniprot_accession"), on="uniprot_accession", how="left")

# Prefer explicit ProtParam length where present, otherwise sequence length
if "length" in df.columns:
    df["protparam_length"] = df["length"]
    df = df.drop(columns=["length"])
else:
    df["protparam_length"] = np.nan

# ---------------------------------------------------------------------
# 7. Signal peptide proxy
# ---------------------------------------------------------------------
sig = read_csv(SIGNAL)
sig = as_numeric(sig, [
    "nterm_hydrophobic_count_5_25",
    "nterm_positive_count_1_10",
])

sig_keep = [
    "uniprot_accession",
    "signal_peptide_proxy",
    "nterm_hydrophobic_count_5_25",
    "nterm_positive_count_1_10",
]
sig_keep = [c for c in sig_keep if c in sig.columns]

df = df.merge(sig[sig_keep].drop_duplicates("uniprot_accession"), on="uniprot_accession", how="left")

# ---------------------------------------------------------------------
# 8. ML predictions
# ---------------------------------------------------------------------
ml = read_csv(ML)
ml = as_numeric(ml, [
    "predicted_foldx_energy_per_residue",
    "prediction_error",
])

ml_keep = [
    "uniprot_accession",
    "predicted_foldx_energy_per_residue",
    "prediction_error",
]
ml_keep = [c for c in ml_keep if c in ml.columns]

df = df.merge(ml[ml_keep].drop_duplicates("uniprot_accession"), on="uniprot_accession", how="left")
df["has_ml_prediction"] = df["predicted_foldx_energy_per_residue"].notna()

# ---------------------------------------------------------------------
# 9. Sequence QC flags
# ---------------------------------------------------------------------
valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
df["sequence_has_X"] = df["sequence"].astype(str).str.contains("X", regex=False)
df["sequence_has_nonstandard_residue"] = df["sequence"].astype(str).apply(
    lambda s: any((aa not in valid_aas) for aa in s if aa.strip())
)
df["is_short_sequence_lt100aa"] = df["sequence_length"] < 100
df["is_long_sequence_gt1200aa"] = df["sequence_length"] > 1200

# ---------------------------------------------------------------------
# 10. Final column ordering
# ---------------------------------------------------------------------
preferred_cols = [
    "uniprot_accession",
    "primary_id",
    "source",
    "organism",
    "organism_type",
    "gh_family",
    "xref_cazy",
    "has_cazy_xref",
    "cazy_status",
    "sequence",
    "sequence_length",
    "protparam_length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "isoelectric_point",
    "gravy",
    "helix_fraction",
    "turn_fraction",
    "sheet_fraction",
    "thermal_label",
    "refseq_acc",
    "refseq_acc_all",
    "refseq_count",
    "has_refseq_link",
    "brenda_temperature_optimum",
    "brenda_temperature_range",
    "brenda_temperature_stability",
    "brenda_ph_optimum",
    "brenda_ph_range",
    "has_brenda_temp",
    "has_brenda_ph",
    "has_any_brenda",
    "signal_peptide_proxy",
    "nterm_hydrophobic_count_5_25",
    "nterm_positive_count_1_10",
    "has_any_structure",
    "has_experimental_structure",
    "has_modelled_structure",
    "structure_count",
    "structure_sources",
    "structure_ids",
    "representative_structure_id",
    "representative_structure_path",
    "modeller_model_count",
    "modeller_template_codes",
    "modeller_best_template_identity",
    "modeller_best_template_coverage",
    "modeller_run_status",
    "has_foldx_result",
    "foldx_structure_source",
    "foldx_structure_id",
    "foldx_input_pdb",
    "foldx_ready",
    "chain_length",
    "foldx_wt_total_energy",
    "foldx_energy_per_residue",
    "hbond_proxy_count",
    "disulfide_count",
    "sasa_total",
    "hbond_per_res",
    "disulfide_per_res",
    "sasa_per_res",
    "has_ml_prediction",
    "predicted_foldx_energy_per_residue",
    "prediction_error",
    "sequence_has_X",
    "sequence_has_nonstandard_residue",
    "is_short_sequence_lt100aa",
    "is_long_sequence_gt1200aa",
]

# Keep BRENDA comment/literature fields at the end
extra_cols = [c for c in df.columns if c not in preferred_cols]
df = df[[c for c in preferred_cols if c in df.columns] + extra_cols]

# Sort for readability
df = df.sort_values(["organism_type", "gh_family", "uniprot_accession"]).reset_index(drop=True)

OUT.parent.mkdir(parents=True, exist_ok=True)
SUMMARY.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(OUT, index=False)

# ---------------------------------------------------------------------
# 11. Summary report
# ---------------------------------------------------------------------
lines = []
lines.append("Xylanase Thesis Master Table Summary")
lines.append("=" * 50)
lines.append(f"Output file: {OUT}")
lines.append(f"Rows: {len(df)}")
lines.append(f"Unique UniProt accessions: {df['uniprot_accession'].nunique()}")
lines.append("")
lines.append("Organism type counts:")
lines.append(str(df["organism_type"].value_counts(dropna=False)))
lines.append("")
lines.append("GH family counts:")
lines.append(str(df["gh_family"].value_counts(dropna=False)))
lines.append("")
lines.append("Organism type x GH family:")
lines.append(str(pd.crosstab(df["organism_type"], df["gh_family"])))
lines.append("")
lines.append("Annotation/computational coverage:")
for col in [
    "has_refseq_link",
    "has_brenda_temp",
    "has_brenda_ph",
    "has_any_brenda",
    "signal_peptide_proxy",
    "has_any_structure",
    "has_experimental_structure",
    "has_modelled_structure",
    "has_foldx_result",
    "has_ml_prediction",
    "sequence_has_X",
    "sequence_has_nonstandard_residue",
    "is_short_sequence_lt100aa",
    "is_long_sequence_gt1200aa",
]:
    if col in df.columns:
        lines.append(f"{col}:")
        lines.append(str(df[col].value_counts(dropna=False)))
        lines.append("")

lines.append("Numeric summaries:")
for col in [
    "sequence_length",
    "protparam_length",
    "chain_length",
    "foldx_wt_total_energy",
    "foldx_energy_per_residue",
    "hbond_proxy_count",
    "disulfide_count",
    "sasa_total",
    "hbond_per_res",
    "disulfide_per_res",
    "sasa_per_res",
    "predicted_foldx_energy_per_residue",
    "prediction_error",
]:
    if col in df.columns:
        lines.append(f"{col}:")
        lines.append(str(pd.to_numeric(df[col], errors='coerce').describe()))
        lines.append("")

SUMMARY.write_text("\n".join(lines))

print(f"Saved: {OUT}")
print(f"Saved: {SUMMARY}")
print(f"Rows: {len(df)}")
print(f"Unique UniProt accessions: {df['uniprot_accession'].nunique()}")
print("\nOrganism x GH family:")
print(pd.crosstab(df["organism_type"], df["gh_family"]))
