#!/usr/bin/env python3

from pathlib import Path
import argparse
import re
import pandas as pd
import numpy as np

ROOT = Path(".").resolve()

METADATA_SOURCES = [
    "data/curated/xylanase_master_uniprot.csv",
    "data/curated/xylanase_master_all_curated_universal_full.csv",
    "data/curated/uniprot_xylanase_curated_with_cazy.csv",
    "data/curated/uniprot_xylanase_curated.csv",
    "data/curated/xylanase_master_all_curated_plus_brenda_missing.csv",
    "data/curated/xylanase_master_all_curated_with_brenda.csv",
    "data/curated/xylanase_master_all_curated.csv",
    "data/curated/xylanase_master_deduplicated.csv",
    "data/curated/xylanase_structured_subset.csv",
    "data/curated/xylanase_structured_subset_with_foldx.csv",
    "data/curated/xylanase_structured_subset_with_foldx_norm.csv",
    "data/curated/pdb_inventory.csv",
    "results/optionC_original_only/pdb60_final_candidate_table/pdb60_final_corrected_candidate_table.csv",
    "results/optionC_original_only/pdb60_final_ranked_no_tiers/pdb60_final_ranked_no_tiers.csv",
]

TABLES = [
    "final_evidence_tiered_summary_annotated.csv",
    "pdb_branch_integrated_annotated.csv",
    "homology_branch_integrated_annotated.csv",
    "top15_pdb_branch_annotated.csv",
    "top15_homology_branch_annotated.csv",
    "final_evidence_tiered_summary.csv",
    "pdb_branch_integrated.csv",
    "homology_branch_integrated.csv",
    "top15_pdb_branch.csv",
    "top15_homology_branch.csv",
]

DEFAULT_DIR = "results/integration/stratified_pdb_vs_homology_new_ml_pdb_foldx_docking"


def read_csv(path):
    path = Path(path)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def find_col(df, patterns):
    for pat in patterns:
        for c in df.columns:
            if re.search(pat, c, flags=re.I):
                return c
    return None


def clean_accession(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip()

    if x.lower() in {"nan", "none", ""}:
        return np.nan

    if "_" in x:
        first = x.split("_")[0]
        if re.match(r"^[A-Z0-9]{5,10}$", first):
            return first

    return x


def first_nonnull(series):
    for v in series:
        if pd.notna(v) and str(v).strip() not in {"", "nan", "None"}:
            return v
    return np.nan


def derive_enzyme_label(gh):
    if pd.isna(gh):
        return "xylanase"

    gh = str(gh).upper()

    if "GH10" in gh:
        return "GH10 endo-beta-1,4-xylanase"

    if "GH11" in gh:
        return "GH11 endo-beta-1,4-xylanase"

    return "xylanase"


def build_metadata_mapping():
    records = []
    logs = []

    for rel in METADATA_SOURCES:
        path = ROOT / rel

        if not path.exists():
            logs.append((rel, "missing", "file not found"))
            continue

        try:
            df = read_csv(path)
        except Exception as e:
            logs.append((rel, "failed", str(e)))
            continue

        acc_col = find_col(df, [
            r"^uniprot_accession$",
            r"uniprot.*accession",
            r"^accession$",
            r"entry",
            r"^protein$",
        ])

        if not acc_col:
            logs.append((rel, "skipped", "no accession column"))
            continue

        enzyme_col = find_col(df, [
            r"protein.*name",
            r"enzyme.*name",
            r"recommended.*name",
            r"description",
            r"product",
            r"function",
        ])

        organism_col = find_col(df, [
            r"^organism$",
            r"organism.*name",
            r"scientific.*name",
            r"species",
            r"query_organism",
            r"reference_organism",
        ])

        organism_type_col = find_col(df, [
            r"organism.*type",
            r"query_organism_type",
            r"reference_organism_type",
            r"kingdom",
        ])

        gh_col = find_col(df, [
            r"^gh_family$",
            r"gh.*family",
            r"query_gh_family",
            r"reference_gh_family",
            r"family",
        ])

        ec_col = find_col(df, [
            r"^ec$",
            r"ec_number",
            r"ec.*number",
        ])

        temp = pd.DataFrame()
        temp["candidate_accession"] = df[acc_col].apply(clean_accession)
        temp["metadata_source"] = rel
        temp["candidate_enzyme_name_source"] = df[enzyme_col] if enzyme_col else np.nan
        temp["candidate_organism_source"] = df[organism_col] if organism_col else np.nan
        temp["candidate_organism_type_source"] = df[organism_type_col] if organism_type_col else np.nan
        temp["candidate_gh_family_source"] = df[gh_col] if gh_col else np.nan
        temp["candidate_ec_number_source"] = df[ec_col] if ec_col else np.nan

        temp = temp[temp["candidate_accession"].notna()]
        records.append(temp)

        logs.append((
            rel,
            "used",
            f"rows={len(df)}, accession={acc_col}, enzyme={enzyme_col}, organism={organism_col}, organism_type={organism_type_col}, gh={gh_col}, ec={ec_col}"
        ))

    if not records:
        raise SystemExit("No metadata source files could be used.")

    pool = pd.concat(records, ignore_index=True)

    mapping = (
        pool.groupby("candidate_accession", dropna=False)
        .agg({
            "candidate_enzyme_name_source": first_nonnull,
            "candidate_organism_source": first_nonnull,
            "candidate_organism_type_source": first_nonnull,
            "candidate_gh_family_source": first_nonnull,
            "candidate_ec_number_source": first_nonnull,
            "metadata_source": first_nonnull,
        })
        .reset_index()
    )

    return mapping, logs


def get_candidate_accession(df):
    for c in ["candidate_accession", "std_uniprot", "uniprot_accession", "protein", "accession"]:
        if c in df.columns:
            return df[c].apply(clean_accession)

    return pd.Series(np.nan, index=df.index)


def fill_table(df, mapping):
    df = df.copy()

    df["candidate_accession"] = get_candidate_accession(df)

    out = df.merge(mapping, on="candidate_accession", how="left")

    if "candidate_enzyme_name" not in out.columns:
        out["candidate_enzyme_name"] = np.nan

    if "candidate_organism" not in out.columns:
        out["candidate_organism"] = np.nan

    if "candidate_organism_type" not in out.columns:
        out["candidate_organism_type"] = np.nan

    if "candidate_gh_family" not in out.columns:
        out["candidate_gh_family"] = np.nan

    if "candidate_ec_number" not in out.columns:
        out["candidate_ec_number"] = np.nan

    out["candidate_enzyme_name"] = out["candidate_enzyme_name"].combine_first(out["candidate_enzyme_name_source"])
    out["candidate_organism"] = out["candidate_organism"].combine_first(out["candidate_organism_source"])
    out["candidate_organism_type"] = out["candidate_organism_type"].combine_first(out["candidate_organism_type_source"])
    out["candidate_gh_family"] = out["candidate_gh_family"].combine_first(out["candidate_gh_family_source"])
    out["candidate_ec_number"] = out["candidate_ec_number"].combine_first(out["candidate_ec_number_source"])

    out["candidate_enzyme_label"] = out["candidate_gh_family"].apply(derive_enzyme_label)
    out["candidate_enzyme_display"] = out["candidate_enzyme_name"].combine_first(out["candidate_enzyme_label"])

    if "mutation" in out.columns:
        mut = out["mutation"].astype(str)
    elif "std_mutation" in out.columns:
        mut = out["std_mutation"].astype(str)
    else:
        mut = pd.Series("", index=out.index)

    out["candidate_label_clean"] = out["candidate_accession"].fillna("").astype(str)

    has_mut = mut.notna() & ~mut.isin(["nan", "None", ""])
    out.loc[has_mut, "candidate_label_clean"] = (
        out.loc[has_mut, "candidate_accession"].astype(str)
        + "_"
        + mut.loc[has_mut].astype(str)
    )

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--integration-dir", default=DEFAULT_DIR)
    args = parser.parse_args()

    integration_dir = ROOT / args.integration_dir

    if not integration_dir.exists():
        raise SystemExit(f"Integration directory not found: {integration_dir}")

    mapping, logs = build_metadata_mapping()

    report = []
    report.append("# Metadata fallback fill report\n")
    report.append("This fills enzyme/protein name, organism, organism type, GH family, and EC metadata from broad curated UniProt/xylanase files.\n")

    report.append("## Metadata sources\n")
    for rel, status, msg in logs:
        report.append(f"- `{rel}`: {status} — {msg}")

    report.append("\n## Output tables\n")

    written = set()

    for table in TABLES:
        in_path = integration_dir / table

        if not in_path.exists():
            continue

        # Avoid processing both raw and annotated into same output name accidentally.
        out_name = table.replace(".csv", "_metadata_filled.csv")
        out_path = integration_dir / out_name

        if str(out_path) in written:
            continue

        df = read_csv(in_path)
        filled = fill_table(df, mapping)

        filled.to_csv(out_path, index=False)
        written.add(str(out_path))

        n = len(filled)

        enzyme_name_n = filled["candidate_enzyme_name"].notna().sum()
        enzyme_display_n = filled["candidate_enzyme_display"].notna().sum()
        organism_n = filled["candidate_organism"].notna().sum()
        organism_type_n = filled["candidate_organism_type"].notna().sum()
        gh_n = filled["candidate_gh_family"].notna().sum()
        ec_n = filled["candidate_ec_number"].notna().sum()

        report.append(f"\n### `{table}`")
        report.append(f"- rows: {n}")
        report.append(f"- exact enzyme/protein names mapped: {enzyme_name_n}")
        report.append(f"- enzyme display available: {enzyme_display_n}")
        report.append(f"- organism mapped: {organism_n}")
        report.append(f"- organism type mapped: {organism_type_n}")
        report.append(f"- GH family mapped: {gh_n}")
        report.append(f"- EC mapped: {ec_n}")
        report.append(f"- output: `{out_path.relative_to(ROOT)}`")

        missing_org = filled[
            filled["candidate_accession"].notna() & filled["candidate_organism"].isna()
        ]["candidate_accession"].dropna().unique().tolist()

        if missing_org:
            report.append(f"- still missing organism preview: `{missing_org[:30]}`")

    report_path = integration_dir / "METADATA_FALLBACK_FILL_REPORT.md"

    with open(report_path, "w") as f:
        f.write("\n".join(report))
        f.write("\n")

    print(f"Saved fallback report: {report_path}")
    print()
    print("\n".join(report[:160]))


if __name__ == "__main__":
    main()
