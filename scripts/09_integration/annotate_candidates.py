#!/usr/bin/env python3
# Purpose: Annotate integrated candidates.

from pathlib import Path
import argparse
import re
import pandas as pd
import numpy as np

ROOT = Path(".").resolve()

DEFAULT_MASTER_CANDIDATES = [
    "data/curated/xylanase_master_all_curated_with_brenda.csv",
    "data/curated/xylanase_master_all_curated.csv",
    "data/curated/xylanase_master_deduplicated.csv",
    "data/curated/xylanase_master.csv",
]

DEFAULT_INTEGRATION_DIR = "results/integration/stratified_pdb_vs_homology_new_ml_pdb_foldx_docking"

TABLES_TO_ANNOTATE = [
    "final_evidence_tiered_summary.csv",
    "pdb_branch_integrated.csv",
    "homology_branch_integrated.csv",
    "top15_pdb_branch.csv",
    "top15_homology_branch.csv",
]


def read_table(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")

    return pd.read_csv(path)


def find_existing_master():
    for p in DEFAULT_MASTER_CANDIDATES:
        path = ROOT / p
        if path.exists():
            return path
    return None


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

    # Remove common candidate mutation suffix if accidentally present.
    # Example: D5UGW9_LA23A -> D5UGW9
    if "_" in x:
        first = x.split("_")[0]
        if re.match(r"^[A-Z0-9]{5,10}$", first):
            x = first

    return x


def first_nonnull(series):
    for v in series:
        if pd.notna(v) and str(v).strip() not in {"", "nan", "None"}:
            return v
    return np.nan


def build_master_mapping(master):
    master = master.copy()

    acc_col = find_col(master, [
        r"^uniprot_accession$",
        r"uniprot.*accession",
        r"^accession$",
        r"entry",
        r"protein_accession",
    ])

    if not acc_col:
        raise SystemExit("Could not identify accession column in master dataset.")

    enzyme_col = find_col(master, [
        r"protein.*name",
        r"enzyme.*name",
        r"recommended.*name",
        r"full.*name",
        r"description",
        r"product",
        r"function",
    ])

    organism_col = find_col(master, [
        r"^organism$",
        r"organism.*name",
        r"scientific.*name",
        r"taxonomic.*name",
        r"species",
    ])

    organism_type_col = find_col(master, [
        r"organism.*type",
        r"source.*type",
        r"taxonomy.*group",
        r"kingdom",
    ])

    gh_col = find_col(master, [
        r"^gh_family$",
        r"gh.*family",
        r"query_gh_family",
        r"family",
    ])

    ec_col = find_col(master, [
        r"ec_number",
        r"ec.*number",
        r"enzyme.*commission",
    ])

    gene_col = find_col(master, [
        r"gene.*name",
        r"gene$",
    ])

    master["candidate_accession"] = master[acc_col].apply(clean_accession)

    cols = ["candidate_accession"]

    rename = {}

    if enzyme_col:
        cols.append(enzyme_col)
        rename[enzyme_col] = "candidate_enzyme_name"

    if organism_col:
        cols.append(organism_col)
        rename[organism_col] = "candidate_organism"

    if organism_type_col:
        cols.append(organism_type_col)
        rename[organism_type_col] = "candidate_organism_type"

    if gh_col:
        cols.append(gh_col)
        rename[gh_col] = "candidate_gh_family"

    if ec_col:
        cols.append(ec_col)
        rename[ec_col] = "candidate_ec_number"

    if gene_col:
        cols.append(gene_col)
        rename[gene_col] = "candidate_gene_name"

    mapping = master[cols].rename(columns=rename)

    # Collapse duplicates by accession.
    agg = {}
    for c in mapping.columns:
        if c != "candidate_accession":
            agg[c] = first_nonnull

    mapping = (
        mapping[mapping["candidate_accession"].notna()]
        .groupby("candidate_accession", dropna=False)
        .agg(agg)
        .reset_index()
    )

    detected = {
        "accession_col": acc_col,
        "enzyme_col": enzyme_col,
        "organism_col": organism_col,
        "organism_type_col": organism_type_col,
        "gh_col": gh_col,
        "ec_col": ec_col,
        "gene_col": gene_col,
        "master_rows": len(master),
        "mapping_rows": len(mapping),
    }

    return mapping, detected


def annotate_table(df, mapping):
    df = df.copy()

    # Prefer standardised accession if present.
    acc_source_col = None

    for c in [
        "candidate_accession",
        "std_uniprot",
        "uniprot_accession",
        "protein",
        "accession",
    ]:
        if c in df.columns:
            acc_source_col = c
            break

    if not acc_source_col:
        df["candidate_accession"] = np.nan
    else:
        df["candidate_accession"] = df[acc_source_col].apply(clean_accession)

    out = df.merge(mapping, on="candidate_accession", how="left")

    # Create a cleaner display label.
    if "mutation" in out.columns:
        mut = out["mutation"].astype(str)
    elif "std_mutation" in out.columns:
        mut = out["std_mutation"].astype(str)
    else:
        mut = pd.Series("", index=out.index)

    out["candidate_label_clean"] = out["candidate_accession"].fillna("").astype(str)

    has_mut = mut.notna() & ~mut.isin(["nan", "None", ""])
    out.loc[has_mut, "candidate_label_clean"] = (
        out.loc[has_mut, "candidate_accession"].astype(str) + "_" + mut.loc[has_mut].astype(str)
    )

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", default=None)
    parser.add_argument("--integration-dir", default=DEFAULT_INTEGRATION_DIR)
    args = parser.parse_args()

    if args.master:
        master_path = ROOT / args.master
    else:
        master_path = find_existing_master()

    if master_path is None or not master_path.exists():
        raise SystemExit("No master dataset found. Use --master path/to/master.csv")

    integration_dir = ROOT / args.integration_dir

    if not integration_dir.exists():
        raise SystemExit(f"Integration directory not found: {integration_dir}")

    master = read_table(master_path)
    mapping, detected = build_master_mapping(master)

    report_lines = []
    report_lines.append("# Integrated candidate accession annotation report\n")
    report_lines.append(f"Master file: `{master_path.relative_to(ROOT)}`\n")
    report_lines.append("## Detected master columns\n")
    for k, v in detected.items():
        report_lines.append(f"- `{k}`: `{v}`")

    report_lines.append("\n## Annotated integration tables\n")

    for table in TABLES_TO_ANNOTATE:
        in_path = integration_dir / table

        if not in_path.exists():
            report_lines.append(f"- Missing table, skipped: `{in_path.relative_to(ROOT)}`")
            continue

        df = read_table(in_path)
        annotated = annotate_table(df, mapping)

        out_path = integration_dir / table.replace(".csv", "_annotated.csv")
        annotated.to_csv(out_path, index=False)

        n = len(annotated)
        acc_n = annotated["candidate_accession"].notna().sum()

        enzyme_n = annotated["candidate_enzyme_name"].notna().sum() if "candidate_enzyme_name" in annotated.columns else 0
        organism_n = annotated["candidate_organism"].notna().sum() if "candidate_organism" in annotated.columns else 0
        gh_n = annotated["candidate_gh_family"].notna().sum() if "candidate_gh_family" in annotated.columns else 0

        report_lines.append(f"\n### `{table}`")
        report_lines.append(f"- rows: {n}")
        report_lines.append(f"- accession mapped: {acc_n}")
        report_lines.append(f"- enzyme/protein name mapped: {enzyme_n}")
        report_lines.append(f"- organism mapped: {organism_n}")
        report_lines.append(f"- GH family mapped: {gh_n}")
        report_lines.append(f"- output: `{out_path.relative_to(ROOT)}`")

        missing = annotated[
            annotated["candidate_accession"].notna()
            & (
                annotated["candidate_organism"].isna()
                if "candidate_organism" in annotated.columns
                else True
            )
        ]["candidate_accession"].dropna().unique().tolist()

        if missing:
            report_lines.append(f"- unmapped organism accessions preview: `{missing[:20]}`")

    report_path = integration_dir / "ACCESSION_ENZYME_ORGANISM_ANNOTATION_REPORT.md"

    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
        f.write("\n")

    print(f"Master file used: {master_path}")
    print(f"Saved annotation report: {report_path}")
    print()
    print("\n".join(report_lines[:80]))


if __name__ == "__main__":
    main()
