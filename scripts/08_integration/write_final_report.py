#!/usr/bin/env python3

from pathlib import Path
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Write final workflow report."
    )
    parser.add_argument("--curated-master", required=True)
    parser.add_argument("--sequence-features", required=True)
    parser.add_argument("--structure-features", required=True)
    parser.add_argument("--foldx-wt", required=True)
    parser.add_argument("--foldx-mutations", required=True)
    parser.add_argument("--docking", required=True)
    parser.add_argument("--ml-summary", required=True)
    parser.add_argument("--integration", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    curated = pd.read_csv(args.curated_master)
    sequence = pd.read_csv(args.sequence_features)
    structure = pd.read_csv(args.structure_features)
    foldx_wt = pd.read_csv(args.foldx_wt)
    foldx_mutations = pd.read_csv(args.foldx_mutations)
    docking = pd.read_csv(args.docking)
    ml_summary = pd.read_csv(args.ml_summary)
    integration = pd.read_csv(args.integration)

    docking_ready = docking[docking["status"] == "ready"]
    top = integration.head(10)

    lines = [
        "# Final workflow report",
        "",
        "This report was generated automatically from the Snakemake workflow.",
        "",
        "## Workflow mode",
        "",
        "The current run is a test-mode reproducibility run. FoldX and docking were executed on a small balanced subset of bacterial GH10, bacterial GH11, fungal GH10, and fungal GH11 enzymes.",
        "",
        "## Output summary",
        "",
        f"- Curated enzyme records: {len(curated)}",
        f"- Sequence-feature rows: {len(sequence)}",
        f"- Structure-feature rows: {len(structure)}",
        f"- FoldX WT stability rows: {len(foldx_wt)}",
        f"- FoldX mutation rows: {len(foldx_mutations)}",
        f"- Docking rows: {len(docking)}",
        f"- Successful docking rows: {len(docking_ready)}",
        f"- ML summary rows: {len(ml_summary)}",
        f"- Integrated candidate rows: {len(integration)}",
        "",
        "## Group coverage in integration table",
        "",
        integration.groupby(["organism_type", "gh_family"]).size().reset_index(name="rows").to_markdown(index=False),
        "",
        "## Docking score summary",
        "",
        docking_ready["vina_affinity"].describe().to_markdown(),
        "",
        "## Top integrated candidates",
        "",
        top[
            [
                "rank",
                "uniprot_accession",
                "organism_type",
                "gh_family",
                "pdb_id",
                "foldx_mutation_code",
                "foldx_ddg",
                "mean_vina_affinity",
                "best_vina_affinity",
                "foldx_class",
                "integrated_score",
                "candidate_category",
            ]
        ].to_markdown(index=False),
        "",
        "## ML branch note",
        "",
        "The ML branch executed successfully, but in test mode the dataset is too small for thesis-scale interpretation. Final ML interpretation should use the full FoldX and docking tables.",
        "",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Saved final report: {output_path}")


if __name__ == "__main__":
    main()
