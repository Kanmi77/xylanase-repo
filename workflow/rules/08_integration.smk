rule integrate_candidates:
    input:
        foldx_mutations=config["outputs"]["foldx_mutations"],
        docking=config["outputs"]["docking_summary"],
        ml_summary=config["outputs"]["ml_summary"]
    output:
        config["outputs"]["integration_scores"]
    log:
        "logs/snakemake/08_integrate_candidates.log"
    shell:
        """
        python scripts/08_integration/integrate_candidates.py \
            --foldx-mutations {input.foldx_mutations} \
            --docking {input.docking} \
            --ml-summary {input.ml_summary} \
            --output {output} \
            > {log} 2>&1
        """

rule write_final_report:
    input:
        curated_master=config["outputs"]["curated_master"],
        sequence_features=config["outputs"]["sequence_features"],
        structure_features=config["outputs"]["structure_features"],
        foldx_wt=config["outputs"]["foldx_wt"],
        foldx_mutations=config["outputs"]["foldx_mutations"],
        docking=config["outputs"]["docking_summary"],
        ml_summary=config["outputs"]["ml_summary"],
        integration=config["outputs"]["integration_scores"]
    output:
        config["outputs"]["final_report"]
    log:
        "logs/snakemake/08_write_final_report.log"
    shell:
        """
        python scripts/08_integration/write_final_report.py \
            --curated-master {input.curated_master} \
            --sequence-features {input.sequence_features} \
            --structure-features {input.structure_features} \
            --foldx-wt {input.foldx_wt} \
            --foldx-mutations {input.foldx_mutations} \
            --docking {input.docking} \
            --ml-summary {input.ml_summary} \
            --integration {input.integration} \
            --output {output} \
            > {log} 2>&1
        """
