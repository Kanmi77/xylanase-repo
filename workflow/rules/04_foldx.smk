rule run_foldx_wt_stability:
    input:
        structure_features=config["outputs"]["structure_features"]
    output:
        config["outputs"]["foldx_wt"]
    log:
        "logs/snakemake/04_run_foldx_wt_stability.log"
    shell:
        """
        python scripts/04_foldx/run_foldx_wt_stability.py \
            --config config/workflow_config.yaml \
            --input {input.structure_features} \
            --output {output} \
            > {log} 2>&1
        """


rule run_foldx_mutation_screen:
    input:
        wt=config["outputs"]["foldx_wt"]
    output:
        config["outputs"]["foldx_mutations"]
    log:
        "logs/snakemake/04_run_foldx_mutation_screen.log"
    shell:
        """
        python scripts/04_foldx/run_foldx_mutation_screen.py \
            --config config/workflow_config.yaml \
            --input {input.wt} \
            --output {output} \
            > {log} 2>&1
        """
