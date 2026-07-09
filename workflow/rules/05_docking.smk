rule run_docking_screen:
    input:
        mutations=config["outputs"]["foldx_mutations"]
    output:
        config["outputs"]["docking_summary"]
    log:
        "logs/snakemake/05_run_docking_screen.log"
    shell:
        """
        python scripts/05_docking/run_docking_screen.py \
            --config config/workflow_config.yaml \
            --input {input.mutations} \
            --output {output} \
            > {log} 2>&1
        """
