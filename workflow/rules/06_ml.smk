rule run_machine_learning:
    input:
        foldx_wt=config["outputs"]["foldx_wt"],
        foldx_mutations=config["outputs"]["foldx_mutations"],
        docking=config["outputs"]["docking_summary"]
    output:
        config["outputs"]["ml_summary"]
    log:
        "logs/snakemake/06_run_machine_learning.log"
    shell:
        """
        python scripts/06_ml/run_machine_learning.py \
            --config config/workflow_config.yaml \
            --foldx-wt {input.foldx_wt} \
            --foldx-mutations {input.foldx_mutations} \
            --docking {input.docking} \
            --output {output} \
            > {log} 2>&1
        """
