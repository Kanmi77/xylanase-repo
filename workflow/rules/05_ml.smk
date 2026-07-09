rule run_machine_learning_analysis:
    input:
        docking="results/workflow/04_docking/docking_analysis_summary.csv"
    output:
        "results/workflow/05_ml/machine_learning_analysis_summary.csv"
    log:
        "logs/snakemake/05_machine_learning_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage machine_learning \
            --output {output} \
            > {log} 2>&1
        """
