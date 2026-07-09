rule run_molecular_dynamics_analysis:
    input:
        ml="results/workflow/05_ml/machine_learning_analysis_summary.csv"
    output:
        "results/workflow/06_md/molecular_dynamics_analysis_summary.csv"
    log:
        "logs/snakemake/06_molecular_dynamics_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage molecular_dynamics \
            --output {output} \
            > {log} 2>&1
        """
