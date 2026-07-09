rule run_integration_analysis:
    input:
        ml="results/workflow/05_ml/machine_learning_analysis_summary.csv",
        md="results/workflow/06_md/molecular_dynamics_analysis_summary.csv"
    output:
        "results/workflow/07_integration/integration_analysis_summary.csv"
    log:
        "logs/snakemake/07_integration_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage integration \
            --output {output} \
            > {log} 2>&1
        """
