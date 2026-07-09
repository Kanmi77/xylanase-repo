rule run_docking_analysis:
    input:
        foldx="results/workflow/03_foldx/foldx_analysis_summary.csv"
    output:
        "results/workflow/04_docking/docking_analysis_summary.csv"
    log:
        "logs/snakemake/04_docking_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage docking \
            --output {output} \
            > {log} 2>&1
        """
