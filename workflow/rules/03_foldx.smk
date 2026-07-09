rule run_foldx_analysis:
    input:
        structure="results/workflow/02_structure/structure_analysis_summary.csv"
    output:
        "results/workflow/03_foldx/foldx_analysis_summary.csv"
    log:
        "logs/snakemake/03_foldx_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage foldx \
            --output {output} \
            > {log} 2>&1
        """
