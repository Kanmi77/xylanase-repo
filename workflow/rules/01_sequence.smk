rule run_sequence_analysis:
    input:
        manifest="results/workflow/00_inputs/input_table_manifest.csv"
    output:
        "results/workflow/01_sequence/sequence_analysis_summary.csv"
    log:
        "logs/snakemake/01_sequence_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage sequence \
            --output {output} \
            > {log} 2>&1
        """
