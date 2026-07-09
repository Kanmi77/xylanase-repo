rule run_structure_analysis:
    input:
        sequence="results/workflow/01_sequence/sequence_analysis_summary.csv"
    output:
        "results/workflow/02_structure/structure_analysis_summary.csv"
    log:
        "logs/snakemake/02_structure_analysis.log"
    shell:
        """
        python scripts/workflow/run_analysis_stage.py \
            --config config/workflow_config.yaml \
            --stage structure \
            --output {output} \
            > {log} 2>&1
        """
