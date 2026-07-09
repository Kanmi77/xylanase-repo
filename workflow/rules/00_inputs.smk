rule build_input_manifest:
    input:
        config_file="config/workflow_config.yaml"
    output:
        "results/workflow/00_inputs/input_table_manifest.csv"
    log:
        "logs/snakemake/00_build_input_manifest.log"
    shell:
        """
        python scripts/workflow/build_input_manifest.py \
            --config {input.config_file} \
            --output {output} \
            > {log} 2>&1
        """
