rule download_uniprot_records:
    output:
        "data/raw/uniprot_records.tsv"
    log:
        "logs/snakemake/01_fetch_uniprot_records.log"
    shell:
        """
        python scripts/01_curation/fetch_uniprot_records.py \
            --config config/workflow_config.yaml \
            --output {output} \
            > {log} 2>&1
        """


rule curate_xylanase_dataset:
    input:
        raw_records="data/raw/uniprot_records.tsv"
    output:
        config["outputs"]["curated_master"]
    log:
        "logs/snakemake/01_curate_master_dataset.log"
    shell:
        """
        python scripts/01_curation/curate_master_dataset.py \
            --config config/workflow_config.yaml \
            --input {input.raw_records} \
            --output {output} \
            > {log} 2>&1
        """
