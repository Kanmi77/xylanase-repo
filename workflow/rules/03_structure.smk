rule build_structure_inventory:
    input:
        master=config["outputs"]["curated_master"]
    output:
        "results/03_structure/structure_inventory.csv"
    log:
        "logs/snakemake/03_build_structure_inventory.log"
    shell:
        """
        python scripts/03_structure/build_structure_inventory.py \
            --config config/workflow_config.yaml \
            --input {input.master} \
            --output {output} \
            > {log} 2>&1
        """


rule compute_structure_features:
    input:
        inventory="results/03_structure/structure_inventory.csv"
    output:
        config["outputs"]["structure_features"]
    log:
        "logs/snakemake/03_compute_structure_features.log"
    shell:
        """
        python scripts/03_structure/compute_structure_features.py \
            --input {input.inventory} \
            --output {output} \
            > {log} 2>&1
        """
