rule calculate_sequence_features:
    input:
        master=config["outputs"]["curated_master"]
    output:
        config["outputs"]["sequence_features"]
    log:
        "logs/snakemake/02_compute_sequence_features.log"
    shell:
        """
        python scripts/02_sequence/compute_sequence_features.py \
            --input {input.master} \
            --output {output} \
            > {log} 2>&1
        """


rule prepare_group_fastas:
    input:
        master=config["outputs"]["curated_master"]
    output:
        directory("results/02_sequence/group_fastas")
    log:
        "logs/snakemake/02_prepare_group_fastas.log"
    shell:
        """
        python scripts/02_sequence/prepare_group_fastas.py \
            --config config/workflow_config.yaml \
            --input {input.master} \
            --output-dir {output} \
            > {log} 2>&1
        """


rule align_sequences:
    input:
        fasta_dir="results/02_sequence/group_fastas"
    output:
        directory("results/02_sequence/alignments")
    threads:
        config["sequence_analysis"]["alignment_threads"]
    log:
        "logs/snakemake/02_run_mafft_alignments.log"
    shell:
        """
        python scripts/02_sequence/run_mafft_alignments.py \
            --input-dir {input.fasta_dir} \
            --output-dir {output} \
            --threads {threads} \
            > {log} 2>&1
        """


rule build_phylogenetic_trees:
    input:
        alignment_dir="results/02_sequence/alignments"
    output:
        directory("results/02_sequence/trees")
    log:
        "logs/snakemake/02_run_fasttree.log"
    shell:
        """
        python scripts/02_sequence/run_fasttree.py \
            --input-dir {input.alignment_dir} \
            --output-dir {output} \
            > {log} 2>&1
        """


rule summarize_conserved_positions:
    input:
        alignment_dir="results/02_sequence/alignments",
        trees="results/02_sequence/trees"
    output:
        config["outputs"]["conservation_summary"]
    log:
        "logs/snakemake/02_compute_conservation.log"
    shell:
        """
        python scripts/02_sequence/compute_conservation.py \
            --config config/workflow_config.yaml \
            --input-dir {input.alignment_dir} \
            --output {output} \
            > {log} 2>&1
        """
