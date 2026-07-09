rule write_full_analysis_report:
    input:
        manifest="results/workflow/00_inputs/input_table_manifest.csv",
        sequence="results/workflow/01_sequence/sequence_analysis_summary.csv",
        structure="results/workflow/02_structure/structure_analysis_summary.csv",
        foldx="results/workflow/03_foldx/foldx_analysis_summary.csv",
        docking="results/workflow/04_docking/docking_analysis_summary.csv",
        machine_learning="results/workflow/05_ml/machine_learning_analysis_summary.csv",
        molecular_dynamics="results/workflow/06_md/molecular_dynamics_analysis_summary.csv",
        integration="results/workflow/07_integration/integration_analysis_summary.csv"
    output:
        "results/workflow/08_report/full_analysis_workflow_report.md"
    log:
        "logs/snakemake/08_full_analysis_report.log"
    shell:
        """
        python scripts/workflow/write_full_analysis_report.py \
            --config config/workflow_config.yaml \
            --input-manifest {input.manifest} \
            --sequence {input.sequence} \
            --structure {input.structure} \
            --foldx {input.foldx} \
            --docking {input.docking} \
            --machine-learning {input.machine_learning} \
            --molecular-dynamics {input.molecular_dynamics} \
            --integration {input.integration} \
            --output {output} \
            > {log} 2>&1
        """
