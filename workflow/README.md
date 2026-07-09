# Full Snakemake analysis workflow

This workflow reproduces the thesis analysis at the organised result-table level.

## Stages

| Stage | Rule file | Output |
|---|---|---|
| Inputs | `rules/00_inputs.smk` | Input table manifest |
| Sequence | `rules/01_sequence.smk` | Sequence-analysis summary |
| Structure | `rules/02_structure.smk` | Structure-analysis summary |
| FoldX | `rules/03_foldx.smk` | FoldX stability and mutation summary |
| Docking | `rules/04_docking.smk` | Docking-analysis summary |
| Machine learning | `rules/05_ml.smk` | ML-analysis summary |
| Molecular dynamics | `rules/06_md.smk` | MD-analysis summary |
| Integration | `rules/07_integration.smk` | Candidate-prioritisation summary |
| Report | `rules/08_report.smk` | Full workflow report |

## Changing enzyme type

Edit:

`config/workflow_config.yaml`

Important fields:

- `project.enzyme_name`
- `project.enzyme_key`
- `project.families`
- `project.organism_groups`
- `tables`

## Run

Dry run:

`snakemake -s workflow/Snakefile -n -p --cores 4`

Real run:

`snakemake -s workflow/Snakefile -p --cores 4`
