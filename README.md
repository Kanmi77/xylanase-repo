# Xylanase Structural Thermostability Assessment

This repository contains the source code, workflow configuration files, and supplementary result tables associated with the thesis:

**In Silico Structural Thermostability Assessment of Bacterial and Fungal Thermostable Xylanases through Bioinformatics and Machine Learning Predictions**

## Overview

This project investigates bacterial and fungal GH10 and GH11 xylanases using a fully in-silico bioinformatics workflow. The analysis includes data curation, sequence analysis, structural feature extraction, FoldX stability assessment, mutation screening, molecular docking, machine learning, and integrated candidate prioritisation.

## Repository Structure

| Folder | Description |
|---|---|
| `config/` | Configuration files for the Snakemake workflow |
| `workflow/` | Main Snakemake workflow files |
| `workflow/rules/` | Modular Snakemake rule files |
| `scripts/` | Python and shell scripts used by the workflow |
| `results/tables/` | Main thesis result tables |
| `results/source_tables/` | Supporting source result tables |
| `Supplementary Material/` | Supplementary thesis data arranged for external access |

## Supplementary Material

The supplementary material associated with this thesis is available in the `Supplementary Material` folder of this repository.

It includes the main thesis result tables and supporting source tables used for the downstream analyses presented in this study.

The main result tables are accessible here:

[Main_Thesis_Result_Tables](Supplementary%20Material/Main_Thesis_Result_Tables)

The supporting source result tables are accessible here:

[Source_Result_Tables](Supplementary%20Material/Source_Result_Tables)

## Code and Data Availability

The source code for the xylanase structural thermostability assessment workflow developed as part of this thesis is available in this GitHub repository.

The repository includes source code, configuration files, Snakemake workflow files, supplementary result tables, and documentation required to understand and reproduce the computational workflow.

Main workflow files:

- `workflow/Snakefile`
- `workflow/Snakefile.full`
- `workflow/rules/`
- `scripts/`
- `config/workflow_config.yaml`

## Workflow Check

After setting up the required tools and paths, the workflow can be checked with:

snakemake -s workflow/Snakefile -n -p --cores 4

Large generated runtime files, downloaded PDB structures, raw downloads, FoldX working directories, and Vina intermediate receptor/output files are excluded from version control.
