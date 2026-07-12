# Xylanase Structural Thermostability Workflow

This repository contains the reproducible Snakemake workflow developed for the thesis:

**In Silico Structural Thermostability Assessment of Bacterial and Fungal Thermostable Xylanases through Bioinformatics and Machine Learning Predictions**

The workflow supports in-silico thermostability assessment of bacterial and fungal GH10/GH11 xylanases using sequence analysis, structural feature extraction, FoldX stability assessment, mutation screening, molecular docking, machine learning, and integrated candidate ranking.

The final thesis result tables, full datasets, figures, visualizations, and supplementary materials are maintained separately in the results repository:

https://github.com/Kanmi77/xylanase-thermostability-results

## Repository Purpose

This repository is intended for the computational workflow, scripts, configuration files, and environment definition required to reproduce or adapt the xylanase thermostability analysis.

It contains:

- Snakemake workflow files
- modular workflow rules
- Python and shell scripts used by the workflow
- workflow configuration files
- Conda environment definition
- tool-path configuration files
- documentation for running the workflow

Generated result outputs, large runtime files, raw downloads, FoldX working directories, and docking intermediate files are excluded from version control.

## Workflow Overview

The workflow is organized into the following analysis stages:

1. UniProt xylanase record retrieval and dataset curation
2. sequence-feature calculation
3. multiple sequence alignment using MAFFT
4. phylogenetic reconstruction using FastTree
5. conserved-position analysis
6. structure inventory and structural-feature extraction
7. FoldX wild-type stability calculation
8. FoldX mutation screening
9. receptor preparation and AutoDock Vina docking
10. machine-learning analysis
11. integrated candidate ranking

## Repository Structure

| Folder/File | Description |
|---|---|
| `config/` | Workflow configuration files |
| `environments/` | Conda environment file for active workflow dependencies |
| `workflow/` | Main Snakemake workflow files |
| `workflow/rules/` | Modular Snakemake rule files |
| `scripts/` | Python and shell scripts executed by the workflow |
| `data/` | Lightweight input templates or manually supplied input files |
| `docs/` | Workflow notes and supporting documentation |
| `.gitignore` | Excludes generated outputs, logs, raw downloads, and runtime files |

## Tools, Packages & Versions

| Tool / Package     | Version / Configuration | Role / Purpose                                      |
|--------------------|-------------------------|-----------------------------------------------------|
| **Python**         | 3.10                    | Main scripting language                             |
| **Snakemake**      | 7.32.4                  | Workflow management                                 |
| **pandas**         | Conda-resolved          | Table processing and result parsing                 |
| **numpy**          | Conda-resolved          | Numerical processing                                |
| **scikit-learn**   | Conda-resolved          | Machine-learning models and evaluation              |
| **Biopython**      | Conda-resolved          | Sequence parsing and feature handling               |
| **PyYAML**         | Conda-resolved          | Reading workflow configuration files                |
| **requests**       | Conda-resolved          | UniProt / API data retrieval                        |
| **MAFFT**          | Conda-resolved          | Multiple sequence alignment                         |
| **FastTree**       | Conda-resolved          | Phylogenetic tree reconstruction                    |
| **Meeko**          | Conda-resolved          | Receptor preparation for docking                    |
| **AutoDock Vina**  | 1.2.5                   | Molecular docking                                   |
| **FoldX**          | FoldX4                  | Protein stability calculation and mutation screening |


### Prerequisites

- Conda / Miniconda / Mamba
- Git
- FoldX Suite (external, license required — see Configuration below)

### 1. Create and activate the Conda environment

```bash
conda env create -f environments/thesis.yml
conda activate xylanase-thesis
```

### 2. Configure the workflow

Before running, edit the main configuration file:

```bash
config/workflow_config.yaml
```

**Required manual configuration** (critical paths):

- FoldX executable path
- FoldX `rotabase.txt` path
- xylobiose ligand PDBQT path
- xylotriose ligand PDBQT path
- Output base directories
- Test-mode settings (optional)

### 3. Dry-run validation (recommended first step)

```bash
snakemake -s workflow/Snakefile -n -p --cores 4
```

This builds the directed acyclic graph (DAG) and checks for configuration or rule errors without executing any jobs.

### 4. Run the workflow

```bash
# Standard workflow run
snakemake -s workflow/Snakefile -p --cores 4

# Full workflow definition (all rules)
snakemake -s workflow/Snakefile.full -p --cores 4
```

> **Tip**: The `-p` flag prints the shell commands that will be executed. You can increase `--cores` according to your available resources.

## External Tools

### FoldX
FoldX is **not** included in the Conda environment. It must be obtained separately under the FoldX Suite license.

After installation, provide the paths to the FoldX executable and `rotabase.txt` in `config/workflow_config.yaml`.

### GROMACS
Molecular dynamics simulations were performed outside this Snakemake workflow and are not part of the active pipeline.

## Output Policy & Git Tracking

All generated outputs are intentionally **excluded from Git version control** via `.gitignore`.

This includes:
- Snakemake runtime folders (`.snakemake/`)
- Log files and runtime reports
- Raw downloads and fetched PDB structures
- FoldX working directories and repaired structures
- Docking receptor files and AutoDock Vina intermediate outputs
- Test-mode result folders
- Any large intermediate or final result files

Only source code, Snakefiles, modular rules, configuration templates, scripts, and documentation are tracked in this repository.
