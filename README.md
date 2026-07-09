# xylanase-repo
Computational workflow for in silico thermostability prediction of GH10 and GH11 xylanases. Integrates sequence analysis, structural modelling, FoldX mutagenesis, docking, machine learning, and molecular dynamics validation. Includes Snakemake workflow, Python scripts, and result visualisations.
# Xylanase Thermostability Prediction Pipeline

Pipeline for screening and prioritising thermostabilising mutations in bacterial and fungal GH10 and GH11 xylanases through integrated bioinformatics and machine learning.

This pipeline was developed as part of a Master's thesis:

> Bada, K. (2026). *In Silico Structural Thermostability Assessment of Bacterial and Fungal Thermostable Xylanases through Bioinformatics and Machine Learning Predictions*. Technische Hochschule Deggendorf.

---

## Overview

This Snakemake pipeline processes GH10 and GH11 xylanase sequences to:

- Curate sequence and structural data from UniProt, CAZy, and BRENDA
- Perform physicochemical characterisation and phylogenetic analysis
- Identify conserved motifs and thermal-contrast positions
- Build homology models using MODELLER
- Predict stability effects of mutations using FoldX
- Assess substrate binding retention via molecular docking (AutoDock Vina)
- Train machine learning models for stability prediction
- Validate top candidates through molecular dynamics (GROMACS)
- Integrate all evidence into a source-stratified candidate ranking


---


## System Requirements

| Requirement | Specification |
| :--- | :--- |
| **OS** | Linux (Ubuntu 22.04 LTS) |
| **RAM** | Minimum 32 GB (64 GB recommended for MD) |
| **Storage** | 50 GB free disk space |
| **Conda** | Miniconda or Anaconda |
| **GPU** | Optional but recommended for MD simulations |

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/your-username/xylanase-repo
cd xylanase-repo

# 2. Create conda environments
cd environments
conda env create -f xylanase.yml
conda env create -f foldx.yml
conda env create -f gromacs.yml
cd ..

# 3. Install external tools manually
# FoldX 4.0: https://foldxsuite.crg.eu/
# GROMACS 2025.4: https://manual.gromacs.org/
# MODELLER 10.8: https://salilab.org/modeller/
# AutoDock Vina 1.2.5: pip install vina

Usage
Input Preparation
1. Proteomes to be screened
The proteomes to be screened by the pipeline should be in .faa or .fasta format.

The proteomes directory should be located in data/proteomes/:

bash
cd ~/xylanase-repo
mkdir -p data/proteomes/
cd data/proteomes/

# Place your proteomic data here

cd ~/xylanase-repo
2. GH families
The pipeline is designed for GH10 and GH11 xylanases. A GH_families.txt file should be created with the desired GH families (one per line). Subfamily level can be used.

## Input Files

Before running the pipeline, prepare the required input files.

### GH Family List

Create a list of glycoside hydrolase (GH) families to analyse.

```bash
cd data
touch GH_families.txt

echo "GH10" >> GH_families.txt
echo "GH11" >> GH_families.txt
```

**Example**

```text
GH10
GH11
```

---

### Reference Structures

Reference crystal structures are required for TM-align validation of predicted protein structures.

Store PDB files under:

```text
data/
└── reference_structures/
    ├── GH10/
    │   ├── 1O8S.pdb
    │   └── 1VBR.pdb
    └── GH11/
        ├── 1VBU.pdb
        └── 3NIY.pdb
```

Create the directories:

```bash
mkdir -p data/reference_structures/GH10
mkdir -p data/reference_structures/GH11
```

> **Note**
> All reference structures must be in **PDB (.pdb)** format.

---

### Docking Ligands

By default, the pipeline docks **xylobiose** and **xylotriose**.

Create the ligand directory:

```bash
mkdir -p results/docking/ligands
```

Place ligand files such as:

```text
results/
└── docking/
    └── ligands/
        ├── xylobiose.pdbqt
        └── xylotriose.pdbqt
```

> **Note**
> To use different ligands, edit:
>
> `scripts/docking/run_docking.sh`

---

### Thermal Labels (Optional)

Experimental temperature labels from BRENDA can be included for machine learning.

```bash
mkdir -p data/thermal_labels
```

Expected file:

```text
data/thermal_labels/brenda_temperature.csv
```

---

# Running the Pipeline

Activate the Conda environment.

```bash
conda activate xylanase
```

### Dry Run

Verify the workflow before execution.

```bash
snakemake --cores 8 --dry-run
```

### Run Entire Pipeline

```bash
snakemake --cores 8
```

### Run Individual Modules

| Module | Command |
|---------|---------|
| Sequence analysis | `snakemake results/tables/physicochemical_summary.csv` |
| FoldX mutation screening | `snakemake results/tables/foldx_mutations.csv` |
| Machine learning | `snakemake results/tables/ml_performance.csv` |
| Final candidate ranking | `snakemake results/tables/final_candidates.csv` |

---

# Output

The pipeline generates the following outputs.

| Directory | Description |
|------------|-------------|
| `output_files/` | Final integrated results and candidate rankings |
| `results/tables/` | Summary tables generated throughout the workflow |
| `results/figures/` | Publication-quality figures |
| `results/foldx/` | FoldX repair and mutation results |
| `results/docking/` | Docking scores, receptors and ligand files |
| `results/md/` | Molecular dynamics trajectories and analyses |

> **Note**
> Docking results can be visualised using **PyMOL**, **ChimeraX**, or other molecular visualization software.

---

# Customisation

## Candidate Filtering

Edit:

```text
scripts/filtering/filter_candidates.py
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `foldx_threshold` | `0.0` | Maximum ΔΔG allowed |
| `docking_threshold` | `0.0` | Maximum docking score change |
| `rmsd_threshold` | `1.0` | RMSD cutoff |
| `tm_score_threshold` | `0.75` | Minimum TM-score |
| `residue_length_threshold` | `800` | Maximum protein length |

---

## Docking Parameters

Edit:

```text
scripts/docking/run_docking.sh
```

| Parameter | Default |
|-----------|---------|
| `EXHAUST` | `32` |
| `ENERGY_RANGE` | `4` |
| `NUM_MODES` | `10` |

---

## Molecular Dynamics Parameters

Edit:

```text
scripts/md/run_gromacs.sh
```

| Parameter | Default |
|-----------|---------|
| `SIMULATION_TIME` | `10000 ps` |
| `TEMPERATURE` | `373 K` |
| `PRESSURE` | `1.0 bar` |

---

# Reproducibility

Clone the repository.

```bash
git clone https://github.com/your-username/xylanase-repo.git
cd xylanase-repo
```

Create the environments.

```bash
cd environments

conda env create -f xylanase.yml
conda env create -f foldx.yml
conda env create -f gromacs.yml

cd ..
```

Prepare the required input files as described above.

Run the workflow.

```bash
conda activate xylanase
snakemake --cores 8
```

---

