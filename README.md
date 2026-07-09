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
| **OS** | Linux (tested on Ubuntu 22.04 LTS) |
| **RAM** | Minimum 32 GB (64 GB recommended for MD) |
| **Storage** | 50 GB free disk space |
| **Conda** | Miniconda or Anaconda |
| **Git** | Installed |
| **GPU** | Optional but recommended for MD simulations |

---

## Installation

### 1. Clone the repository

```bash
cd ~
git clone https://github.com/your-username/xylanase-repo
cd xylanase-repo

2. Create conda environments
bash
cd environments
conda env create -f xylanase.yml
conda env create -f foldx.yml
conda env create -f gromacs.yml
cd ..
3. Install external tools
These tools should be installed manually. The pipeline expects them at the locations below with the exact environment names.

Tool	Installation	Expected Location	Environment Name
FoldX 4.0	FoldX download	~/tools/foldx/	foldx
GROMACS 2025.4	GROMACS installation	~/tools/gromacs/	gromacs
AutoDock Vina 1.2.5	pip install vina	~/tools/vina/	vina
MODELLER 10.8	MODELLER installation	~/tools/modeller/	modeller
Important notes:

For FoldX, the rotabase.txt file must be placed in the active directory before running RepairPDB and Stability runs.

For GROMACS, ensure the force field (AMBER99SB-ILDN) is available.

All scripts use conda activate <env_name> — environments must have exactly these names.

Adjust paths in scripts if tools are installed to different locations.

Input Preparation
GH Families
The pipeline is designed for GH10 and GH11 xylanases. Create a GH_families.txt file inside the data/ directory:

bash
cd ~/xylanase-repo/data
touch GH_families.txt
echo "GH10" >> GH_families.txt
echo "GH11" >> GH_families.txt
UniProt Sequence Retrieval
The pipeline will automatically fetch sequences from UniProt for the specified GH families. Alternatively, you can provide your own FASTA files in data/proteomes/:

bash
mkdir -p data/proteomes/
# Place your .faa or .fasta files here
Reference Structures for TM-Align Validation
Reference structures (experimental PDB files) for structural validation should be placed in subdirectories by GH family:

text
data/reference_structures/GH10/1O8S.pdb
data/reference_structures/GH10/1VBR.pdb
data/reference_structures/GH11/1VBU.pdb
data/reference_structures/GH11/3NIY.pdb
Thermal Labels (BRENDA Data)
Experimental temperature data from BRENDA should be placed in:

text
data/thermal_labels/brenda_temperature.csv
Ligand Files for Docking
Xylobiose and xylotriose ligand files (in PDBQT or SDF format) should be placed in:

text
results/docking/ligands/xylobiose.pdbqt
results/docking/ligands/xylotriose.pdbqt
Running the Pipeline
Activate the base environment
bash
conda activate xylanase
Dry run (test the workflow)
bash
snakemake --cores 8 --dry-run
Full run
bash
snakemake --cores 8
Run specific targets
bash
# Run only sequence analysis
snakemake results/tables/physicochemical_summary.csv

# Run only FoldX mutation screening
snakemake results/foldx/modeller_mutations.csv

# Run only machine learning
snakemake results/tables/ml_performance.csv
Output Structure
Directory	Contents
results/tables/	All CSV and TSV summary tables (physicochemical, phylogenetic, conservation, FoldX, ML, MD)
results/figures/	Visualisations (bar charts, scatter plots, phylogenetic trees, heatmaps)
results/foldx/	FoldX repair logs, stability outputs, mutation screening results
results/docking/	Docking scores, binding retention assessments
results/md/	MD trajectories, RMSD/RMSF/Rg/H-bond analyses, DSSP outputs
output_files/	Final integrated candidate rankings and source-stratified summaries
Customisation
Modify Filtering Thresholds
Edit scripts/filtering/filter_candidates.py to change:

Parameter	Default	Description
--foldx_threshold	0.0	Maximum FoldX ΔΔG for stabilising mutations
--docking_threshold	0.0	Maximum mean docking change for retention
--rmsd_threshold	1.0	Maximum RMSD for MD validation
--tm_score_threshold	0.75	Minimum TM-score for homology model validation
Modify Docking Parameters
Edit scripts/docking/run_docking.sh to change:

bash
EXHAUST=8              # Exhaustiveness (default: 8, increase for accuracy)
ENERGY_RANGE=4         # Energy range for output poses
NUM_MODES=10           # Number of binding modes to output
Modify MD Parameters
Edit scripts/md/run_gromacs.sh to change:

bash
SIMULATION_TIME=10000  # Simulation time in ps (10 ns default)
TEMPERATURE=373        # Simulation temperature in K
PRESSURE=1.0           # Pressure in bar
Reproducibility
All scripts, workflow definitions, and intermediate data generated in this study are systematically archived to ensure full computational reproducibility. The complete pipeline, including parameter specifications and processing logs, has been preserved.

To reproduce the analysis:

bash
# 1. Clone the repository
git clone https://github.com/your-username/xylanase-repo

# 2. Set up environments
cd xylanase-repo/environments
conda env create -f xylanase.yml
conda env create -f foldx.yml
conda env create -f gromacs.yml

# 3. Run the workflow
cd ..
conda activate xylanase
snakemake --cores 8

# 4. All results will be in results/ and output_files/

Citation
If you use this pipeline, please cite:


Bada, K. (2026). In Silico Structural Thermostability Assessment of 
Bacterial and Fungal Thermostable Xylanases through Bioinformatics 
and Machine Learning Predictions. Master's Thesis, Technische 
Hochschule Deggendorf.
