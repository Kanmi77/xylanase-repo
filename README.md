# Xylanase Thermostability Assessment Repository

This repository contains the curated scripts, result tables, manifests and documentation for the thesis:

**In Silico Structural Thermostability Assessment of Bacterial and Fungal Thermostable Xylanases through Bioinformatics and Machine Learning Predictions**

The repository is organised as a cleaned thesis archive. It contains the scripts and final/supporting result files required to describe, reproduce and inspect the computational workflow. Large intermediate files such as molecular dynamics trajectories, docking intermediates, FoldX working folders and full model-building outputs are not included.

---

## Repository structure

| Folder | Description |
|---|---|
| `data/` | Curated and processed input datasets used in the thesis workflow |
| `docs/methods/` | Method notes and workflow documentation |
| `docs/manifests/` | File manifests with paths, checksums and descriptions |
| `results/` | Final and supporting result tables used for thesis reporting |
| `scripts/` | Thesis-aligned workflow scripts grouped by method stage |

---

## Script organisation

Each folder in `scripts/` corresponds to one major thesis workflow stage.  
Each stage contains one public runner script named `run_*.py` and the supporting scripts used within that stage.

| Stage folder | Public runner | Thesis workflow stage |
|---|---|---|
| `scripts/01_data_curation/` | `run_data_curation.py` | Data acquisition, UniProt/CAZy/BRENDA/PDB/RefSeq curation |
| `scripts/02_sequence_phylogeny/` | `run_sequence_phylogeny.py` | Sequence features, alignment, phylogeny and conservation |
| `scripts/03_structure_analysis/` | `run_structure_analysis.py` | PDB/MODELLER structural features and TM-align validation |
| `scripts/04_foldx/` | `run_foldx_analysis.py` | FoldX wild-type stability and mutation ΔΔG analysis |
| `scripts/05_mutation_screening/` | `run_mutation_screening.py` | Mutation panel preparation and parsing |
| `scripts/06_docking/` | `run_docking_analysis.py` | AutoDock Vina docking preparation, parsing and WT-mutant comparison |
| `scripts/07_machine_learning/` | `run_machine_learning.py` | Experimental-label machine learning and feature importance |
| `scripts/08_molecular_dynamics/` | `run_molecular_dynamics.py` | MD metric extraction and WT-mutant summary tables |
| `scripts/09_integration/` | `run_integration_and_ranking.py` | Integrated candidate scoring and final ranking |

---

## Main result folders

| Folder | Contents |
|---|---|
| `results/data_curation/` | Dataset preparation outputs |
| `results/sequence_phylogeny/` | Sequence, conservation and thermal-label outputs |
| `results/structure_analysis/` | Structural feature and TM-align outputs |
| `results/foldx/` | FoldX wild-type and mutation summaries |
| `results/mutation_screening/` | Mutation candidate outputs |
| `results/docking/` | Docking score and WT-mutant comparison outputs |
| `results/machine_learning/` | ML datasets, model reports and feature-importance outputs |
| `results/molecular_dynamics/` | MD summary metrics and comparison tables |
| `results/integration/` | Integrated candidate scoring and final ranking tables |

---

## Manifests

| Manifest | Description |
|---|---|
| `docs/manifests/scripts_manifest.tsv` | List of included scripts with workflow stage, role, path and SHA256 checksum |
| `docs/manifests/results_manifest.tsv` | List of included result files with path, size and SHA256 checksum |
| `docs/manifests/pdb_only_foldx_docking_manifest.txt` | PDB-only FoldX and docking supporting file list |
| `docs/manifests/pdb_reference_md_methodology_manifest.txt` | PDB reference and molecular dynamics methodology file list |

---

## Workflow summary

The thesis workflow contains the following computational stages:

1. Data acquisition and curation of GH10/GH11 bacterial and fungal xylanases.
2. Sequence feature extraction, multiple sequence alignment and phylogenetic analysis.
3. PDB and MODELLER structural-feature extraction with TM-align validation.
4. FoldX wild-type stability and mutation ΔΔG analysis.
5. Mutation candidate preparation and screening.
6. AutoDock Vina docking of wild-type and mutant structures with xylo-oligosaccharide ligands.
7. Machine-learning prediction using experimental thermal labels and sequence/structure features.
8. Molecular dynamics validation of selected wild-type and mutant systems.
9. Integrated candidate scoring and final thermostability/function-retention ranking.

---

## Excluded files

The following files were intentionally excluded from the public repository because they are large, intermediate, software-generated or not required for thesis inspection:

- GROMACS trajectories and binary simulation files: `.xtc`, `.trr`, `.tpr`, `.edr`, `.cpt`, `.gro`
- AutoDock/Vina intermediate receptor and ligand files: `.pdbqt`, docking logs and temporary folders
- FoldX executables, databases and working directories
- MODELLER temporary output folders and large model-generation intermediates
- Raw exploratory branches not used in the final thesis reporting

---

## Reproducibility notes

The repository provides a thesis-level reproducibility archive rather than a single-command full rerun of every computationally expensive step. Some stages require external tools and local installations, including:

- Python
- pandas, NumPy and scikit-learn
- Biopython
- MAFFT
- FastTree
- MODELLER
- FoldX4
- AutoDock Vina
- GROMACS
- TM-align
- FreeSASA

The included manifests provide file paths and SHA256 checksums for traceability of the archived scripts and result tables.

---

## Author

Kamaldeen Olasunkanmi Bada  
Master's Thesis, Life Science Informatics  
Deggendorf Institute of Technology
