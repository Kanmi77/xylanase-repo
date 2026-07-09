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
