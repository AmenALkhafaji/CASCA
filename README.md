# C-ASCA for Metagenomics Data

This repository contains MATLAB code for applying **Consensus ANOVA–Simultaneous Component Analysis (C-ASCA)** to metagenomics count data using multiple normalization outputs from the same microbiome dataset.

C-ASCA extends the classical ASCA framework by integrating several preprocessing or normalization methods into a multiblock model. The method separates normalization-dependent disagreement from shared biological consensus and exports a consensus-normalized representation for downstream analysis.

---

## Overview

Normalization choice is a major source of variability in microbiome data analysis. Different preprocessing methods can substantially alter multivariate structure, feature rankings, and biological interpretation.

This repository provides a consensus framework that:

- Integrates multiple normalization outputs jointly  
- Quantifies agreement and disagreement across methods  
- Extracts shared biological signal  
- Produces a consensus data representation  
- Supports simulation and real-world benchmarking

Applications included:

- Simulation benchmarking with known ground truth  
- Real-world *Schubert CDI* dataset analysis

---

## Contact

**Amen Al Khafaji**  
Email: amen.a.khabeer@uotechnology.edu.iq

Last document update: **16/05/2026**

---

## Requirements

The following software is required:

- MATLAB
- MEDA Toolbox version 1.12

MEDA Toolbox repository:

https://github.com/josecamachop/MEDA-Toolbox

---

## Installation

### 1. Download MEDA Toolbox

Download MEDA Toolbox version 1.12 from:

https://github.com/josecamachop/MEDA-Toolbox

### 2. Add MEDA Toolbox to MATLAB Path

addpath(genpath('path_to/MEDA-Toolbox-1.12'));
savepath;


3. Add This Repository to MATLAB Path
addpath(genpath('path_to/CASCA_repository'));
savepath;

Replace path_to/CASCA_repository with the location of this repository.

---

### Repository Structure
.
├── README.md
├── LICENSE
├── CASCA.m
├── CASCA_panel.m
├── Confmat.m
├── schubert.m
├── schubert_overlap.m
├── Data/
└── Figures/

---

## Folder Description
| Folder / File      | Description                               |
| ------------------ | ----------------------------------------- |
| CASCA.m            | Main simulation analysis script           |
| CASCA_panel.m      | Generates simulation visualization panels |
| Confmat.m          | Computes confusion-matrix metrics         |
| schubert.m         | Main real-data analysis script            |
| schubert_overlap.m | Overlap and comparison analysis           |
| Data/              | Input simulation and real-world datasets  |
| Figures/           | Exported output figures                   |


---

## General Workflow

Load raw microbiome count data
        ↓
Generate multiple normalization outputs
        ↓
Build multiblock matrix
        ↓
Run C-ASCA decomposition
        ↓
Separate consensus and disagreement structure
        ↓
Interpret scores, loadings, and contributions
        ↓
Export consensus-normalized matrix

---

## Usage

Option 1: Simulation Analysis

Use this workflow for controlled benchmarking with known ground truth.

Step 1: Prepare Simulation Files

Download or generate the required simulation datasets.

Place all files in the same working directory as:

CASCA.m
Step 2: Run Main Simulation Analysis
CASCA
Step 3: Generate Visualization Panels
CASCA_panel
Step 4: Compute Classification Metrics
Confmat

---

## Outputs may include:

True positives
False positives
False negatives
Precision
Recall
F1 score

---


### Output

Typical outputs include:

Consensus score plots
Loading plots
Method disagreement structure
Simulation recovery metrics
Overlap statistics
Exported consensus-normalized matrix
Supported Normalization Methods

---

## The framework is designed for integrating multiple preprocessing outputs, including:

Raw counts
Total-sum scaling (TSS)
Rarefaction
Centered log-ratio (CLR)
Other normalization strategies supplied by the user
Important Notes
MEDA Toolbox must be added to the MATLAB path before running any script.
If MATLAB reports missing functions, confirm recursive path installation.
If input files are not found, verify dataset locations.
If scripts contain local hard-coded paths, update them before execution.

---

### Citation

If you use this repository, please cite:

Al Khafaji A, Gómez-Llorente C, Camacho J. Consensus ANOVA–Simultaneous Component Analysis for the Normalization of Metagenomics Data. Manuscript in preparation.

Update this citation once the final journal reference becomes available.




--- 

## License

This repository is released under the MIT License.

See the LICENSE file for full details.





