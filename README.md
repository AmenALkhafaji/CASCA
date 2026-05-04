# CASCA
New versión of ASCA for consensus diffrent method togather
# C-ASCA for Metagenomics Data

This repository contains MATLAB code for applying **Consensus ANOVA--Simultaneous Component Analysis (C-ASCA)** to metagenomics count data.

C-ASCA integrates multiple preprocessing/normalization outputs from the same microbiome dataset into a multiblock ASCA model. The framework separates normalization-dependent disagreement from shared biological consensus and exports a consensus representation for downstream analysis.

## Requirements

This code requires:

- MATLAB
- MEDA Toolbox version 1.12

Download MEDA-Toolbox-1.12 from:

https://github.com/josecamachop/MEDA-Toolbox

After downloading MEDA, add it to your MATLAB path:

```matlab
addpath(genpath('path_to/MEDA-Toolbox-1.12'));
savepath;

Simulation analysis

To run the simulation analysis:

Download the required simulation input files.
Place the simulation files in the same working directory as CASCA.m, unless you modify the file paths inside the scripts.
Run the main C-ASCA simulation script:
CASCA

To generate the simulation visualization panels, run:

CASCA_panel

To compute confusion-matrix results, including true positives, false positives, false negatives, precision, recall, and F1 score, run:

Confmat
Real-world Schubert CDI dataset analysis

To run the real-data analysis:

Download the Schubert CDI dataset files.
Place the Schubert files in the same working directory as schubert.m, unless you modify the file paths inside the scripts.
Run:
schubert

To compute overlap and method-comparison results for the Schubert dataset, run:

schubert_overlap
Important notes

The MEDA Toolbox must be added to the MATLAB path before running any C-ASCA script.

If MATLAB reports missing functions, check that both MEDA-Toolbox-1.12 and this repository have been added recursively to the MATLAB path.

If MATLAB reports missing input files, check that the dataset files are located in the directory expected by the scripts.

If the scripts contain local file paths, update those paths before running the analysis.

If you use this code, please cite the associated manuscript:

Al Khafaji A, Gómez-Llorente C, Camacho J. Consensus ANOVA--Simultaneous Component Analysis for the Normalization of Metagenomics Data.

Update this citation with the final journal reference once the manuscript is published.

Contact

For questions about the code, please contact:

Amen Al Khafaji
amen.a.khabeer@uotechnology.edu.iq

