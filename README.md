# C-ASCA for Metagenomics DataThis repository contains MATLAB code for applying **Consensus ANOVA--Simultaneous Component Analysis (C-ASCA)** to metagenomics count data.C-ASCA extends ASCA to integrate multiple preprocessing or normalization outputs from the same microbiome dataset into a multiblock model. The framework separates normalization-dependent disagreement from shared biological consensus and exports a consensus representation for downstream analysis.## RequirementsThis code requires:- MATLAB- MEDA Toolbox version 1.12Download MEDA-Toolbox-1.12 from:https://github.com/josecamachop/MEDA-ToolboxAfter downloading MEDA, add it to your MATLAB path:```matlabaddpath(genpath('path_to/MEDA-Toolbox-1.12'));savepath;
Replace path_to/MEDA-Toolbox-1.12 with the actual path on your computer.
Then add this C-ASCA repository to the MATLAB path:
addpath(genpath('path_to/CASCA_repository'));savepath;
Replace path_to/CASCA_repository with the actual path of this repository.
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


Run the main real-data analysis script:


schubert
To compute overlap and method-comparison results for the Schubert dataset, run:
schubert_overlap
Important notes
The MEDA Toolbox must be added to the MATLAB path before running any C-ASCA script.
If MATLAB reports missing functions, check that both MEDA-Toolbox-1.12 and this repository have been added recursively to the MATLAB path.
If MATLAB reports missing input files, check that the dataset files are located in the directory expected by the scripts.
If the scripts contain local file paths, update those paths before running the analysis.
Citation
If you use this code, please cite the associated manuscript:
Al Khafaji A, Gómez-Llorente C, Camacho J. Consensus ANOVA--Simultaneous Component Analysis for the Normalization of Metagenomics Data.
Update this citation with the final journal reference once the manuscript is published.
Contact
For questions about the code, please contact:
Amen Al Khafaji
amen.a.khabeer@uotechnology.edu.iq
License
This repository is released under the MIT License. See the LICENSE file for details.
