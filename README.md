This folder contains raw analysis scripts for the following publication:

Visser K, Thomsen L.B, Miedema A, de Koning M.E, van Goor H, van der Naalt J, Nossent A.Y, Jacobs B, van der Horn HJ. Uncovering the utility of circulating microRNAs in acute mild traumatic brain injury, Brain Communications, 2026.

*Please note that all notebooks are published without change with absolute references to accompanying files including helper functions and constants. 
These separated files are included. Please adjust these environment variables where necessary. 

- mirna_analysis_main.ipynb is a notebook containing the machine learning pipeline analysis. It expects the miRNA to have been preprocessed, this notebook contains all the code necessary for recreating graphs which are later processed in Inkscape.
- mirna_qc.ipynb used metrics from srnatoolbox to calculate metrics such as readcounts
- mirna_validation.ipynb contains code for calculating demographics and injury variables as well as various supplementary analyses.
- The files aimtbi_** contains project wide constants which are shared across multiple studies. Therefore certain functions are not implemented in the above notebooks
- enrichment_analysis.Rmd contains an R script to perform enrichment analysis on selected microRNA
- Further files including aim_micro_rna_job.sh, inside_container_v1.sh, populate_hpc.py are used to run sRNAbench inside a dockerfile and carry out differential expression analysis. Refer to sRNABench documentation for a detail guide. The bash files were run inside the University of Groningen High Performance Cluster Habrok. 

Raw data can be obtained upon reasonable request to the corresponding author of the above study. 

September 9, 2026. Dr. K, Visser. University of Groningen, University Medical Center Groningen, Groningen, the Netherlands.
