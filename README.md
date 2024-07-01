This repository contains all the code for processing and analysis of benzonase cross-linking-and-immunoprecipitation (bCLIP) for Integrator complex proteins in mouse embryonic stem cells (mESC) and RNA-sequencing (RNA-seq) data for siRNA-knockdowns of Integrator components in mESC and human HEK293 cell lines, generated in the lab of prof. Dr. Stefanie Jonas by Moes Murielle.
In addition, several public datasets were downloaded and analyzed alongside:

- enhanced cross-linking-and-immunoprecipitation (eCLIP) data for

Briefly, the computations were done in the following steps:

**1)** Running jupyter notebook "bCLIP.ipynb".
It contains the code to download necesary public data, organize the input files into folders, create custom annotation .gtf files (union of GENCODE + RNAcentral), and prepare the configuration for the snakemake workflow.

**2)**  
