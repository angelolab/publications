# Code for Ranek et al., QUICHE reveals structural definitions of anti-tumor responses in triple negative breast cancer

## Introduction
QUICHE is a statistical differential abundance testing method that can be used to discover cellular niches differentially enriched in spatial regions, longitudinal samples, or clinical patient groups. If you'd like to perform spatial enrichment analysis using QUICHE, please see the associated repo: https://github.com/jranek/quiche. Alternatively, if you'd like to reproduce the analysis from the paper, please see below.

## Data access
You can download all of the preprocessed MIBI-TOF datasets (`.h5ad` files) from the [Zenodo](https://zenodo.org/records/14290163) repository. Imaging data and cell segmentation masks can be found in the [BioStudies](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1507) repository.

## Installation
You can clone the git repository by, 
```
git clone https://github.com/angelolab/publications.git
```
Then change the working directory as, 
```
cd publications/2024-Ranek_etal_QUICHE
```

Given that there are a number of python packages for spatial enrichment evaluation, we recommend that you create a conda environment using the provided yml file.

```
conda env create -f venv_quiche_benchmark.yml
```

Once the environment is created, you can activate it by,
```
conda activate venv_quiche_benchmark
```

If you'd like to reproduce some of the analyses from the paper, you'll also need to install the necessary R packages.

```R
## QUICHE analysis
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#edger v3.40.2
BiocManager::install("edgeR")
#statmod v1.5.0
install.packages('statmod')

## Logistic regression analysis
#glmnet v4.1-8
install.packages('glmnet')

## Cox proportional hazards regression analysis
#survival v3.5-8
install.packages('survival')

## Simulation analysis
#splatter v1.18.2
BiocManager::install("splatter")
```

## Description of files
* `quiche` contains scripts for performing differential spatial enrichment analysis with QUICHE, benchmarking, and plotting. 
* `scripts_TNBC_preprocessing`: contains preprocessing scripts for Spain, Stanford, and NeoTRIP cohorts, including cell segmentation, cell phenotyping, tumor region identification, and collagen1 fiber segmentation.
* `scripts_benchmarking`: contains scripts for benchmarking spatial enrichment methods, including in silico data generation and evaluation.
* `abundance_prediction.R`: performs logistic regression analysis using cell type abundances.
* `utils.R`: contains helper functions for logistic regression analysis.
* `cox_model.R`: performs Cox proportional hazards regression analysis.
* `figure01-05` reproduce the main figures in the paper.
* `supplementary_figure01-24` reproduce the supplementary figures in the paper.
* `supplementary_plot_helpers.py`: contains helper functions for figure generation.
