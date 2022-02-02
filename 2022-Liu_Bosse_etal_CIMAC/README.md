# Code for Liu and Bosse et al., Reproducible, high-dimensional imaging in archival human tissue by Multiplexed Ion Beam Imaging by Time-of-Flight (MIBI-TOF)

Description of files:
1. mpi_ppp_ihc_regression.ipynb: Calculates normalization coefficients and normalizes mean pixel intensity (MPI) and percent positive pixel (PPP), performs regression for MPI/PPP, performs comparison with IHC data. Run this notebook first.
- Generates plots in Figure 3C-D, Figure 5, Supplementary Figure 4B-C, Supplementary Figure 9
2. mpi_ppp_scatterplots.R: Creates MPI/PPP scatterplots and plots CV of each core
- Generates plots in Figure 3B, Supplementary Figure 4A, Supplementary Figure 5
3. cellClustering.R: Clusters cells using mean expression of each cell (after segmentation done in Mesmer)
- Generates heatmap in Figure 4A
4. cellCluster_correlations.R: Calculates Spearman correlation between serial sections of same core using cell phenotype frequencies
- Generates plots in Figure 4B, Supplementary Figure 8
5. cell_tsne.R: Creates tSNE of all cells
- Generates plots in Supplementary Figure 7C
6. plot_norm_coef.R: Creates plots analyzing normalization coefficients
- Generates plots in Supplementary Figure 2

All data tables needed to run these scripts available on Zenodo at https://doi.org/10.5281/zenodo.5945388 

