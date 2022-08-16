# Code for Liu et al., Robust phenotyping of highly multiplexed tissue imaging data using pixel-level clustering

**User-friendly Jupyter notebooks for running Pixie are available at github.com/angelolab/ark-analysis. Please see that repository to use Pixie for your own data.** The scripts here were used to generate the figures in the paper.

Description of files:
1. 0_pixelClustering_functions.R: Functions used for pixel clustering and cell clustering
2. 1_imageToMatrix.py: Converts single-channel TIFs into pixel x markers matrix, necessary before all pixel clustering
3. 2_pixelClustering.R: Performs pixel clustering
    - Figure 2b-c,3,4; Supp. Figure 2-12
4. 3_pixelCustering_heatmaps.R: Creates heatmaps of pixel cluster expression
    - Figure 2b, Supp. Figure 2-12
5. 4_save_pixelClusterOverlay_TIFs.py: Saves pixel phenotype maps
    - Figure 2c, Supp. Figure 2-12
6. 5_cluster_consistency_score_pixel.R: Calculates cluster consistency score for 5 replicates of a dataset (pixel)
    - Figure 2d-e, Supp. Figure 4-10
7. 6_cluster_consistency_score_pixel.py: Completes calculation of cluster consistency score and creates colored overlays
    - Figure 2d-e, Supp. Figure 4-10
8. 7_compare_cluster_consistency_score.R: Compares cluster consistency score between different parameter choices
    - Figure 2f
9. 8_single_cell_cytof.R: Uses FlowSOM to cluster a single cell CyTOF dataset (Hartmann, et al. 2020)
    - Supp. Figure 4c,d
10. 9_single_cell_scrnaseq.R: Uses FlowSOM to cluster a single-cell RNA-sequencing dataset (Seurat)
    - Supp. Figure 4c,e
11. 10_get_pixels_outside_cells.py: Counts percentage of pixels outside of cellular masks
    - Figure 3a-c
12. 11_compare_pixels_outside_cells.R: Compares of percentage of pixels outside of cellular masks across datasets
    - Figure 3d
13. 12_pixel_cluster_corelations.R: Calculates corrleation of pixel cluster frequency between replicate serial sections
    - Supp. Figure 12
14. 13_cellClustering_pixelComposition.R: Performs cell clustering using pixel cluster composition
    - Figure 5c, Supp. Figure 14
15. 14_cellClustering_integratedExpression.R: Performs cell clustering using integrated expression (as a comparison)
    - Figure 5b
16. 15_save_cellClusterOverlay_TIFs.py: Saves cell phenotype maps
    - Figure 5e, Supp. Figure 13-14
17. 16_cluster_consistency_score_cell.R: Calculates cluster consistency score for 5 replicates of a datset (cell)
    - Supp. Figure 13
18. 17_cluster_consistency_score_cell.py: Completes calculation of cluster consistency score and creates colored overlays
    - Supp. Figure 13
19. 18_compare_silhouette_scores.R: Computes Silhouette scores for cell clustering
    - Figure 5d
20. 19_som_leiden_phenograph_test.R: Performs time test for different clustering algorithms
    - Supp. Figure 15


