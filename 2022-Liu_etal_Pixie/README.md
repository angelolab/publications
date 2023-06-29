# Code for Liu et al., Robust phenotyping of highly multiplexed tissue imaging data using pixel-level clustering

### User-friendly Jupyter notebooks for running Pixie are available at [github.com/angelolab/pixie](https://github.com/angelolab/pixie). Please see that repository to use Pixie for your own data.

The scripts here were used to generate the figures in the paper. Description of files:
1. 0_pixelClustering_functions.R: Functions for pixel clustering and cell clustering
2. 1_imageToMatrix.py: Converts single-channel TIFs into pixel x markers matrix, necessary before pixel clustering
3. 2_pixelClustering.R: Performs pixel clustering
    - Figure 2b-c,3,4; Supp. Figure 2-15
4. 3_pixelCustering_heatmaps.R: Creates heatmaps of pixel cluster expression
    - Figure 2b, Supp. Figure 2-15
5. 4_save_pixelClusterOverlay_TIFs.py: Saves pixel phenotype maps
    - Figure 2c, Supp. Figure 2-15
6. 5_cluster_consistency_score_pixel.R: Calculates cluster consistency score for 5 replicates of a dataset (pixel)
    - Figure 2d-e, Supp. Figure 4-6, 8-12
7. 6_cluster_consistency_score_pixel.py: Completes calculation of cluster consistency score and creates colored overlays
    - Figure 2d-e, Supp. Figure 4-6, 8-12
8. 7_compare_cluster_consistency_score.R: Compares cluster consistency score between different parameter choices
    - Figure 2f
9. 8_single_cell_cytof.R: Uses Leiden and FlowSOM to cluster a single cell CyTOF dataset (Hartmann, et al. 2020)
    - Supp. Figure 4c,d
10. 9_single_cell_scrnaseq.R: Uses Leiden and FlowSOM to cluster a single-cell RNA-sequencing dataset (Seurat)
    - Supp. Figure 4c,e
11. 10_get_pixels_outside_cells.py: Counts percentage of pixels outside of cellular masks
    - Figure 3a-c
12. 11_compare_pixels_outside_cells.R: Compares percentage of pixels outside of cellular masks across datasets
    - Figure 3d
13. 12_pixel_cluster_corelations.R: Calculates corrleation of pixel cluster frequency between replicate serial sections
    - Supp. Figure 14
14. 13_cellClustering_pixelComposition.R: Performs cell clustering using pixel cluster composition
    - Figure 5c, Supp. Figure 16-19, 21
15. 14_cellClustering_integratedExpression.R: Performs cell clustering using integrated expression (as a comparison)
    - Figure 5b, Supp. Figure 16-19
16. 15_save_cellClusterOverlay_TIFs.py: Saves cell phenotype maps
    - Figure 5e, Supp. Figure 16-19, 21
17. 16_cluster_consistency_score_cell.R: Calculates cluster consistency score for 5 replicates of a datset (cell)
    - Supp. Figure 16-17
18. 17_cluster_consistency_score_cell.py: Completes calculation of cluster consistency score and creates colored overlays
    - Supp. Figure 16-17
19. 18_compare_silhouette_scores.R: Computes Silhouette scores for cell clustering
    - Figure 5d
20. 19_imageToMatrix_preprocessed.py: Performs pre-processing steps (pixel normalization and channel normalization) on individual images
    - Supp. Figure 17
21. 20_otsu_thresholding.py: Performs Otsu thresholding on single-channel images
    - Supp. Figure 20
22. 21_otsu_thresholding_count_combinations.R: Count the number of different phenotypes resulting from Otsu thresholding
    - Supp. Figure 20
23. 22_pixie_pixels_runtime_test.py: Performs runtime test for pixel clustering in Pixie
    - Supp. Figure 22a
24. 23_pixie_cells_runtime_test.py: Performs runtime test for cell clustering in Pixie
    - Supp. Figure 22b
25. 24_som_leiden_phenograph_test.R: Performs runtime test for different clustering algorithms
    - Supp. Figure 22c
