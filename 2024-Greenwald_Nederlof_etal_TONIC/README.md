# Code for Greenwald, Nederlof, et al., Temporal and spatial composition of the tumor microenvironment predicts response to immune checkpoint inhibition

This repo contains scripts for analyzing the TONIC MIBI data. Below is a description of how to navigate the TONIC datasets, with specific information regarding the data file formats, as well as the scripts used to generate the data.

## Table of Contents
- [Scripts](#scripts)
- [Directory Structure](#directory-structure)
- [Data Structures](#data-structures)
- [Analysis Files](#analysis-files)
- [Output Files](#output-files)



## Scripts 
### image_processing
[4_compensate_image_data_working.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/image_processing/4_compensate_image_data_working.ipynb): This notebook will guide you through the Rosetta algorithm, which is used to remove background and noise from image data prior to analysis. It is based on the compensation matrix approach that has been used for correcting flow-cytometry data. The Rosetta matrix contains rows for each of the sources of noise, and columns for each of the output channels. Each entry in the matrix represents the proportional contamination from a given noise channel to a given output channel.

[4b_normalize_image_data_working.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/image_processing/4b_normalize_image_data_working.ipynb): This notebook will walk you through the process of normalizing the image data. Changes in detector sensitivity over a run can result in different image intensities, even when there are no actual difference in biological signal. To correct for this, we use the median pulse height (MPH) to measure detector sensitivity. We then combine this estimate of sensitivity with an instrument tuning curve to determine the normalization coefficient for each FOV.

[5_rename_and_reorganize_working.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/image_processing/5_rename_and_reorganize_working.ipynb): This notebook wa used to organize the image data following processing so that it is ready to be analyzed.


### single_cell_analysis

[Segment_Image_Data_working.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/Segment_Image_Data_working.ipynb): This notebook contains the process of using [Mesmer](https://www.nature.com/articles/s41587-021-01094-0) to segment your image data. This includes selecting the appropriate channel(s) for segmentation, running your data through the network, and then extracting single-cell statistics from the resulting segmentation mask.

[TONIC_pixel_clustering.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/TONIC_pixel_clustering.ipynb): The first step in the [Pixie](https://www.nature.com/articles/s41467-023-40068-5) pipeline is to run the pixel clustering notebook. The notebook walks you through the process of generating pixel clusters for the data, and lets you specify what markers to use for the clustering, train a model, use it to classify the entire dataset, and generate pixel cluster overlays.

[TONIC_cell_clustering.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/TONIC_cell_clustering.ipynb): The second step in the [Pixie](https://www.nature.com/articles/s41467-023-40068-5) pipeline is to run the cell clustering notebook. This notebook will use the pixel clusters generated in the first notebook to cluster the cells in the dataset. The notebook walks you through generating cell clusters for your data and generates cell cluster overlays.

[cell_neighbors_analysis-TONIC.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/cell_neighbors_analysis-TONIC.ipynb): There are two components of this notebook: neighborhood diversity and cell distance analysis. The diversity analysis can be used to determine the homogeneity/diversity of the neighbors surrounding each of the cells in our images. The cell distance section can be used to analyze the proximty/distance between cell phenotypes in samples.

[Calculate_Mixing_Scores-TONIC.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/Calculate_Mixing_Scores-TONIC.ipynb): The purpose of this notebook is to calculate the mixing score of a sample. Defined originally as a way to quantify the degree of mixing between specifically tumor and immune cells [Keren et al , Cell 2018](https://www.cell.com/cell/fulltext/S0092-8674(18)31100-0), the mixing score in this script has been adapted to analyze any two cell populations provided.

[example_neighborhood_analysis_script-TONIC.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/example_neighborhood_analysis_script-TONIC.ipynb): This neighborhood analysis notebook sheds light on neighborhoods made of micro-environments which consist of a collection of cell phenotypes, classified using kmeans clustering.

[example_fiber_segmentation-TONIC.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/example_fiber_segmentation-TONIC.ipynb): This is a notebook for segmenting fiber objects from collagen (or other) channel images. This notebook also provides visualizations and post-processing to extract relevant fiber features.

[ECM_Pixie_Cluster_Pixels.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/ECM_Pixie_Cluster_Pixels.ipynb): This notebook uses the pixel clustering algorithm from [Pixie](https://www.nature.com/articles/s41467-023-40068-5) to cluster ECM regions in the samples.

[ECM_pixel_clustering_stats.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/ECM_pixel_clustering_stats.ipynb): This notebook contains code to extract statistics from the ECM pixel clusters.


### feature_analysis

[1_postprocessing_cell_table_updates.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/1_postprocessing_cell_table_updates.py): This file takes the cell table generated by Pixie, and transforms it for plotting. Some of this functionality is 
has now been incorporated into [notebook 4](https://github.com/angelolab/ark-analysis/blob/main/templates/4_Post_Clustering.ipynb) in the [ark-analysis](https://github.com/angelolab/ark-analysis) repo. Other parts, however, have not yet been put into `ark`, such as aggregating cell populations. It also creates simplified cell tables
with only the necessary columns for specific plotting tasks.

[2_postprocessing_metadata.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/2_postprocessing_metadata.py): This file transforms the metadata files for analysis. It creates annotations in the metadata files that need to be computed from the
data, such as which patients have data from multiple timepoints.

[3_create_image_masks.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/3_create_image_masks.py): This file creates masks for each image based on supplied criteria. It identifies background based on the gold channel and tumor compartments based on ECAD staining patterns. It then takes these masks, and assigns each cell each image to the mask that it overlaps most with.

[4_ecm_preprocessing.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/4_ecm_preprocessing.py): This file creates ECM masks for each image based on the expression level of Collagen, Fibronectin, and FAP. We classified sections of each image as either Cold Collagen, Hot Collagen, or non-ECM, and then calculated the proportion of these classification in the image. 

[5_create_dfs_per_core.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/5_create_dfs_per_core.py): This file creates the dfs which will be used for plotting core-level information. It transforms the cell table into
a series of long-format dfs which can be easily used for data visualization. It creates separate dfs for cell population evaluations, functional marker
evaluation, etc.

[6_create_fov_stats.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/6_create_fov_stats.py): This file aggregates the various fov features and timepoint features into separate files, and  additionally filters out any unnecessary features based on their correlation within compartments.

[7_create_evolution_df.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/7_create_evolution_df.py): This file compares features across various timepoints and treatments. 

[8_genomics_postprocessing.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/8_genomics_postprocessing.py): This file processes genomics data and creates the feature table.


### multivariate_modeling

[process_feature_response_final.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/multivariate_modeling/process_feature_response_final.py): reformats the timepoint combined MIBI features file to prepare for per timepoint prediction modeling.
The same is done for [RNA](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/multivariate_modeling/process_feature_response_rna_signature_final.py) and [DNA](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/multivariate_modeling/process_feature_response_dna_final.py) feaures.

[all_timepoints.R](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/multivariate_modeling/all_timepoints.R): this script creates a Lasso model with 10 different random seeds for each individual timepoint. A 3-fold stratified cross-validation is used to select the level of sparsity in the Lasso model. The resulting AUROC was used as the primary metric for evaluating model performance.
The same is done for [RNA](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/multivariate_modeling/all_timepoints_rna_signature_final.R) and [DNA](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/multivariate_modeling/all_timepoints_dna_final.R) feaures.


### figure_generation
These scripts are used to generate Figures 1-5, as well as Extended Data Figures 1-12.


## Directory Structure
### Top Level Folders
`image_data`: Contains the single channel images for each FOV. 

`segmentation_data`: Contains the whole cell and nuclear segmentation masks for each FOV.

`analysis_files`: This directory should initially contain a cell table (generated with [Mesmer](https://www.nature.com/articles/s41587-021-01094-0) and annotated by [Pixie](https://www.nature.com/articles/s41467-023-40068-5)). The scripts expect a column named 
"cell_meta_cluster" containing the cell clusters, as well "fov" with the specific image name. 
This folder will also contain the final data tables generated by the TNBC scripts.

`output_files`: This directory will be created in [5_create_dfs_per_core.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/5_create_dfs_per_core.py) and store the per core and per timepoint data files for each feature. These will be aggregated to form the final data tables stored in *analysis_files*.

`intermediate_files`: This directory should contain subfolders storing any fov and cell level feature analysis done on the data. In addition, there should be a subdirectory containing the metadata
about each fov, each timepoint, and each patient, as appropriate for your study.

### Directory Tree
* TONIC_Cohort (base directory)
  * image_data 
  * segmentation_data
    * deepcell_output
  * analysis_files
  * output_files
  * intermediate_files
    * metadata
    * post_processing - contains specifications for the filtering of the data tables in *output_files*  
    * mask_dir - contains the compartment masks generated in [3_create_image_masks.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/3_create_image_masks.py)
    * fiber_segmentation_processed_data - image level fiber analysis ([code](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/example_fiber_segmentation-TONIC.ipynb))
      * tile_stats_512 - 512x512 tile analysis
    * spatial_analysis
      * dist_mats
      * neighborhood_mats - neighboring cell count/frequency at specified pixel radius and cell cluster level 
      * mixing_score - image level mixing score of various cell population combinations ([code](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/Calculate_Mixing_Scores-TONIC.ipynb))
      * cell_neighbor_analysis - data detailing cell diversity and linear distance between cell populations in an image ([code](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/cell_neighbors_analysis-TONIC.ipynb))
      * neighborhood_analysis - kmeans neighborhood analysis ([code](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/example_neighborhood_analysis_script-TONIC.ipynb))
    * ecm - generated in [4_ecm_preprocessing.py](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/feature_analysis/4_ecm_preprocessing.py)
    * ecm_pixel_clustering - generated in [ECM_Pixie_Cluster_Pixels.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/ECM_Pixie_Cluster_Pixels.ipynb) and [ECM_pixel_clustering_stats.ipynb](https://github.com/angelolab/publications/blob/main/2024-Greenwald_Nederlof_etal_TONIC/single_cell_analysis/ECM_pixel_clustering_stats.ipynb)



## Data Structures
In order to facilitate different analyses, there are a small number of distinct formats for storing data. 

*cell table*: This is the lowest level representation of the data, from which almost all other data formats are derived. Each row represents a single cell from a single image. Columns represent the different features for each cell. For example, the unique ID for each cell is located in the `label` column. The image that the cell came from is noted in the `fov` column, and the intensity of staining for CD68 protein is indicated by the `CD68` column. 
In addition, there are often multiple levels of granularity in the clustering scheme, which are represented here as different columns. For example, `cell_cluster` has more fine-grained assignments, with more distinct cell types, than `cell_cluster_broad`, which has a simpler schema. 

|  label  |    fov    | Ecadherin | CD68   | CD3   | cell_cluster | cell_cluster_broad |
|:-------:|:---------:|:---------:|:------:|:-----:|:------------:|:------------------:| 
|    1    | TMA1_FOV1 |    0.4    |  0.01  | 0.01  |    Cancer    |       Cancer       |
|    2    | TMA1_FOV1 |   0.01    |  0.0   |  0.8  |    T cell    |       Immune       | 
|   19    | TMA2_FOV4 |   0.01    |  0.8   | 0.01  |  Macrophage  |       Immune       | 

*segmentation mask*: This is the lowest level spatial representation of the data, from which most other spatial data formats are derived. Each image has a single segmentation mask, which has the locations of each cell. Cells are represented on a per-pixel basis, based on their `label` in the `cell_table`. For example, all of the pixels belonging to cell 1 would have a value of 1, all of the pixels belonging to cell 2 would have a value of 2, etc etc. Shown below is a simplified example, with cell 1 on the left and cell 2 on the right. 
```
0 0 0 0 0 0 0 0 0 0 
0 1 1 0 0 0 0 2 2 0 
1 1 1 1 0 0 2 2 2 2 
1 1 1 1 0 0 2 2 2 0 
1 1 0 0 0 0 0 2 2 0 
1 0 0 0 0 0 0 2 0 0 
0 0 0 0 0 0 0 0 0 0 
```

*distance_matrix.xr*: this data structure represents the distances between all cells in an image. The rows and columns are labeled according to the cell ID of each cell in an image, with the value at `ij`th position representing the euclidian distance, in pixels, between cell `i` and cell `j`.

|      | 1   |  3  |  6  |  8  | 
|:----:|:---:|:---:|:---:|:---:| 
|  1   |  0  | 200 | 30  | 21  | 
|  3   | 200 |  0  | 22  | 25  | 
|  6   | 30  | 22  |  0  | 300 | 
|  8   | 21  | 25  | 300 |  0  | 



*neighborhood_matrix*: This data structures summarizes information about the composition of a cell's neighbors. Each row represents an individual cell, with the columns representing the neighboring cells. For example, the first row would represent the number of cells of each cell type present within some pre-determined distance around the first cell in the image. 


|   fov     | label | cell_cluster | T cell | B cell | Macrophage | Treg | 
|:---------:|:-----:|:------------:|:------:|:------:|:----------:|:----:| 
| TMA1_FOV1 |   1   |    B cell    |   9    |   0    |     3      |  1   | 
| TMA1_FOV1 |   2   |     Treg     |   5    |   2    |     0      |  5   | 
| TMA2_FOV4 |   5   |    T cell    |   4    |   0    |     4      |  6   | 

## Analysis Files

*harmonized_metadata*: This data frame details the various FOVs and their associated tissue and patient IDs, timepoint, etc.

*feature_metadata*: This file gives more detailed information about the specifications that make up each of the features in the fov and timepoint feature tables. The columns include, general feature name, unique feature name, compartment, cell population, cell population level, and feature type details.

*timepoint_combined_features*: This dataframe details feature data for patients at various timepoints and includes the relevant metadata. It also includes evolution features, which describe the difference in feature values between two timepoints.

|         feature_name_unique          | raw_mean | normalized_mean | Patient_ID |            Timepoint            |                   combined_name                   |
|:------------------------------------:|:--------:|:---------------:|:----------:|:-------------------------------:|:-------------------------------------------------:|
|             area_Cancer              |   0.1    |       2.6       |     1      |  pre_treatement__on_treatment   |     area_Cancer__pre_treatement__on_treatment      |
| cluster_broad_diversity__cancer_core |  -0.01   |      -0.6       |     2      |          on_treatment           | cluster_broad_diversity_cancer_core__on_treatment |
|   max_fiber_density__stroma_border   |   -1.8   |      -0.7       |     3      |         pre_treatement          | max_fiber_density__stroma_border__pre_treatement  |


*combined_cell_table_normalized_cell_labels_updated*: The original cell table with all cell level data included. See the cell table description in [Data Structures](#Data-Structures) for more information. 

*cell_table_clusters*: Subset of the cell table containing just the FOV name, cell label, and different cluster labels.

*cell_table_counts*: Consolidated cell table with only marker count data.

*cell_table_morph*: Subset of the cell table containing only the morphological data for each cell (area, perimeter, major_axis_length, etc.). 

*cell_table_func_single_positive*: A cell table containing only the functional marker positivity data.

*cell_table_func_all*: A cell table containing all possible pairwise marker positivity data.

*fov_features*: This file is a combination of all feature metrics calculated on a per image basis. The file *fov_features_filtered* is also produced, which is the entire feature file with any highly correlated features removed.

The fov_features table aggregates features of many different types together, all of which are detailed in [Ouput Files](#Output-Files).

| Tissue_ID | fov | raw_value | normalized_value |   feature_name    |      feature_name_unique       |  compartment  | cell_pop |   feature_type   |
|:---------:|:---:|:---------:|:----------------:|:-----------------:|:------------------------------:|:-------------:|:--------:|:----------------:|
|    T1     |  1  |    0.1    |       2.6        | B__Cancer__ratio  |  B__Cancer__ratio_cancer_core  |  cancer_core  | multiple |  density_ratio   |
|    T2     |  2  |   -0.01   |       -0.6       | cancer_diversity  | cancer_diversity_cancer_border | cancer_border |  Cancer  | region_diversity |
|    T3     |  5  |   -1.8    |       -0.7       | max_fiber_density |       max_fiber_density        |  stroma_core  |   all    |      fiber       |

In the example table above, we see there are multiple columns that contain descriptive information about the statistics contained in each row. While `feature_name_unique` obviously gives the most granular description of the value, we can also use the other columns to quickly subset the data for specific analysis. 
For example, to look at all features within one region type across every image, we simply filter the `compartment` for only "cancer_core". 
Alternatively, we could compare the granular cell type diversity of all immune classified cells across regions by filtering both the `feature_type` as "cell_diversity" and `cell_pop` as "immune".


*timepoint_features*: While the data table above is aggregated *per_core*, this data is a combination of all feature metrics calculated on a per sample timepoint basis.  The file *timepoint_features_filtered* is also produced, which is the entire feature file with any highly correlated features removed.

| Tissue_ID |   feature_name    |      feature_name_unique       |  compartment  | cell_pop | raw_mean | raw_std | normalized_mean | normalized_std |
|:---------:|:-----------------:|:------------------------------:|:-------------:|:--------:|:-------:|:-------:|:---------------:|:--------------:|
|    T1     | B__Cancer__ratio  |  B__Cancer__ratio_cancer_core  |  cancer_core  | multiple |   0.1    |   1.3   |       2.6       |      0.3       |
|    T2     | cancer_diversity  | cancer_diversity_cancer_border | cancer_border |  Cancer  |  -0.01   |   0.3   |      -0.6       |      1.1       |
|    T3     | max_fiber_density |       max_fiber_density        |  stroma_core  |   all    |   -1.8   |   -16   |      -0.7       |      0.2       |

The file *timepoint_evolution_features* details the difference in feature values between two distinct timepoints from the same patient.


## Output Files

The individual feature data that combines into *fov_features* and *timepoint_features* can be found in the corresponding files detailed below.
Each of the data frames in this section can be further stratified based on the feature relevancy and redundancy. The files below can have any of the following suffixes:
* *_filtered*: features removed if there are less than 5 cells of the specified type
* *_deduped*: redundant features removed
* *_filtered_deduped*: both of the above filtering applied


1. *cluster_df*: This data structure summarizes key informaton about cell clusters on a per-image basis, rather than a per-cell basis. Each row represents a specific summary observation for a specific image of a specific cell type. For example, the number of B cells in a given image. The key columns are `fov`, which specifies the image the observation is from; `cell_type`, which specifies the cell type the observation is from; `metric`, which describes the specific summary statistic that was calculated; and `value`, which is the actual value of the summary statistic. For example, one statistic might be `cell_count_broad`, which would represent the number of cells per image, enumerated according the cell types in the `broad` clustering scheme. Another might be `cell_freq_detail`, which would be the frequency of the specified cell type out of all cells in the image, enumerated based on the detailed clustering scheme.

|   fov     | cell_type   | value |      metric       |    Timepoint    | 
|:---------:|:-----------:|:-----:|:-----------------:|:---------------:| 
| TMA1_FOV1 |   Immune    |  100  | cell_count_broad  |  pre_treatment  | 
| TMA1_FOV1 |    Treg     |  0.1  | cell_freq_detail  |  pre_treatment  | 
| TMA2_FOV4 | Macrophage  |  20   | cell_count_detail |  on_treatement  |  

In addition to these core columns, metadata can be added to facilitate easy analysis, such as disease stage, prognosis, anatomical location, or other information that is useful for plotting purposes. 


2. *functional_df*: This data structure summarizes information about the functional marker status of cells on a per-image basis. Each row represents the functional marker status of a single functional marker, in a single cell type, in a single image. The columns are the same as above, but with an additional `functional_marker` column which indicates which functional marker is being summarized. For example, one row might show the number of Tregs in a given image which are positive for Ki67, while another shows the proportion of cancer cells in an image that are PDL1+. 


|    fov    | cell_type  | value |      metric       | functional marker |     Timepont      | 
|:---------:|:----------:|:-----:|:-----------------:|:-----------------:|:-----------------:| 
| TMA1_FOV1 |   Immune   |  100  | cell_count_broad  |       Ki67        |   pre_treatment   | 
| TMA1_FOV1 |    Treg    |  0.4  | cell_freq_detail  |       PDL1        |   pre_treatment   | 
| TMA2_FOV4 | Macrophage |  20   | cell_count_detail |       TIM3        |   on_treatement   | 

3. *morph_df*: This data structure summarizes information about the morphology of cells on a per-image basis.  Each row represents the morphological statistic, in a single cell type, in a single image.

|    fov    | cell_type  |    value     | metric | functional marker |   Timepont    | 
|:---------:|:----------:|:------------:|:------:|:-----------------:|:-------------:| 
| TMA1_FOV1 |   Immune   |     area     |  100   | cell_count_broad  | pre_treatment | 
| TMA1_FOV1 |    Treg    | area_nuclear |  0.4   | cell_freq_detail  | pre_treatment | 
| TMA2_FOV4 | Macrophage |   nc_ratio   |   20   | cell_count_detail | on_treatement | 

4. *distance_df*: This data structure summarizes information about the closest linear distance between cell types on a per-image basis.

|   fov     |    cell_type     | linear_distance | value |      metric        |   Timepoint    | 
|:---------:|:----------------:|:---------------:|:-----:|:------------------:|:--------------:| 
| TMA1_FOV1 |      Immune      |     Immune      |  100  | cluster_broad_freq | pre_treatement | 
| TMA1_FOV1 |      Immune      |      Treg       |  0.4  | cluster_broad_freq | pre_treatement | 
| TMA2_FOV4 | MacImmunerophage |   Macrophage    |  20   | cluster_broad_freq | on_treatement  | 


5. *diversity_df*: This data structure summarizes information about the diversity of cell types on a per-image basis.

|   fov     |    cell_type     |      diversity_feature       | value |      metric        |     Timepoint      | 
|:---------:|:----------------:|:----------------------------:|:-----:|:------------------:|:------------------:| 
| TMA1_FOV1 |      Immune      | diversity_cell_cluster_broad |  1.1  | cluster_broad_freq |   pre_treatement   | 
| TMA1_FOV1 |      Immune      |    diversity_cell_cluster    |  0.4  | cluster_broad_freq |   pre_treatement   | 
| TMA2_FOV4 | MacImmunerophage | diversity_cell_cluster_broad |   2   | cluster_broad_freq |   on_treatement    | 

6. *fiber_df / fiber_df_per_tile*: This data structure summarizes statistics about the collagen fibers at an image-level and also within 512x512 sized pixel crops of the image.

| Tissue_ID |      fiber_metric       | mean | std |    Timepoint     | 
|:---------:|:-----------------------:|:----:|:---:|:----------------:| 
| TMA1_FOV1 |  fiber_alignment_score  | 2.2  | 0.5 |  pre_treatement  | 
| TMA1_FOV1 |        fiber_are        | 270  | 30  |  pre_treatement  | 
| TMA2_FOV4 | fiber_major_axis_length |  35  |1.9  |  on_treatement   | 

7. *neighborhood_image_proportions / neighborhood_compartment_proportions*: These data files detail the proportion of cells assigned to each kmeans cluster in the image / in each compartment in each image. 

8. *formatted_mixing_scores*: This file contains the mixing scores calculated per image for various cell population combinations.