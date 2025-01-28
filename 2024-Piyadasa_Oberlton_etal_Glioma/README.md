# Piyadasa_Obelton_etal_Glioma

## Supplemental Material:
- Low Level Image Processing: [Toffy](https://github.com/angelolab/toffy)
- Cell Segmentation, Cell Clustering: [Ark](https://github.com/angelolab/ark-analysis)
- Updated QUICHE: QUICHE
- Tables and Metadata: [Data Portal](https://www.bruce.parkerici.org/pages/raw-data-access.html)

## Directories:
- Figure Generation: Code Provided for generating all Main and Supplemental Figures
- Data Analysis:
  - [1_Tumor_Immune_Spatial_Alignment](https://github.com/benoberlton/Piyadasa_Obelton_etal_Glioma/blob/main/Data_Analysis/1_Tumor_Immune_Spatial_Alignment.R): Takes the tumor and immune tables and spatially aligns each of the unique FOVs, keeping the tumor cells from the tumor table and the immune cells from the immune table to create a final combined cell table.
  - [4a_Spatial_Notebook_Prep](https://github.com/benoberlton/Piyadasa_Obelton_etal_Glioma/blob/main/Data_Analysis/4a_Spatial_Notebook_Prep.R): Renames and simplifies the combined cell table to be used in [4b_Cell_neighbors_analysis](https://github.com/benoberlton/Piyadasa_Obelton_etal_Glioma/blob/main/Data_Analysis/4b_Cell_neighbors_analysis.ipynb).
  - [4b_Cell_neighbors_analysis](https://github.com/benoberlton/Piyadasa_Obelton_etal_Glioma/blob/main/Data_Analysis/4b_Cell_neighbors_analysis.ipynb): Used to determine the homogeneity/diversity of the neighbors surrounding each of the cells within each image. Additionally calculates the proximty/distance between cell phenotypes in each FOV.
  - [6_Tumor_Antigen_GSEA](https://github.com/benoberlton/Piyadasa_Obelton_etal_Glioma/blob/main/Data_Analysis/6_Tumor_Antigen_GSEA.R): Performs GSEA (Gene Set Encrichment Analysis) to calculate the enrichment of tumor antigens across various tumor conditions.
  - [7_QUICHE](https://github.com/benoberlton/Piyadasa_Obelton_etal_Glioma/tree/main/Data_Analysis/7_QUICHE): Contains all necessary code to run the version of QUICHE used for this study. QUICHE executes creates spatial cellular neighborhoods followed by differential abundance of these spatiall neighborhoods between research conditions.
  - 8_Random_Forest_Model: Performs Random Forest Modelling and PCA analysis on proteomic, genomic, glycomic, and/or spatial features.
