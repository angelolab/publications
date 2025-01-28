library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(tidyr)
library(data.table)
library(arrow)
combined_data <- fread("/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv")

filtered_data <- combined_data %>% 
  select(sample_id, label, combined_label, `centroid-0`,`centroid-1`, cell_meta_cluster_final, cell_meta_cluster_IT, cell_meta_cluster_ML, table_origin) %>% 
  rename(fov = sample_id,
         original_label = label,
         label = combined_label,
         cell_meta_cluster = cell_meta_cluster_final)

# file_path_combined <- "H://Shared drives//BRUCE_data//cell_tables//combined_panel//20240521_filtered_spatial_combined_panel.csv"
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_spatial_notebook_prepped_combined_panel.csv"
# Save data_Immune as a CSV file
fwrite(filtered_data, file_path_combined, row.names = FALSE)

