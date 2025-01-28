# 20250121
# Authors - Hadeesha Piyadasa

# load required packages -----
library(tidyverse)
library(forcats)
library(data.table)
library(patchwork)
library(emmeans)
library(pheatmap)
library(gplots)
library(viridis)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)
library(ggbreak)
library(combinat)
library(tiff)
library(arrow)
library(readxl)
library(RColorBrewer)
library(umap)
library(stringr)
library(reshape2)
library(GGally)
library(tidyr)
library(scales)
library(ggh4x)
library(webr)
library(ggforce)
library(plotly)
library(progress)
library(MASS)
library(tictoc)

options(arrow.unsafe_metadata = TRUE)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### functions ######
combine_data_metadata <- function(data, metadata, merge_by, remove_duplicates = TRUE) {
  # Ensure that metadata only contains rows that exist in the data
  metadata <- metadata %>% filter(!!sym(merge_by) %in% data[[merge_by]])
  
  # Find duplicate columns in metadata that are also in data except for merge_by
  duplicate_columns <- setdiff(intersect(colnames(metadata), colnames(data)), merge_by)
  
  # Remove the duplicate columns from metadata
  if (length(duplicate_columns) > 0) {
    metadata <- metadata %>% dplyr::select(-one_of(duplicate_columns))
  }
  
  # Merge data and metadata using the column specified by merge_by
  data_merged <- data %>% left_join(metadata, by = merge_by)
  
  # Optionally remove duplicate rows based on the merge_by column
  if (remove_duplicates) {
    data_merged <- data_merged %>% distinct(!!sym(merge_by), .keep_all = TRUE)
  }
  
  return(data_merged)
}

calculate_all_cell_type_count <- function(data, cell_meta_cluster_col, cluster_value, summary_by) {
  all_cell_type_count <- data %>%
    group_by(!!sym(summary_by)) %>%
    summarise(all_cell_type_count = sum(!!sym(cell_meta_cluster_col) == cluster_value, na.rm = TRUE)) %>%
    ungroup()
  
  return(all_cell_type_count)
}

calculate_positive_counts <- function(data, columns, summary_by, cell_meta_cluster_column, cell_meta_cluster_value = NULL, summary_by_2 = NULL) {
  # Filter for rows where cell_meta_cluster_column matches cell_meta_cluster_value, if provided
  if (!is.null(cell_meta_cluster_value)) {
    filtered_data <- data %>%
      filter(!!sym(cell_meta_cluster_column) == cell_meta_cluster_value)
  } else {
    filtered_data <- data
  }
  
  # Check if summary_by_2 is provided
  if (is.null(summary_by_2)) {
    # Group by only the first summary variable
    grouped_data <- filtered_data %>%
      group_by(!!sym(summary_by))
  } else {
    # Group by both summary variables
    grouped_data <- filtered_data %>%
      group_by(!!sym(summary_by), !!sym(summary_by_2))
  }
  
  # Summarise the counts for each group
  non_empty_counts <- grouped_data %>%
    summarise(across(all_of(columns), ~sum(. == 1, na.rm = TRUE), .names = "{.col}_counts"))
  
  return(non_empty_counts)
}

cleanup_data_by_metadata <- function(data, metadata, join_by_col, summary_by) {
  # Replace NA with 0
  data[is.na(data)] <- 0
  
  # Remove duplicates from tumor_metadata
  metadata_unique <- metadata %>%
    distinct(summary_by, .keep_all = TRUE)
  
  # Merge with cell_prop_data
  merged_data <- data %>%
    left_join(metadata_unique, by = join_by_col)
  
  # Apply filters and mutations
  processed_data <- merged_data %>%
    filter(!(final_diagnosis %in% c("Meningioma", "mFTC_brain","Breast_CA","Colon_CA","FTC","Placenta","Tonsil"))) %>%
    filter(!is.na(site)) %>%
    mutate(final_diagnosis = if_else(site == "CHOP" & final_diagnosis == "GBM", "pGBM", final_diagnosis))
  
  return(processed_data)
}

compare_columns_multi <- function(columns, cell_table_tumor_thresholded) {
  # Create a new column name based on the combination
  new_col_name <- paste(columns, collapse = "_")
  
  # Construct the logical condition for the combination with binary values
  condition_expr <- paste(sapply(columns, function(col) paste0("get('", col, "') == 1")), collapse = " & ")
  
  # Create the new column based on the condition
  expression_string <- paste0("(new_col_name) := ifelse(", condition_expr, ", 1, 0)")
  eval(parse(text = paste0("cell_table_tumor_thresholded[, ", expression_string, "]")))
}

generate_density_feature <- function(fov_area_df, data, metadata, feature_type, feature_source, merge_column_metadata, merge_column_data) {
  # Add fov info back
  data <- merge(data, metadata[, c(merge_column_metadata, "sample_id")], 
                by.x = "sample_id", by.y = "sample_id", all.x = TRUE)
  
  # Rename the merge column to 'fov' in the merged data
  colnames(data)[which(names(data) == merge_column_metadata)] <- "fov"
  # Merge with fov_area_df - Assuming fov_area_df is available in the environment
  area_data <- merge(fov_area_df, data, by.x = merge_column_data, by.y = "fov")
  # Rename the merge column to 'fov' in the merged data
  colnames(area_data)[which(names(area_data) == merge_column_data)] <- "fov"

  # Identify columns for density calculation
  ignored_columns <- c(merge_column_data, "sample_id", "height", "width", "area_in_microns", "cell_meta_cluster_final")
  area_colnames <- setdiff(colnames(area_data), ignored_columns)
  denominator <- "area_in_microns"

  # Loop to find the density and remove old columns
  for (num_col in area_colnames) {
    new_col_name <- paste(sub("_counts$", "", num_col), "density", sep = "_")
    area_data[[new_col_name]] <- area_data[[num_col]] / area_data[[denominator]]
    area_data[[num_col]] <- NULL
  }
  
  # Appending "_density" to each element of area_colnames
  feature_columns <- setdiff(colnames(area_data), ignored_columns)
  
  # Update the data frame with feature_type and feature_source
  area_data <- area_data %>%
    mutate(feature_type = feature_type,
           feature_source = feature_source)
  # Combine with metadata
  area_data <- merge(area_data, metadata[, c(merge_column_metadata, "patient_id")], by.x = "fov", by.y = merge_column_metadata, all.x = TRUE)
  
  # Convert to data.table format
  setDT(area_data)
  # Determine id.vars based on feature_source
  id_vars_base <- c('sample_id', "fov", 'patient_id', 'feature_source', 'feature_type')
  id_vars <- if (feature_source == "cell_meta_cluster_final") {
    c(id_vars_base, 'cell_meta_cluster_final')
  } else {
    id_vars_base
  }
  
  # Melt the data
  melted_area_data <- data.table::melt(data = area_data,
                                       id.vars = id_vars,
                                       measure.vars = feature_columns,
                                       variable.name = 'feature_variable',
                                       value.name = 'feature_value')
  
  # Return the melted data frame
  return(melted_area_data)
}

calculate_proportion_immune_func <- function(sample_id_data, cell_type_1, marker_1, cell_type_2, marker_2) {
  # Get the counts for the first cell type and marker combination
  count_1 <- sum(as.numeric(sample_id_data[sample_id_data$cell_meta_cluster_final == cell_type_1, paste0(marker_1, "_func_counts")]), na.rm = TRUE)
  
  # Get the counts for the second cell type and marker combination
  count_2 <- sum(as.numeric(sample_id_data[sample_id_data$cell_meta_cluster_final == cell_type_2, paste0(marker_2, "_func_counts")]), na.rm = TRUE)
  
  # Calculate the proportions
  proportion_A = ifelse((count_1 + count_2) > 10, count_1 / (count_1 + count_2), NA)
  proportion_B = ifelse((count_1 + count_2) > 10, count_2 / (count_1 + count_2), NA)
  # Calculate the proportions only if count_1 and count_2 are both greater than or equal to 5
  # proportion_A <- ifelse(count_1 >= 5 & count_2 >= 5, count_1 / (count_1 + count_2), NA)
  # proportion_B <- ifelse(count_1 >= 5 & count_2 >= 5, count_2 / (count_1 + count_2), NA)
  
  return(list(proportion_A, proportion_B))
}

# Function to check if the value exists in the given column in master_feature_table
value_exists_in_column <- function(value, column_name, data_table) {
  value %in% data_table[[column_name]]
}

# Updated build_row function to check existence for all specified columns and give feedback
build_row <- function(number, filter_col1, filter_col1_value, filter_col2, filter_col2_value, filter_col3 = NA, filter_col3_value = NA, summary_by = NA, group_feature_col = NA, group_1 = NA, group_2 = NA, data_table) {
  
  # A helper function to check value existence and provide feedback if necessary
  check_value_existence_and_provide_feedback <- function(col_value, col_name) {
    if (!is.na(col_value) && !value_exists_in_column(col_value, col_name, data_table)) {
      message(paste("Warning: The value", col_value, "does not exist in the column", col_name, ". Available values are:", paste(unique(data_table[[col_name]]), collapse = ", ")))
      return(TRUE)
    }
    return(TRUE)
  }
  
  # Perform the checks
  if (!check_value_existence_and_provide_feedback(filter_col1_value, filter_col1)) return(NULL)
  if (!check_value_existence_and_provide_feedback(filter_col2_value, filter_col2)) return(NULL)
  if (!check_value_existence_and_provide_feedback(filter_col3_value, filter_col3)) return(NULL)
  if (!check_value_existence_and_provide_feedback(group_1, group_feature_col)) return(NULL)
  if (!check_value_existence_and_provide_feedback(group_2, group_feature_col)) return(NULL)
  
  # Return a named vector that represents a row in the new data frame
  return(c(number, filter_col1, filter_col1_value, filter_col2, filter_col2_value, filter_col3, filter_col3_value, summary_by, group_feature_col, group_1, group_2))
}

# Modify the new_row function to replace an existing row if the same row_number is used
new_row <- function(row_number, filter_col1, filter_col1_value, filter_col2, filter_col2_value, filter_col3, filter_col3_value, summary_by, group_feature_col, group_1, group_2, data_table, new_feature_table, new_df_cols) {
  # Check if the column names provided exist in master_feature_table
  row <- tryCatch({
    build_row(
      number = row_number,
      filter_col1 = filter_col1,
      filter_col1_value = filter_col1_value,
      filter_col2 = filter_col2,
      filter_col2_value = filter_col2_value,
      filter_col3 = filter_col3,
      filter_col3_value = filter_col3_value,
      summary_by = summary_by,
      group_feature_col = group_feature_col,
      group_1 = group_1,
      group_2 = group_2,
      data_table = data_table
    )
  }, error = function(e) {
    message(e$message)
    NULL  # Return NULL if there's an error
  })
  
  # Check if the new row was created successfully
  if (!is.null(row)) {
    # If 'new_feature_table' does not exist in the global environment, create it
    if (!exists("new_feature_table", envir = .GlobalEnv)) {
      assign("new_feature_table", setNames(data.frame(matrix(ncol = length(new_df_cols), nrow = 0)), new_df_cols), envir = .GlobalEnv)
    }
    
    # Check if the row_number already exists
    if (row_number <= nrow(get("new_feature_table", envir = .GlobalEnv))) {
      # Replace the existing row
      new_feature_table <- get("new_feature_table", envir = .GlobalEnv)
      new_feature_table[row_number, ] <- setNames(as.data.frame(t(row)), new_df_cols)
      assign("new_feature_table", new_feature_table, envir = .GlobalEnv)
    } else {
      # Add the new row if the row_number does not exist
      new_feature_table <- get("new_feature_table", envir = .GlobalEnv)
      new_feature_table <- rbind(new_feature_table, setNames(as.data.frame(t(row)), new_df_cols))
      assign("new_feature_table", new_feature_table, envir = .GlobalEnv)
    }
    print(get("new_feature_table", envir = .GlobalEnv))
  } else {
    message("The new row could not be created because one or more values do not exist in the specified columns.")
  }
}

calculate_adjusted_proportions <- function(data, marker_columns, cell_meta_cluster_col, cluster_value, summary_by) {
  # Initialize an empty dataframe to store the adjusted proportions
  adjusted_df <- data.frame(sample_id = unique(data[[summary_by]]))
  
  # Loop through each marker column to adjust proportions
  for (specific_marker in marker_columns) {
    print(specific_marker)
    # Create a subset of data excluding the specific marker
    data_without_specific_marker <- data %>%
      filter(!(!!sym(specific_marker) == 1 & !!sym(cell_meta_cluster_col) == cluster_value))
    
    # Calculate the total tumor cells for each sample excluding the specific marker
    all_tumor_count <- data_without_specific_marker %>%
      group_by_at(vars(!!sym(summary_by))) %>%
      summarise(all_tumor_count = sum(!!sym(cell_meta_cluster_col) == cluster_value, na.rm = TRUE)) %>%
      ungroup()
    
    # Loop through each marker to calculate adjusted proportions
    for (marker in marker_columns) {
      if (marker != specific_marker) {
        # Calculate the sum of each marker for each sample, excluding the specific marker
        specific_marker_count <- data_without_specific_marker %>%
          filter(!!sym(cell_meta_cluster_col) == cluster_value) %>%
          group_by_at(vars(!!sym(summary_by))) %>%
          summarise(marker_count = sum(!!sym(marker), na.rm = TRUE)) %>%
          ungroup()
        #print(specific_marker_count)
        # Merge the specific marker count with the all tumor count
        merged_counts <- merge(all_tumor_count, specific_marker_count, by = summary_by, all.x = TRUE)
        
        # Calculate the adjusted proportion
        adjusted_prop_col <- paste(marker, specific_marker, "neg_prop", sep = "_")
        merged_counts[[adjusted_prop_col]] <- merged_counts$marker_count / merged_counts$all_tumor_count
        
        # Select only the necessary columns
        merged_counts <- merged_counts[, c(summary_by, adjusted_prop_col)]
        
        # Merge this into the adjusted_df dataframe
        adjusted_df <- merge(adjusted_df, merged_counts, by = summary_by, all = TRUE)
      }
    }
  }
  
  # Replace any NAs with zero
  adjusted_df[is.na(adjusted_df)] <- 0
  
  # Return the adjusted dataframe
  return(adjusted_df)
}

calculate_adjusted_counts <- function(data, marker_columns, cell_meta_cluster_col, cluster_value, summary_by) {
  list_counts_dfs <- list()
  
  # Loop through each marker column to get counts
  for (specific_marker in marker_columns) {
    # Create a subset of data excluding the specific marker
    data_without_specific_marker <- data %>%
      filter(!(!!sym(specific_marker) == 1 & !!sym(cell_meta_cluster_col) == cluster_value))
    
    # Calculate the total tumor cells for each sample excluding the specific marker
    all_tumor_count <- data_without_specific_marker %>%
      group_by_at(vars(!!sym(summary_by))) %>%
      summarise(all_tumor_count = sum(!!sym(cell_meta_cluster_col) == cluster_value, na.rm = TRUE)) %>%
      ungroup()
    
    # Initialize a dataframe to store counts including all_tumor_count
    counts_df <- all_tumor_count
    
    # Loop through each marker to get counts
    for (marker in marker_columns) {
      # Create a column name based on the marker
      count_col_name <- paste(marker, "count", sep = "_")
      
      # Calculate the sum of each marker for each sample, excluding the specific marker
      specific_marker_count <- data_without_specific_marker %>%
        filter(!!sym(cell_meta_cluster_col) == cluster_value) %>%
        group_by_at(vars(!!sym(summary_by))) %>%
        summarise(!!count_col_name := sum(!!sym(marker), na.rm = TRUE)) %>%
        ungroup()
      
      # Merge the specific marker count with the all tumor count
      counts_df <- merge(counts_df, specific_marker_count, by = summary_by, all.x = TRUE)
    }
    
    # Replace any NAs with zero
    counts_df[is.na(counts_df)] <- 0
    
    # Assign the dataframe to the list with the name of the specific marker
    list_counts_dfs[[specific_marker]] <- counts_df
  }
  
  # Return the list of dataframes
  return(list_counts_dfs)
}
#### import metadata ######
metadata <- read.csv("your_path...../metadata_complete.csv") # download from www.bruce.parkerici.org
#### load all features calculated below (skip if files have not been generated) ###### 

cell_table_tumor_metadata <- read_parquet("../tables/cell_table_tumor_metadata.parquet")

cell_table_immune_metadata <- read_parquet("../tables/cell_table_immune_metadata.parquet")

cell_table_all_merged <- read_parquet("../tables/cell_table_all_merged.parquet")

cell_table_all_merged_thresholded <- read_parquet("../tables/cell_table_all_merged_thresholded.parquet")

cell_table_tumor_thresholded <- read_parquet("../tables/cell_table_tumor_thresholded.parquet")

load("../tables/complete_tumor_columns.RData")

cell_table_immune_thresholded <- read_parquet("../tables/cell_table_immune_thresholded.parquet")

load("../tables/complete_immune_func_columns.RData")

fov_area_df_tumor <- read_parquet("../tables/fov_area_df_tumor.parquet")

fov_area_df_immune <- read_parquet("../tables/fov_area_df_immune.parquet")

tumor_count_result <- read_parquet("../tables/tumor_count_result.parquet")

immune_count_result <- read_parquet("../tables/immune_count_result.parquet")

all_cell_count_tumor_FOV_result <- read_parquet("../tables/all_cell_count_tumor_FOV_result.parquet")

all_cell_count_immune_FOV_result <- read_parquet("../tables/all_cell_count_immune_FOV_result.parquet")

positive_count_result_immune_cells <- read_parquet("../tables/positive_count_result_immune_cells.parquet")

positive_count_result_tumor <- read_parquet("../tables/positive_count_result_tumor.parquet")

tumor_antigen_segment_df <- read_parquet("../tables/tumor_antigen_segment_df.parquet")

melted_cell_table_marker_intensity <- read_parquet("../tables/melted_cell_table_marker_intensity.parquet")

melted_cell_table_tumor_cell_ratios <- read_parquet("../tables/melted_cell_table_tumor_cell_ratios.parquet")

tumor_cell_type_count <- read_parquet("../tables/tumor_cell_type_count.parquet")

melted_cell_table_immune_cell_ratios <- read_parquet("../tables/melted_cell_table_immune_cell_ratios.parquet")

immune_cell_type_count <- read_parquet("../tables/immune_cell_type_count.parquet")

melted_tumor_antigen_coverage <- read_parquet("../tables/melted_tumor_antigen_coverage.parquet")

melted_tumor_antigen_coverage_amounts <- read_parquet("../tables/melted_tumor_antigen_coverage_amounts.parquet")

melted_tumor_antigen_coexp_counts <-read_parquet("../tables/melted_tumor_antigen_coexp_counts.parquet")

melted_immune_func_ratios <- read_parquet("../tables/melted_immune_func_ratios.parquet")

immune_func_type_count <- read_parquet("../tables/immune_func_type_count.parquet")

melted_immune_func_to_all <- read_parquet("../tables/melted_immune_func_to_all.parquet")

melted_tumor_cell_type_area <- read_parquet("../tables/melted_tumor_cell_type_area.parquet")

melted_immune_func_cell_type_area <- read_parquet("../tables/melted_immune_func_cell_type_area.parquet")

melted_immune_cell_type_area <- read_parquet("../tables/melted_immune_cell_type_area.parquet")

melted_all_tumor_area <- read_parquet("../tables/melted_melted_all_tumor_area.parquet")

melted_all_immune_area <- read_parquet("../tables/melted_all_immune_area.parquet")

melted_tumor_antigen_segment_df <- read_parquet("../tables/melted_tumor_antigen_segment_df.parquet")

melted_cell_table_sptial_no_metadata <- read_parquet("../tables/melted_cell_table_sptial_no_metadata.parquet")

NS_melted <- read_parquet("../tables/NS_melted.parquet")

MALDI_melted <- read_parquet("../tables/20241112_MALDI_BRUCE_melted.parquet")

melted_neighborhood_freqs_50 <- read_parquet("../tables/melted_neighborhood_freqs_50.parquet")

melted_neighborhood_freqs_ML_radius50 <- read_parquet("../tables/melted_neighborhood_freqs_ML_radius50.parquet")

melted_neighborhood_freqs_IT_radius50 <- read_parquet("../tables/melted_neighborhood_freqs_IT_radius50.parquet")

master_feature_table_filtered <- readRDS("../tables/20240810_master_feature_table_na_removed.rds")

spatial_density_melted <- read_parquet("../tables/spatial_density_melted.parquet")

melted_cell_table_immune_cell_ratios_broad <- read_parquet("../tables/melted_cell_table_immune_cell_ratios_broad.parquet")

percentage_positive_melted <- read_parquet("../tables/percentage_positive_melted.parquet")

#### lists for filtering ######

tumor_columns <- setdiff(colnames(cell_table_tumor_metadata)[which(colnames(cell_table_tumor_metadata) == "ApoE"):which(colnames(cell_table_tumor_metadata) == "VISTA")], c("Au", "Fe", "Noodle", "Nuclear","GD2","IL13RA2"))

#targetabble tumor antigens
tumor_columns_core <- c("B7H3", 
                        "EGFR",  
                        "GM2_GD2", "GPC2",  
                        "HER2", 
                        "NG2", "VISTA")


# different combinations of tumor marker columns that were thresholded manually. These are binary columns indicating positive or negative status
tumor_func_columns <- c("ApoE_func", "B7H3_func", "EGFR_func", "EGFRvIII_func", "GFAP_func", "GM2_GD2_func", "GPC2_func", "H3K27M_func", "HER2_func", "IDH1_R132H_func", 
                        "NG2_func", "Olig2_func", "VISTA_func", "H3K27me3_func","Ki67_func")

tumor_func_columns_core <-  c("ApoE_func", "B7H3_func", "EGFR_func" , "GFAP_func", "GM2_GD2_func" , "GPC2_func" ,  "HER2_func", "IDH1_R132H_func","NG2_func","Olig2_func" ,"VISTA_func")

tumor_func_columns_core_simple = c("B7H3_func", "EGFR_func", "GM2_GD2_func" , "GPC2_func" ,  "HER2_func","NG2_func","VISTA_func")

tumor_func_columns_core_simple_segment = c("B7H3_segment", "EGFR_segment", "GM2_GD2_segment" , "GPC2_segment" ,  "HER2_segment","NG2_segment","VISTA_segment")

tumor_columns_core_simple <- sub("_func$", "", tumor_func_columns_core_simple)

selected_tumor_columns = c("B7H3_func","EGFR_func","GM2_GD2_func")


# immune columns
immune_columns <- setdiff(colnames(cell_table_immune_metadata)[which(colnames(cell_table_immune_metadata) == "Arginase1"):which(colnames(cell_table_immune_metadata) == "Tox")], c("Au", "Fe", "CD208","chan_39","chan_48","chan_70","Noodle", "Nuclear","GD2","IL13RA2"))
immune_columns

cell_type_immune <- cell_table_all_merged %>%
  filter(source == "cell_table_immune") %>%
  distinct(cell_meta_cluster_final) %>%
  pull(cell_meta_cluster_final)

# Exclude specific values
cell_type_immune <- setdiff(cell_type_immune, c("Tumor_cells", "Unassigned","Neurons","Endothelial_cells"))

cell_type_immune_list <- c(cell_type_immune)
cell_type_immune

#immune functional markers
immune_func_marker <- c("TIM3", "CD38", "Tox", "iNOS", "CD86", "Ki67", "PDL1", "LAG3", "PD1", "ICOS", "IDO1", "GLUT1","Arginase1")

#immune markers that were thresholded manually
immune_func_columns = paste0(immune_func_marker, "_func")
immune_func_columns <- c("TIM3_func", "CD38_func", "Tox_func", "iNOS_func", "CD86_func", "Ki67_func", "PDL1_func", "LAG3_func", "PD1_func", "ICOS_func", "IDO1_func", "GLUT1_func","CD31_func","Arginase1_func")
immune_func_columns_core <- c("CD38_func", "iNOS_func", "CD86_func", "PDL1_func", "PD1_func", "ICOS_func", "IDO1_func","Arginase1_func")

# dictionary for calculating phenotype marker intensity  by cell type....

cell_marker_dictionary <- list(
  Tumor_cells = c("GFAP", "Olig2", "CD133"),  # Updated from tumor_cells
  Tcell_CD4 = c("CD45","CD3", "CD4"),  # Updated from CD4_Tcells
  Tcell_CD8 = c("CD45","CD3", "CD8"),  # Updated from CD8_Tcells
  Myeloid_CD141 = c("CD45","CD14", "CD68", "CD163", "CD141", "CD11b", "HLADR"),  # Updated from CD141_Myeloid
  Myeloid_CD14 = c("CD45","CD14", "CD68", "CD11b", "CD163"),  # Updated from CD14_Myeloid
  Myeloid_CD14_CD163 = c("CD45","CD14", "CD68", "CD11b", "CD163"),  # Updated from CD14_Myeloid
  Microglia = c("CD45","CD14", "CD68", "CD163", "CD11b", "TMEM119", "HLADR"),
  Microglia_CD163 = c("CD45","CD14", "CD68", "CD163", "CD11b", "TMEM119", "HLADR"),
  APC = c("CD45","CD14", "CD68", "CD163", "CD141", "CD11b", "HLADR"),
  Macrophage_CD206 = c("CD45","CD14", "CD68", "CD163", "CD206", "CD11b", "HLADR"),  # Updated from CD206_Macrophages
  Tcell_FoxP3 = c("CD45","CD3", "CD4", "FoxP3"),  # Updated from FoxP3_Tcells
  Neutrophils = c("CD45","CD11b", "Calprotectin"),
  Neurons = c("NeuN"),
  DC_Mac_CD209 = c("CD45","CD14", "CD68", "CD163", "CD209", "CD11b", "HLADR"),  # Updated from CD209_DC_Mac
  Macrophage_CD68_CD163 = c("CD45","CD14", "CD68", "CD163", "CD11b", "HLADR"),  # Updated from CD68_Myeloid
  Macrophage_CD68 = c("CD45","CD14", "CD68", "CD163", "CD11b", "HLADR"),  # Updated from CD68_Myeloid
  Immune_unassigned = c("CD45"),  # Updated from Immune_other
  Bcells = c("CD45","CD20", "HLADR"),
  Mast_Cells = c("CD45","Chym_Tryp"),
  Myeloid_CD11b_HLADRplus = c("CD45","CD14", "CD68", "CD163", "CD11b", "HLADR"),  # Newly added
  Myeloid_CD11b_HLADRminus = c("CD45","CD14", "CD68", "CD163", "CD11b", "HLADR"),  # Newly added
  Endothelial_cells = c("CD31")
  
)

# dictionary for calculating functional marker intensity  by cell type....

universal_markers <- c("GLUT1", "Ki67", "PDL1")

cell_func_dictionary <- list(
  Immune_unassigned = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Macrophage_CD68 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Tcell_CD8 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Myeloid_CD141 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Tumor_cells = c("Tox", "IDO1", "PDL1",  "HLA1", "ApoE", "B7H3", "CD133", "EGFR", "EGFRvIII",  "GFAP", "GM2_GD2", "GPC2", "H3K27M", "H3K27me3", "HER2","IDH1_R132H", "Ki67", "NG2", "Olig2", "VISTA"),  
  Myeloid_CD14 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Microglia = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  APC = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Macrophage_CD206 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Tcell_CD4 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Tcell_FoxP3 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Neutrophils = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),   
  Neurons = c(),
  DC_Mac_CD209 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Bcells = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),   
  Mast_Cells = c(),
  Endothelial_cells = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67","CD31"), 
  Myeloid_CD14_CD163 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Microglia_CD163 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Myeloid_CD11b_HLADRplus = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Myeloid_CD11b_HLADRminus = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Macrophage_CD68_CD163 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67")
)


# Combine the dictionaries and add universal markers
combined_cell_marker_dictionary <- cell_marker_dictionary

for (key in names(cell_func_dictionary)) {
  if (key %in% names(combined_cell_marker_dictionary)) {
    combined_cell_marker_dictionary[[key]] <- unique(c(combined_cell_marker_dictionary[[key]], cell_func_dictionary[[key]], universal_markers))
  } else {
    combined_cell_marker_dictionary[[key]] <- unique(c(cell_func_dictionary[[key]], universal_markers))
  }
}

# Print the combined dictionary
print(combined_cell_marker_dictionary)


# Creating a dictionary with cell types as keys and functional markers as values (only for immune panel)

func_immune_dictionary <- list(
  Immune_unassigned = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Macrophage_CD68 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Macrophage_CD68_CD163 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Tcell_CD8 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),    
  Myeloid_CD141 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Tumor_cells = c("Tox", "IDO1", "PDL1"),  
  Myeloid_CD14 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Myeloid_CD14_CD163 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Microglia = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Microglia_CD163 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  APC = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Macrophage_CD206 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Tcell_CD4 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Myeloid_CD11b_HLADRminus = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Myeloid_CD11b_HLADRplus = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1"), 
  Tcell_FoxP3 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  DC_CD123 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Neutrophils = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),  
  Neurons = c(),
  DC_Mac_CD209 = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"),   
  Bcells = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Mast_cells = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67"), 
  Unassigned = c(),
  Endothelial_cells = c("TIM3", "CD38", "Tox", "iNOS", "CD86", "LAG3", "PD1", "ICOS", "IDO1","Arginase1", "GLUT1","Ki67","CD31") 
)


# all "functional" markers
all_func_cols <- c(tumor_func_columns,immune_func_columns)
all_func_cols <- unique(all_func_cols)
all_func_cols


all_marker_cols <- c(
  "CD11b", "CD123", "CD133", "CD14", "CD141", "CD163", "CD20", "CD206", "CD208", 
  "CD209", "CD3", "CD31", "CD38", "CD4", "CD40", "CD45", "CD47", "CD68", "CD8", 
  "CD86", "Calprotectin", "Chym_Tryp", "FoxP3", "GFAP", "GLUT1", "HLA1", "HLADR", 
  "ICOS", "IDO1", "Ki67", "LAG3", "NeuN", "Olig2", "PD1", "PDL1", "TIM3", "TMEM119", 
  "Tox", "ApoE", "B7H3","CD70", "EGFR", "EGFRvIII", "EphA2","GFAP", "GM2_GD2", 
  "GPC2", "H3K27M", "H3K27me3", "HER2", "IDH1_R132H", 
  "NG2", "VISTA","Arginase1"
)

immune_list <- c(
  "CD11b", "CD14", "CD141", "CD163", "CD206", 
  "CD209", "CD3","CD38", "CD4", "CD45", "CD68", "CD8", 
  "CD86", "Calprotectin", "Chym_Tryp", "FoxP3", "GLUT1", "HLA1", "HLADR", 
  "ICOS", "IDO1", "Ki67", "LAG3", "PD1", "PDL1", "TIM3", "TMEM119", 
  "Tox", "Arginase1")

tumor_list <- c("HLA1",
                "ApoE", "B7H3", "CD133", "EGFR", "EGFRvIII", "GFAP", "GM2_GD2", 
                "GPC2", "H3K27M", "H3K27me3", "HER2","IDH1_R132H", "Ki67", 
                "NG2", "Olig2", "VISTA")


cell_type_to_all_marker_dict <- list(
  Myeloid_CD14 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Tumor_cells = list(
    list(markers = tumor_list, source = "cell_table_tumor"),
    list(markers = c("Tox", "GLUT1", "PDL1"), source = "cell_table_immune")),
  APC = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Macrophage_CD68_CD163 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Tcell_CD8 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Unassigned = list(
    list(markers = c(), source = "cell_table_tumor")),
  Endothelial_cells = list(
    list(markers = c("CD31", "GLUT1", "CD38", "CD141","IDO1"), source = "cell_table_immune")),
  Tcell_FoxP3 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Microglia = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Neurons = list(
    list(markers = c("NeuN", "Ki67", "HLA1", "Tox"), source = "cell_table_immune"),
    list(markers = c("ApoE"), source = "cell_table_tumor")),
  Immune_unassigned = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Macrophage_CD68 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Myeloid_CD141 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Macrophage_CD206 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Myeloid_CD11b = list(
    list(markers = immune_list, source = "cell_table_immune")),
  DC_CD123 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Neutrophils = list(
    list(markers = c(), source = "cell_table_immune")),
  DC_Mac_CD209 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  DC_CD206 = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Bcells = list(
    list(markers = immune_list, source = "cell_table_immune")),
  Mast_Cells = list(
    list(markers = immune_list, source = "cell_table_immune"))
)

#### import cell tables ######

cell_table_all_merged_thresholded <- arrow::read_parquet(file.path(output_dir, "../tables/cell_table_all_merged_thresholded.parquet"))

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### PREP DATA FOR FEATURE GENERATION 
##### generate columns to identify co expression of tumor markers ######

#filter only the tumor panel
cell_table_tumor_thresholded <- cell_table_all_merged_thresholded %>% 
  filter(source == "cell_table_tumor")

# make sure its a data table
setDT(cell_table_tumor_thresholded)

tumor_func_columns_core_simple
# Generate combinations and update the data table
for (i in 2) {
  combos <- combn(tumor_func_columns_core_simple, i, simplify = FALSE)
  for (combo in combos) {
    cat("Processing combination:", paste(combo, collapse = ", "), "\n")
    cell_table_tumor_thresholded <- compare_columns_multi(combo, cell_table_tumor_thresholded)
  }
}

# Print column names after operation
cat("Column names after operation:\n")
print(names(cell_table_tumor_thresholded))

write_parquet(cell_table_tumor_thresholded, "../tables/cell_table_tumor_thresholded.parquet")

# Select columns starting from "B7H3_func_EGFR_func" to the end that end with "_func"
tumor_combination_columns <- colnames(cell_table_tumor_thresholded) %>%
  .[which(. == "B7H3_func_EGFR_func"):length(.)] %>%
  .[grepl("_func$", .)]

# Combine tumor_func_columns_core_simple and tumor_combination_columns
complete_tumor_columns <- c(tumor_func_columns_core_simple, tumor_combination_columns)
complete_tumor_columns
save(complete_tumor_columns, file = "../tables/complete_tumor_columns.RData")

##### generate columns to identify co expression of immune functional markers ######

#filter only the immune panel
cell_table_immune_thresholded <- cell_table_all_merged_thresholded %>% 
  filter(source == "cell_table_immune")

immune_func_columns
# make sure its a data table
setDT(cell_table_immune_thresholded)

counter <- 0

# Generate combinations and update the data table
for (i in 2) {
  combos <- combn(immune_func_columns, i, simplify = FALSE)
  for (combo in combos) {
    counter <- counter + 1
    cat("Processing combination", counter, ":", paste(combo, collapse = ", "), "\n")
    cell_table_immune_thresholded <- compare_columns_multi(combo, cell_table_immune_thresholded)
  }
}

# Print column names after operation
cat("Column names after operation:\n")
print(names(cell_table_immune_thresholded))

write_parquet(cell_table_immune_thresholded, "../tables/cell_table_immune_thresholded.parquet")

# Select columns starting from "TIM3_func_CD38_func" to the end that end with "_func"
immune_combination_columns <- colnames(cell_table_immune_thresholded) %>%
  .[which(. == "TIM3_func_CD38_func"):length(.)] %>%
  .[grepl("_func$", .)]

# Combine immune_func_columns and immune_combination_columns
complete_immune_func_columns <- c(immune_func_columns, immune_combination_columns)
complete_immune_func_columns
save(complete_immune_func_columns, file = "../tables/complete_immune_func_columns.RData")

##### get size of each fov to calculate cell density #####
# 
# base_dir <- "/Volumes/Extreme Pro/GBM/cohorts/20230324/tumor/image_data"
# # List all subdirectories
# subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
# # Initialize a data frame to store results
# fov_area_df_tumor <- data.frame(fov = character(), height = integer(), width = integer(), stringsAsFactors = FALSE)
# # Loop through each subdirectory
# for (subdir in subdirs) {
#   print(subdir)
#   # Construct the path to the TIFF file
#   tiff_file_path <- file.path(subdir, "Nuclear.tiff")
# 
#   # Check if the file exists
#   if (file.exists(tiff_file_path)) {
#     # Read TIFF file
#     img <- readTIFF(tiff_file_path, info = TRUE)
# 
#     # Get dimensions
#     height <- dim(img)[1]
#     width <- dim(img)[2]
# 
#     # Append to the result data frame
#     fov_area_df_tumor <- rbind(fov_area_df_tumor, data.frame(fov = basename(subdir), height = height, width = width))
#   }
# }
# # Add area_in_microns column based on the height
# fov_area_df_tumor$area_in_microns <- ifelse(fov_area_df_tumor$height == 2048, 640000,
#                                             ifelse(fov_area_df_tumor$height == 1024, 160000, NA))
# # Print the updated result
# print(fov_area_df_tumor)
# write_parquet(fov_area_df_tumor, "../tables/fov_area_df_tumor.parquet")
# 
# 
# base_dir <- "/Volumes/Extreme Pro/GBM/cohorts/20230324/immune/image_data"
# # List all subdirectories
# subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
# # Initialize a data frame to store results
# fov_area_df_immune <- data.frame(fov = character(), height = integer(), width = integer(), stringsAsFactors = FALSE)
# # Loop through each subdirectory
# for (subdir in subdirs) {
#   print(subdir)
#   # Construct the path to the TIFF file
#   tiff_file_path <- file.path(subdir, "Nuclear.tiff")
# 
#   # Check if the file exists
#   if (file.exists(tiff_file_path)) {
#     # Read TIFF file
#     img <- readTIFF(tiff_file_path, info = TRUE)
# 
#     # Get dimensions
#     height <- dim(img)[1]
#     width <- dim(img)[2]
# 
#     # Append to the result data frame
#     fov_area_df_immune <- rbind(fov_area_df_immune, data.frame(fov = basename(subdir), height = height, width = width))
#   }
# }
# # Add area_in_microns column based on the height
# fov_area_df_immune$area_in_microns <- ifelse(fov_area_df_immune$height == 2048, 640000,
#                                              ifelse(fov_area_df_immune$height == 1024, 160000, NA))
# # Print the updated result
# print(fov_area_df_immune)
# 
# write_parquet(fov_area_df_immune, "../tables/fov_area_df_immune.parquet")


#####  generate a data table where each row is a sample_id and we have total tumor or immune  cell counts ##### 

tumor_count_result <- calculate_all_cell_type_count(cell_table_tumor_thresholded, "cell_meta_cluster_final_broad", "Tumor", "sample_id")
tumor_count_result <- dplyr::rename(tumor_count_result, all_tumor_count = all_cell_type_count)
write_parquet(tumor_count_result, "../tables/tumor_count_result.parquet")

tumor_count_result_patient_id <- calculate_all_cell_type_count(cell_table_tumor_thresholded, "cell_meta_cluster_final_broad", "Tumor", "patient_id")
tumor_count_result_patient_id <- dplyr::rename(tumor_count_result_patient_id, all_tumor_count = all_cell_type_count)
write_parquet(tumor_count_result_patient_id, "../tables/tumor_count_result_patient_id.parquet")

immune_count_result <- calculate_all_cell_type_count(cell_table_immune_thresholded, "cell_meta_cluster_final_broad", "Immune", "sample_id")
immune_count_result <- dplyr::rename(immune_count_result, all_immune_count = all_cell_type_count)
write_parquet(immune_count_result, "../tables/immune_count_result.parquet")

all_cell_count_immune_FOV_result <- calculate_all_cell_type_count(cell_table_immune_thresholded, "cell_meta_cluster_all", "all", "sample_id")
all_cell_count_immune_FOV_result <- dplyr::rename(all_cell_count_immune_FOV_result, all_cell_count_FOV_immune = all_cell_type_count)
write_parquet(all_cell_count_immune_FOV_result, "../tables/all_cell_count_immune_FOV_result.parquet")

all_cell_count_tumor_FOV_result <- calculate_all_cell_type_count(cell_table_tumor_thresholded, "cell_meta_cluster_all", "All", "sample_id")
all_cell_count_tumor_FOV_result <- dplyr::rename(all_cell_count_tumor_FOV_result, all_cell_count_tumor_FOV = all_cell_type_count)
write_parquet(all_cell_count_tumor_FOV_result, "../tables/all_cell_count_tumor_FOV_result.parquet")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### FEATURE GENERATION 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### marker intensity per sample ######

# Empty list to store results
results_list <- list()

# Loop through each cell type in the dictionary
for (cell_type in names(cell_type_to_all_marker_dict)) {
  # Retrieve list of configurations for the current cell type
  configurations <- cell_type_to_all_marker_dict[[cell_type]]
  
  for (config in configurations) {
    markers <- config$markers
    source_filter <- config$source
    print(config)
    # Filter data for the current cell type, source, and exclude 100.199 from markers
    # Note: Assuming `cell_table_all_merged_thresholded` already has a `source` column to filter on
    summary_data <- cell_table_all_merged_thresholded %>%
      filter(cell_meta_cluster_final == cell_type, source == source_filter) %>%
      group_by(sample_id, source, patient_id,cell_meta_cluster_final) %>%
      summarise(across(all_of(markers), ~mean(.x[.x != 0 & .x != 100.199], na.rm = TRUE)), .groups = 'drop')

    # Store summarised data in the list
    results_list[[paste(cell_type, source_filter, sep = "_")]] <- summary_data
  }
}

# Combine all results into a single data frame
final_result <- bind_rows(results_list)

# Create a function to apply the marker mapping
apply_marker_mapping <- function(df, dictionary) {
  marker_columns <- names(df)[sapply(df, is.numeric)]
  
  df <- df %>%
    rowwise() %>%
    mutate(across(
      all_of(marker_columns),
      ~ if_else(cur_column() %in% dictionary[[cell_meta_cluster_final]], .x, NA_real_)
    ))
  return(df)
}

# Apply the function to final_result
final_result_mapped <- apply_marker_mapping(final_result, combined_cell_marker_dictionary)

# Print the first few rows of the modified dataframe
print(head(final_result_mapped))


final_result_mapped[sapply(final_result_mapped, is.nan)] <- 0

# generate long format feature table
cell_table_marker_intensity <- final_result_mapped %>%
  mutate(feature_type = "Protein_Intensity",
         feature_source = "cell_meta_cluster_final",
         bio_feature_type = "protein_intensity")


all_cols_used <- colnames(cell_table_marker_intensity)[which(colnames(cell_table_marker_intensity) == "CD11b"):which(colnames(cell_table_marker_intensity) == "NeuN")]
all_cols_used <- c("CD11b", "CD14", "CD141", "CD163", "CD206", "CD209", "CD3", "CD38", "CD4",  "CD45", 
                   "CD68", "CD8", "CD86", "Calprotectin", "Chym_Tryp", "FoxP3", "GLUT1", "HLA1", "HLADR", "ICOS", "IDO1", "Ki67", 
                  "LAG3", "PD1", "PDL1", "TIM3", "TMEM119", "Tox", "ApoE", "B7H3", "CD133", "EGFR", "EGFRvIII", "GFAP", "GM2_GD2", "GPC2",
                  "H3K27M", "H3K27me3", "HER2", "IDH1_R132H", "NG2", "Olig2", "VISTA", "CD31", "NeuN", "Arginase1")
all_cols_used

melted_cell_table_marker_intensity <- melt(data = cell_table_marker_intensity,
                                               id.vars = c('sample_id','patient_id','source','cell_meta_cluster_final','feature_source','feature_type','bio_feature_type'),
                                               measure.vars = all_cols_used,
                                               variable.name = 'feature_variable',
                                               value.name = 'feature_value')



melted_cell_table_marker_intensity <- na.omit(melted_cell_table_marker_intensity, cols = "feature_value")


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_marker_intensity$Feature_Class <- "Non_Spatial"

# Update Feature_Type based on the feature_variable
melted_cell_table_marker_intensity$feature_type <- with(melted_cell_table_marker_intensity, 
                                                        ifelse(feature_variable %in% c("CD11b", "CD14", "CD141", "CD163", "CD206", "CD209", 
                                                                                       "CD3", "CD38", "CD4", "CD45", "CD68", "CD8", 
                                                                                       "Calprotectin", "Chym_Tryp", "FoxP3", "HLADR", 
                                                                                       "TMEM119", "GLUT1", "CD31", "NeuN", "HLA1", "Olig2", 
                                                                                       "GFAP", "CD133"), 
                                                               "Phenotype_marker_intensity",
                                                               ifelse(feature_variable %in% c("ICOS", "IDO1", "Ki67", "CD86", "H3K27me3", 
                                                                                              "LAG3", "PD1", "PDL1", "TIM3", "Tox", "ApoE","Arginase1"), 
                                                                      "Functional_marker_intensity",
                                                                      ifelse(feature_variable %in% c("B7H3", "EGFR", "EGFRvIII", "GM2_GD2", "GPC2", 
                                                                                                     "H3K27M", "HER2", "IDH1_R132H", "NG2", "VISTA"), 
                                                                             "Tumor_Antigens_Intensity", 
                                                                             feature_type))))  # Keep the original Feature_Type if none of the conditions are met

melted_cell_table_marker_intensity$Broad_Feature_Type <- "Protein"

write_parquet(melted_cell_table_marker_intensity, "../tables/melted_cell_table_marker_intensity.parquet")

##### marker intensity median ----

# Empty list to store results
results_list <- list()

# Loop through each cell type in the dictionary
for (cell_type in names(cell_type_to_all_marker_dict)) {
  # Retrieve list of configurations for the current cell type
  configurations <- cell_type_to_all_marker_dict[[cell_type]]
  
  for (config in configurations) {
    markers <- config$markers
    source_filter <- config$source
    print(config)
    # Filter data for the current cell type, source, and exclude 100.199 from markers
    # Note: Assuming `cell_table_all_merged_thresholded` already has a `source` column to filter on
    summary_data <- cell_table_all_merged_thresholded %>%
      filter(cell_meta_cluster_final == cell_type, source == source_filter) %>%
      group_by(sample_id, source, patient_id,cell_meta_cluster_final) %>%
      summarise(across(all_of(markers), ~median(.x[.x != 0 & .x != 100.199], na.rm = TRUE)), .groups = 'drop')
    
    # Store summarised data in the list
    results_list[[paste(cell_type, source_filter, sep = "_")]] <- summary_data
  }
}

# Combine all results into a single data frame
final_result <- bind_rows(results_list)

# Create a function to apply the marker mapping
apply_marker_mapping <- function(df, dictionary) {
  marker_columns <- names(df)[sapply(df, is.numeric)]
  
  df <- df %>%
    rowwise() %>%
    mutate(across(
      all_of(marker_columns),
      ~ if_else(cur_column() %in% dictionary[[cell_meta_cluster_final]], .x, NA_real_)
    ))
  return(df)
}

# Apply the function to final_result
final_result_mapped <- apply_marker_mapping(final_result, combined_cell_marker_dictionary)

# Print the first few rows of the modified dataframe
print(head(final_result_mapped))

final_result_mapped[sapply(final_result_mapped, is.nan)] <- 0

# generate long format feature table
cell_table_marker_intensity_median <- final_result_mapped %>%
  mutate(feature_type = "Protein_Intensity",
         feature_source = "cell_meta_cluster_final",
         bio_feature_type = "protein_intensity")

# Define the columns used
all_cols_used <- c("CD11b", "CD14", "CD141", "CD163", "CD206", "CD209", "CD3", "CD38", "CD4",  "CD45", 
                   "CD68", "CD8", "CD86", "Calprotectin", "Chym_Tryp", "FoxP3", "GLUT1", "HLA1", "HLADR", "ICOS", "IDO1", "Ki67", 
                   "LAG3", "PD1", "PDL1", "TIM3", "TMEM119", "Tox", "ApoE", "B7H3", "CD133", "EGFR", "EGFRvIII", "GFAP", "GM2_GD2", "GPC2",
                   "H3K27M", "H3K27me3", "HER2", "IDH1_R132H", "NG2", "Olig2", "VISTA", "CD31", "NeuN","Arginase1")

# Melt the data
melted_cell_table_marker_intensity_median <- melt(data = cell_table_marker_intensity_median,
                                           id.vars = c('sample_id','patient_id','source','cell_meta_cluster_final','feature_source','feature_type','bio_feature_type'),
                                           measure.vars = all_cols_used,
                                           variable.name = 'feature_variable',
                                           value.name = 'feature_value')

# Remove rows with NA values in the 'feature_value' column
melted_cell_table_marker_intensity_median <- na.omit(melted_cell_table_marker_intensity_median, cols = "feature_value")

# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_marker_intensity_median$Feature_Class <- "Non_Spatial"

# Update Feature_Type based on the feature_variable and append '_median'
melted_cell_table_marker_intensity_median$feature_type <- with(melted_cell_table_marker_intensity_median, 
                                                        ifelse(feature_variable %in% c("CD11b", "CD14", "CD141", "CD163", "CD206", "CD209", 
                                                                                       "CD3", "CD38", "CD4", "CD45", "CD68", "CD8", 
                                                                                       "Calprotectin", "Chym_Tryp", "FoxP3", "HLADR", 
                                                                                       "TMEM119", "GLUT1", "CD31", "NeuN", "HLA1", "Olig2", 
                                                                                       "GFAP", "CD133"), 
                                                               "Phenotype_marker_intensity_median",
                                                               ifelse(feature_variable %in% c("ICOS", "IDO1", "Ki67", "CD86", "H3K27me3", 
                                                                                              "LAG3", "PD1", "PDL1", "TIM3", "Tox", "ApoE","Arginase1"), 
                                                                      "Functional_marker_intensity_median",
                                                                      ifelse(feature_variable %in% c("B7H3", "EGFR", "EGFRvIII", "GM2_GD2", "GPC2", 
                                                                                                     "H3K27M", "HER2", "IDH1_R132H", "NG2", "VISTA"), 
                                                                             "Tumor_Antigens_Intensity_median", 
                                                                             feature_type))))  # Keep the original Feature_Type if none of the conditions are met

# Set broad feature type
melted_cell_table_marker_intensity_median$Broad_Feature_Type <- "Protein"

# Save the modified melted table to a parquet file
write_parquet(melted_cell_table_marker_intensity_median, "../tables/melted_cell_table_marker_intensity_median.parquet")

##### tumor antigen intensity segments #####

cell_table_used <- copy(cell_table_tumor_thresholded)  # Preserving the original dataframe
# Filter the dataframe to only include rows where cell_meta_cluster_broad is "Tumor"
cell_table_used <- cell_table_used[cell_table_used$cell_meta_cluster_final_broad == "Tumor", ]

# Assuming 'tumor_columns_core_simple' is already defined
quantile_counts_list <- list()

for (column_name in tumor_columns_core_simple) {
  
  # Filtering out 100.199 and 0 values for calculations
  filtered_values <- cell_table_used[[column_name]][cell_table_used[[column_name]] != 100.199 & cell_table_used[[column_name]] != 0]
  
  # Rank values and create equal bins
  ranked_values <- rank(-filtered_values, ties.method = "min")  # Rank descending
  n_bins <- 3
  bin_width <- ceiling(max(ranked_values) / n_bins)
  bins <- cut(ranked_values, breaks = seq(0, max(ranked_values) + bin_width, by = bin_width), include.lowest = TRUE, labels = FALSE)
  
  # Assign ranks back to the original dataframe
  cell_table_used[[column_name]][cell_table_used[[column_name]] != 100.199 & cell_table_used[[column_name]] != 0] <- bins
  
  # Create a column for binned segments
  segment_column_name <- paste0(column_name, "_segment")
  cell_table_used[[segment_column_name]] <- NA_integer_
  cell_table_used[[segment_column_name]][cell_table_used[[column_name]] != 100.199 & cell_table_used[[column_name]] != 0] <- bins
  
  # Calculate counts for each bin
  quantile_counts <- cell_table_used %>%
    filter(!!sym(column_name) != 100.199 & !!sym(column_name) != 0) %>%
    group_by(!!sym(segment_column_name)) %>%
    summarise(count = n())
  
  # Store counts for each column in a list
  quantile_counts_list[[column_name]] <- quantile_counts
}

# Now, 'cell_table_used' includes binned information for each column as well as their count statistics in 'quantile_counts_list'
quantile_counts_list


# Identify columns ending with '_segment'
segment_columns <- grep("_segment$", names(cell_table_used), value = TRUE)

# Generate a new dataframe with the selected columns
tumor_antigen_segment_df <- cell_table_used %>%
  dplyr::select(label, patient_id, sample_id, cell_meta_cluster_final, cell_meta_cluster_final_broad, all_of(segment_columns), all_of(tumor_func_columns_core_simple))

write_parquet(tumor_antigen_segment_df, "../tables/tumor_antigen_segment_df.parquet")

# Melt the dataframe
melted_tumor_antigen_segment_df <- tumor_antigen_segment_df %>%
  pivot_longer(cols = all_of(segment_columns), names_to = "segment", values_to = "value")

library(dplyr)
library(tidyr)

melted_tumor_antigen_segment_df <- melted_tumor_antigen_segment_df %>%
  group_by(sample_id, patient_id, segment) %>%
  summarise(
    count_1 = sum(value == 1, na.rm = TRUE),
    count_2 = sum(value == 2, na.rm = TRUE),
    count_3 = sum(value == 3, na.rm = TRUE),
    total = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    percent_1 = count_1 / total * 100,
    percent_2 = count_2 / total * 100,
    percent_3 = count_3 / total * 100
  ) %>%
  ungroup() %>%
  mutate(
    total_percent = percent_1 + percent_2 + percent_3,
    percent_1 = percent_1 / total_percent,
    percent_2 = percent_2 / total_percent,
    percent_3 = percent_3 / total_percent
  ) %>%
  dplyr::select(-total_percent)

melted_tumor_antigen_segment_df <- melted_tumor_antigen_segment_df %>%
  pivot_longer(cols = starts_with("percent_"), names_to = "percent_type", values_to = "percent_value") %>%
  mutate(segment_percent = case_when(
    percent_type == "percent_1" ~ paste0(segment, "_high_prop"),
    percent_type == "percent_2" ~ paste0(segment, "_med_prop"),
    percent_type == "percent_3" ~ paste0(segment, "_low_prop")
  )) %>%
  dplyr::select(sample_id, patient_id, segment_percent, percent_value) %>%
  pivot_longer(cols = percent_value, names_to = "value_type", values_to = "value") %>%
  mutate(segment = segment_percent) %>%
  dplyr::select(sample_id, patient_id, segment, value) %>%
  mutate(feature_type = "Tumor_Antigens_Intensity_Segments", feature_source = "whole_sample")

colnames(melted_tumor_antigen_segment_df)[colnames(melted_tumor_antigen_segment_df) == "segment"] <- "feature_variable"
colnames(melted_tumor_antigen_segment_df)[colnames(melted_tumor_antigen_segment_df) == "value"] <- "feature_value"

melted_tumor_antigen_segment_df <- melted_tumor_antigen_segment_df %>%
  mutate(
    feature_type = ifelse(grepl("prop", feature_variable), "Tumor_Antigens_Intensity_Segments", "Tumor_Antigens_Intensity_Segments"),
    feature_source = "whole_sample"
  )

# View the reshaped dataframe
print(melted_tumor_antigen_segment_df)

melted_tumor_antigen_segment_df <- na.omit(melted_tumor_antigen_segment_df, cols = "feature_value")
melted_tumor_antigen_segment_df$bio_feature_type <- "protein_intensity"

# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_tumor_antigen_segment_df$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_tumor_antigen_segment_df$Broad_Feature_Type <- "Protein"

write_parquet(melted_tumor_antigen_segment_df, "../tables/melted_tumor_antigen_segment_df.parquet")

##### ratio of tumor cell types per sample ######

positive_count_result_tumor <- calculate_positive_counts(
  cell_table_tumor_thresholded, 
  complete_tumor_columns, 
  summary_by = "sample_id", 
  cell_meta_cluster_column = "cell_meta_cluster_final_broad", 
  cell_meta_cluster_value = "Tumor"
)


all_complete_tumor_columns_counts <- paste0(complete_tumor_columns, "_counts")

write_parquet(positive_count_result_tumor, "../tables/positive_count_result_tumor.parquet")

# # Assuming `positive_count_result_tumor` is already created
# melted_positive_count_result_tumor <- melt(
#   positive_count_result_tumor,
#   id.vars = "sample_id",
#   measure.vars = all_complete_tumor_columns_counts,
#   variable.name = "feature_variable",
#   value.name = "feature_value"
# )
# 
# # generate long format feature table
# melted_positive_count_result_tumor <- melted_positive_count_result_tumor %>%
#   mutate(feature_type = "tumor_cell_features_counts",
#          feature_source = "whole_sample")
# 
# 
# # View the result
# head(melted_positive_count_result_tumor)
# 
# write_parquet(melted_positive_count_result_tumor, "../tables/melted_positive_count_result_tumor.parquet")

# Merge tumor_count and non_empty_counts based on sample_id
cell_prop_tumor_thresholded_sample_id <- merge(positive_count_result_tumor, tumor_count_result, by = "sample_id")
cell_prop_tumor_thresholded_sample_id <- merge(cell_prop_tumor_thresholded_sample_id, all_cell_count_tumor_FOV_result, by = "sample_id")

# get the names of columns that will be the numerator in our calculation
#proportion_columns <- setdiff(names(cell_prop_tumor_thresholded_sample_id), c("sample_id", "all_tumor_count","all_cell_count_tumor_FOV"))
proportion_columns <- setdiff(names(cell_prop_tumor_thresholded_sample_id), c("sample_id","all_cell_count_tumor_FOV"))
proportion_columns


# Filter elements with more than one '_func_' substring
filtered_proportion_columns <- proportion_columns[grepl("(.*_func_.*){1,}", proportion_columns) | !grepl("_func_", proportion_columns)]

print(filtered_proportion_columns)

# these will be the denominator or the tumor cell types we are interested in finding the ratios with.
tumor_func_columns_core_simple_count <- paste0(tumor_func_columns_core_simple, "_counts")
tumor_func_columns_core_simple_count <- c(tumor_func_columns_core_simple_count, "all_tumor_count","all_cell_count_tumor_FOV")
tumor_func_columns_core_simple_count


for (num_col in filtered_proportion_columns) {
  print(num_col)
  for (denom_col in tumor_func_columns_core_simple_count) {
    # Skip if the numerator and denominator columns are the same
    if (num_col == denom_col) {
      next
    }
    
    # Extract the components from the numerator column
    num_col_components <- unlist(strsplit(sub("_counts$", "", num_col), "_func_"))
    
    if (denom_col %in% c("all_tumor_count", "all_cell_count_tumor_FOV") || any(sapply(num_col_components, function(component) grepl(component, denom_col)))) {
      new_col_name = paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "prop", sep = "_")
      # Initialize the new column with NA or appropriate default value
      cell_prop_tumor_thresholded_sample_id[[new_col_name]] <- NA
      for (i in 1:nrow(cell_prop_tumor_thresholded_sample_id)) {
        if (cell_prop_tumor_thresholded_sample_id[i, denom_col] > 20) {
          cell_prop_tumor_thresholded_sample_id[i, new_col_name] <- cell_prop_tumor_thresholded_sample_id[i, num_col] / cell_prop_tumor_thresholded_sample_id[i, denom_col]
        }
      }
    }
  }
}
tumor_func_columns_count <- paste0(tumor_func_columns_core_simple, "_counts")
tumor_func_columns_count <- c(tumor_func_columns_count, "all_tumor_count")
tumor_func_columns_count

tumor_func_columns_except_all = setdiff(tumor_func_columns_count, "all_tumor_count")
tumor_func_columns_except_all

for (num_index in 1:length(tumor_func_columns_except_all)) {
  num_col <- tumor_func_columns_except_all[num_index]
  if (!num_col %in% names(cell_prop_tumor_thresholded_sample_id)) next
  
  for (denom_index in (num_index + 1):length(tumor_func_columns_except_all)) {
    denom_col <- tumor_func_columns_except_all[denom_index]
    if (!denom_col %in% names(cell_prop_tumor_thresholded_sample_id)) next
    
    # Proceed with the assumption that num_col and denom_col are valid column names
    new_col_name_A = paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", num_col), "plus", sub("_counts$", "", denom_col), "prop", sep = "_")
    new_col_name_B = paste(sub("_counts$", "", denom_col), "over", sub("_counts$", "", num_col), "plus", sub("_counts$", "", denom_col), "prop", sep = "_")
    
    
    # Initialize new columns with NA
    cell_prop_tumor_thresholded_sample_id[[new_col_name_A]] <- NA
    cell_prop_tumor_thresholded_sample_id[[new_col_name_B]] <- NA
    
    for (i in 1:nrow(cell_prop_tumor_thresholded_sample_id)) {
      sum_of_columns <- cell_prop_tumor_thresholded_sample_id[i, num_col] + cell_prop_tumor_thresholded_sample_id[i, denom_col]
      if (length(sum_of_columns) > 0 && sum_of_columns > 20) {
        cell_prop_tumor_thresholded_sample_id[i, new_col_name_A] <- cell_prop_tumor_thresholded_sample_id[i, num_col] / sum_of_columns
        cell_prop_tumor_thresholded_sample_id[i, new_col_name_B] <- cell_prop_tumor_thresholded_sample_id[i, denom_col] / sum_of_columns
      }
    }
  }
}


#grep("EGFR_func_over_B7H3_func_counts_plus_EGFR_func_counts_prop", names(cell_prop_tumor_thresholded_sample_id), value = TRUE)


# 
# # Call the function with your specific data and parameters
# adjusted_proportions_df <- calculate_adjusted_proportions(
#   data = cell_table_tumor_thresholded,
#   marker_columns = tumor_func_columns_core_simple, # The list of marker columns except 'all_tumor_count'
#   cell_meta_cluster_col = "cell_meta_cluster_final_broad", # The metadata column for cell type
#   cluster_value = "Tumor", # The value indicating a tumor cell
#   summary_by = "sample_id" # The column by which to summarize
# )
# 
# 
# 
# # Call the function with your specific data and parameters
# list_counts_dfs <- calculate_adjusted_counts(
#   data = cell_table_tumor_thresholded,
#   marker_columns = tumor_func_columns_core_simple, # The list of marker columns
#   cell_meta_cluster_col = "cell_meta_cluster_final_broad", # The metadata column for cell type
#   cluster_value = "Tumor", # The value indicating a tumor cell
#   summary_by = "sample_id" # The column by which to summarize
# )
# 
# 
# # Initialize the combined data frame with sample_id
# combined_counts_df <- data.frame(sample_id = list_counts_dfs[[1]]$sample_id)
# 
# # Loop through each data frame in the list
# for (marker_name in names(list_counts_dfs)) {
#   # Get the data frame from the list
#   df <- list_counts_dfs[[marker_name]]
#   
#   # Rename all columns except sample_id to include the specific marker name as suffix
#   colnames(df)[-1] <- paste(colnames(df)[-1], "minus", marker_name, sep = "_")
#   
#   # Merge the data frame with the combined data frame
#   if ("sample_id" %in% names(combined_counts_df)) {
#     combined_counts_df <- merge(combined_counts_df, df, by = "sample_id", all = TRUE)
#   } else {
#     combined_counts_df <- df  # If the combined_counts_df is empty, initialize it with the first df
#   }
# }
# 
# # Replace any NAs with zero
# combined_counts_df[is.na(combined_counts_df)] <- 0
# 
# 



####

# Identify the columns ending with '_counts' and '_prop'
tumor_counts_cols <- grep("_counts$", names(cell_prop_tumor_thresholded_sample_id), value = TRUE)
tumor_prop_cols <- grep("_prop$", names(cell_prop_tumor_thresholded_sample_id), value = TRUE)

# Add 'sample_id' column to both lists
tumor_counts_cols_sample_id <- c("sample_id", tumor_counts_cols)
tumor_prop_cols_sample_id <- c("sample_id", tumor_prop_cols)

# Create two subsets, one has the counts and other has the ratios. 
tumor_cell_type_count <- cell_prop_tumor_thresholded_sample_id[, tumor_counts_cols_sample_id]
write_parquet(tumor_cell_type_count, "../tables/tumor_cell_type_count.parquet")

tumor_cell_ratios <- cell_prop_tumor_thresholded_sample_id[, tumor_prop_cols_sample_id]



# # Perform an anti-join to find sample_ids in adjusted_proportions_df not in tumor_cell_ratios
# sample_ids_not_found <- anti_join(adjusted_proportions_df, tumor_cell_ratios, by = "sample_id")
# tumor_cell_ratios_combined <- merge(tumor_cell_ratios,adjusted_proportions_df, by = "sample_id", all = FALSE)
# 
# tumor_prop_cols_adjusted <- grep("_prop$", names(adjusted_proportions_df), value = TRUE)
# tumor_prop_cols_combined <- c(tumor_prop_cols_adjusted, tumor_prop_cols)
# tumor_prop_cols_combined
# 
# # Print out the sample_ids
# print(sample_ids_not_found$sample_id)


# generate long format feature table
tumor_cell_ratios <- tumor_cell_ratios %>%
  mutate(feature_type = "tumor_cell_features",
         feature_source = "whole_sample")



tumor_cell_ratios <- combine_data_metadata(tumor_cell_ratios, metadata, "sample_id")

setDT(tumor_cell_ratios)

write_parquet(tumor_cell_ratios, "../tables/tumor_cell_ratios_combined.parquet")


melted_cell_table_tumor_cell_ratios <- melt(data = tumor_cell_ratios,
                                           id.vars = c('sample_id','patient_id','feature_source','feature_type'),
                                           measure.vars = tumor_prop_cols,
                                           variable.name = 'feature_variable',
                                           value.name = 'feature_value')


melted_cell_table_tumor_cell_ratios <- melted_cell_table_tumor_cell_ratios %>%
  mutate(feature_type = case_when(
    str_count(feature_variable, "func") == 3 & str_detect(feature_variable, "(func.*func.*_over_).*func") ~ "Cell_Ratios",
    str_detect(feature_variable, "plus") ~ "Cell_Ratios",
    str_detect(feature_variable, "over_all_tumor_count") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_cell_count_tumor_FOV_prop") ~ "Cell_Abundance",
    TRUE ~ "not_relevant"
  ))

melted_cell_table_tumor_cell_ratios <- melted_cell_table_tumor_cell_ratios %>%
  mutate(bio_feature_type = case_when(
    str_count(feature_variable, "func") == 3 & str_detect(feature_variable, "(func.*func.*_over_).*func") ~ "Tumor_cells",
    str_detect(feature_variable, "plus") ~ "Tumor_cells",
    str_detect(feature_variable, "over_all_tumor_count") ~ "Relative_to_all_tumor_cells",
    str_detect(feature_variable, "all_cell_count_tumor_FOV_prop") ~ "Relative_to_all_cells",
    TRUE ~ "not_relevant"
  ))

# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_tumor_cell_ratios$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_tumor_cell_ratios$Broad_Feature_Type <- "Cells"


write_parquet(melted_cell_table_tumor_cell_ratios, "../tables/melted_cell_table_tumor_cell_ratios.parquet")
#melted_cell_table_tumor_cell_ratios <- read_parquet("melted_cell_table_tumor_cell_ratios.parquet")


##### calculate tumor antigen coverage  #####

find_top_seven <- function(data, marker_columns, summary_by, cell_meta_cluster_col, cluster_value) {
  data <- data[data[[cell_meta_cluster_col]] == cluster_value, ]
  data <- data[!is.na(data[[summary_by]]), ]
  data[[summary_by]] <- as.character(data[[summary_by]])
  
  result_df <- data.frame(summary = unique(data[[summary_by]]))
  
  for (i in 1:7) {
    result_df[paste0("marker_", i, "_name")] <- NA
    result_df[paste0("marker_", i, "_count")] <- NA
    result_df[paste0("marker_", i, "_prop")] <- NA
  }
  
  result_df$all_tumor_count <- NA
  result_df$no_marker_cells <- NA
  result_df$true_remaining_prop <- NA
  
  for (summary_value in unique(data[[summary_by]])) {
    summary_char <- as.character(summary_value)
    sample_data <- data[data[[summary_by]] == summary_char, ]
    
    total_tumor_cells <- sum(sample_data[[cell_meta_cluster_col]] == cluster_value, na.rm = TRUE)
    result_df[result_df$summary == summary_char, "all_tumor_count"] <- total_tumor_cells
    
    remaining_data <- sample_data
    for (i in 1:7) {
      marker_counts <- sapply(marker_columns, function(marker) sum(remaining_data[[marker]] == 1, na.rm = TRUE))
      most_marker <- names(which.max(marker_counts))
      most_count <- max(marker_counts)
      most_prop <- ifelse(total_tumor_cells > 0, most_count / total_tumor_cells, 0)
      
      result_df[result_df$summary == summary_char, paste0("marker_", i, "_name")] <- most_marker
      result_df[result_df$summary == summary_char, paste0("marker_", i, "_count")] <- most_count
      result_df[result_df$summary == summary_char, paste0("marker_", i, "_prop")] <- most_prop
      
      remaining_data <- remaining_data[remaining_data[[most_marker]] != 1, ]
    }
    
    no_marker_cells <- sum(remaining_data[[cell_meta_cluster_col]] == cluster_value)
    true_remaining_prop <- ifelse(total_tumor_cells > 0, no_marker_cells / total_tumor_cells, 0)
    
    result_df[result_df$summary == summary_char, "no_marker_cells"] <- no_marker_cells
    result_df[result_df$summary == summary_char, "true_remaining_prop"] <- true_remaining_prop
  }
  
  ordered_cols <- c("summary", "all_tumor_count", paste0("marker_", 1:7, "_name"), paste0("marker_", 1:7, "_count"), paste0("marker_", 1:7, "_prop"), "no_marker_cells", "true_remaining_prop")
  result_df <- result_df[ordered_cols]
  
  return(result_df)
}


tumor_antigen_segment_selected <- tumor_antigen_segment_df %>%
  dplyr::select(label, sample_id, ends_with("_segment"))

# Merge the selected columns with cell_table_tumor_thresholded
combined_df <- merge(tumor_antigen_segment_selected, cell_table_tumor_thresholded, by = c("label", "sample_id"))

# filtered_df <- combined_df %>%
#   filter(if_any(ends_with("_segment"), ~ . %in% c(2, 3)))

# FOR FIGURES
# 
# # Example parameters - adjust these according to your specific dataset
# data <- combined_df   # Your input data frame
# marker_columns <- tumor_func_columns_core_simple   # List of marker columns
# summary_by <- "patient_id"   # Column used for grouping
# cell_meta_cluster_col <- "cell_meta_cluster_final_broad"   # Column indicating cell clusters
# cluster_value <- "Tumor"   # Specific cluster to analyze
# 
# # Call the function with the specified parameters
# top_seven_df <- find_top_seven(
#   data = data,
#   marker_columns = marker_columns,
#   summary_by = summary_by,
#   cell_meta_cluster_col = cell_meta_cluster_col,
#   cluster_value = cluster_value
# )
# 
# # Rename the column
# top_seven_df <- top_seven_df %>%
#   dplyr::rename(patient_id = summary)

# FOR MASTER FEATURE table

# Example parameters - adjust these according to your specific dataset
data <- combined_df   # Your input data frame
marker_columns <- tumor_func_columns_core_simple   # List of marker columns
summary_by <- "sample_id"   # Column used for grouping
cell_meta_cluster_col <- "cell_meta_cluster_final_broad"   # Column indicating cell clusters
cluster_value <- "Tumor"   # Specific cluster to analyze

# Call the function with the specified parameters
top_seven_df <- find_top_seven(
  data = data,
  marker_columns = marker_columns,
  summary_by = summary_by,
  cell_meta_cluster_col = cell_meta_cluster_col,
  cluster_value = cluster_value
)

# Rename the column
top_seven_df <- top_seven_df %>%
  dplyr::rename(sample_id = summary)

top_seven_df <- combine_data_metadata(top_seven_df, metadata, "sample_id")


# Function to set duplicate marker names to NA
remove_duplicates <- function(row) {
  unique_markers <- c()
  for (i in 1:7) {
    marker_col <- paste0('marker_', i, '_name')
    if (row[[marker_col]] %in% unique_markers) {
      row[[marker_col]] <- NA
    } else {
      unique_markers <- c(unique_markers, row[[marker_col]])
    }
  }
  return(row)
}

# Apply the function to each row in the dataframe
top_seven_df <- as.data.frame(t(apply(top_seven_df, 1, remove_duplicates)))

# Select columns that start with 'marker' and end with 'count' or 'prop'
marker_columns <- grep("^marker.*(count|prop)$", names(top_seven_df), value = TRUE)

# Convert those columns to numeric
top_seven_df[marker_columns] <- lapply(top_seven_df[marker_columns], as.numeric)

# write_parquet(top_seven_df, "../tables/top_seven_df.parquet")


# Add the new columns
top_seven_df$feature_type <- 'Cell_Abundance'
top_seven_df$feature_source <- 'whole_sample'

# Initialize an empty dataframe to store the final results
melted_tumor_antigen_coverage <- data.frame()
id_vars <- c('sample_id', 'patient_id', 'feature_type', 'feature_source')

# Iterate over marker columns and bind the melted data
for (i in 1:7) {
  name_col <- paste0('marker_', i, '_name')
  print(name_col)
  melted_df <- melt(top_seven_df,
                    id.vars = id_vars,
                    measure.vars = name_col,
                    variable.name = 'temp_variable',
                    value.name = 'feature_variable')
  
  melted_df$feature_value <- i
  
  melted_tumor_antigen_coverage <- rbind(melted_tumor_antigen_coverage, melted_df)
}

# Output the final dataframe
melted_tumor_antigen_coverage

# Remove the temp_variable column
melted_tumor_antigen_coverage$temp_variable <- NULL

melted_tumor_antigen_coverage <- melted_tumor_antigen_coverage[!is.na(melted_tumor_antigen_coverage$feature_variable), ]
melted_tumor_antigen_coverage$feature_variable <- paste0(melted_tumor_antigen_coverage$feature_variable, "_coverage")

melted_tumor_antigen_coverage$bio_feature_type <- "Tumor_cell_coverage"

write_parquet(melted_tumor_antigen_coverage, "../tables/melted_tumor_antigen_coverage.parquet")

# Join the two dataframes by `patient_id`





merged_df <- top_seven_df %>%
  left_join(tumor_count_result %>% dplyr::select(sample_id, all_tumor_count), by = "sample_id", suffix = c("", "_new"))

merged_df <- merged_df %>%
  dplyr::select(-all_tumor_count) %>%  # Remove the old column
  dplyr::rename(all_tumor_count = all_tumor_count_new)  # Rename the new column

# Recalculate the _prop columns
merged_df <- merged_df %>%
  mutate(marker_1_prop = marker_1_count / all_tumor_count,
         marker_2_prop = marker_2_count / all_tumor_count,
         marker_3_prop = marker_3_count / all_tumor_count,
         marker_4_prop = marker_4_count / all_tumor_count,
         marker_5_prop = marker_5_count / all_tumor_count,
         marker_6_prop = marker_6_count / all_tumor_count,
         marker_7_prop = marker_7_count / all_tumor_count)

write_parquet(merged_df, "../tables/top_seven_df.parquet")


# Get the columns that end with _name and the true_remaining_prop column
measure_vars <- c(grep("_prop$", names(top_seven_df), value = TRUE))

# Melt the dataframe
melted_tumor_antigen_coverage_amounts <- melt(top_seven_df, 
                  id.vars = c('sample_id', 'patient_id','feature_source', 'feature_type'),
                  measure.vars = measure_vars,
                  variable.name = 'feature_variable',
                  value.name = 'feature_value')
melted_tumor_antigen_coverage_amounts$feature_value <- as.numeric(melted_tumor_antigen_coverage_amounts$feature_value)

melted_tumor_antigen_coverage_amounts$bio_feature_type <- "cell_abundance"

# Display the melted dataframe
print(melted_tumor_antigen_coverage_amounts)


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_tumor_antigen_coverage_amounts$Feature_Class <- "Non_Spatial"

melted_tumor_antigen_coverage_amounts$feature_type <- "Cell_Abundance"

# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_tumor_antigen_coverage_amounts$Broad_Feature_Type <- "Cells"


write_parquet(melted_tumor_antigen_coverage_amounts, "../tables/melted_tumor_antigen_coverage_amounts.parquet")

##### calculate most co expression of tumor antigens #####

# Assuming cell_table_tumor_threshold is your dataframe
# and tumor_func_columns_core_simple is a vector of column names



# Step 1: Filter for Tumor cells
tumor_cells <- cell_table_tumor_thresholded %>%
  filter(cell_meta_cluster_final_broad == "Tumor")

calculate_overlap <- function(df, markers, group_var) {
  overlaps <- list()
  total_combinations <- (length(markers) * (length(markers) - 1)) / 2
  pb <- progress_bar$new(total = total_combinations, format = "  calculating [:bar] :percent in :elapsed")
  
  for (i in 1:(length(markers) - 1)) {
    for (j in (i + 1):length(markers)) {
      marker1 <- markers[i]
      marker2 <- markers[j]
      overlap_count <- df %>%
        filter(!!sym(marker1) == 1 & !!sym(marker2) == 1) %>%
        group_by(!!sym(group_var)) %>%
        summarise(overlap = n(), .groups = 'drop')
      overlap_count <- overlap_count %>%
        mutate(marker1 = marker1, marker2 = marker2)
      overlaps[[paste(marker1, marker2, sep = "_")]] <- overlap_count

      pb$tick()
    }
  }
  bind_rows(overlaps)
}

tumor_antigen_coexp_counts <- calculate_overlap(tumor_cells, tumor_func_columns_core_simple, "sample_id")

tumor_antigen_coexp_counts

overlap_counts <- calculate_overlap(tumor_cells, tumor_func_columns_core_simple, "patient_id")

# Step 3: Find the top 2 markers with the most overlap for each patient
overlap_counts <- overlap_counts %>%
  group_by(patient_id) %>%
  slice_max(order_by = overlap, n = 1, with_ties = FALSE) %>%
  mutate(combination = paste(sort(c(marker1, marker2)), collapse = "_")) %>%
  ungroup()

overlap_counts

write_parquet(overlap_counts, "../tables/top_coexp_combinations.parquet")

# Assuming metadata dataframe is already loaded
tumor_antigen_coexp_counts <- tumor_antigen_coexp_counts %>%
  left_join(metadata %>% dplyr::select(sample_id, patient_id), by = "sample_id") %>%
  mutate(feature_variable = paste(marker1, marker2, sep = "_"))


# Add the new columns
tumor_antigen_coexp_counts <- tumor_antigen_coexp_counts %>%
  mutate(
    feature_type = 'Cell_Abundance',
    feature_source = 'whole_sample',
    bio_feature_type = "cell_abundance"
  )

# Melt the dataframe
melted_tumor_antigen_coexp_counts <- melt(
  data = tumor_antigen_coexp_counts,
  id.vars = c('sample_id', 'patient_id', 'feature_source', 'feature_type','bio_feature_type'),
  measure.vars = 'overlap',
  variable.name = 'feature_variable',
  value.name = 'feature_value'
)

melted_tumor_antigen_coexp_counts$feature_variable <- tumor_antigen_coexp_counts$feature_variable

# Display the resulting dataframe
head(melted_tumor_antigen_coexp_counts)

# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_tumor_antigen_coexp_counts$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_tumor_antigen_coexp_counts$Broad_Feature_Type <- "Cells"

write_parquet(melted_tumor_antigen_coexp_counts,"../tables/melted_tumor_antigen_coexp_counts.parquet")
  
# Step 3: Find the top 2 markers with the most overlap for each patient
tumor_antigen_coexp_top_comb <- tumor_antigen_coexp_counts %>%
    group_by(sample_id) %>%
    slice_max(order_by = overlap, n = 1, with_ties = FALSE) %>%
    mutate(combination = paste(sort(c(marker1, marker2)), collapse = "_")) %>%
    ungroup()

tumor_antigen_coexp_top_comb
  
write_parquet(tumor_antigen_coexp_top_comb, "../tables/tumor_antigen_coexp_top_comb.parquet")
  

##### ratio of immune cells (based on pixie cell phenotyping) to total immune cells and tumor cell types ######

positive_count_result_immune_cells <- cell_table_immune_thresholded %>%
  filter(source == "cell_table_immune") %>%
  count(sample_id, cell_meta_cluster_final) %>%
  spread(key = cell_meta_cluster_final, value = n, fill = 0)
colnames(positive_count_result_immune_cells)

colnames(positive_count_result_immune_cells)[-which(colnames(positive_count_result_immune_cells) == "sample_id")] <- 
  paste0(colnames(positive_count_result_immune_cells)[-which(colnames(positive_count_result_immune_cells) == "sample_id")], "_counts")
colnames(positive_count_result_immune_cells)

positive_count_result_immune_cells_columns <- colnames(positive_count_result_immune_cells)
positive_count_result_immune_cells_columns <- setdiff(colnames(positive_count_result_immune_cells), "sample_id")
positive_count_result_immune_cells_columns

write_parquet(positive_count_result_immune_cells, "../tables/positive_count_result_immune_cells.parquet")

# 
# # Assuming `positive_count_result_tumor` is already created
# melted_positive_count_result_immune <- melt(
#   positive_count_result_immune_cells,
#   id.vars = "sample_id",
#   measure.vars = positive_count_result_immune_cells_columns,
#   variable.name = "feature_variable",
#   value.name = "feature_value"
# )
# 
# # generate long format feature table
# melted_positive_count_result_immune <- melted_positive_count_result_immune %>%
#   mutate(feature_type = "immune_cell_features_counts",
#          feature_source = "whole_sample")
# 
# 
# # View the result
# head(melted_positive_count_result_immune)

# write_parquet(melted_positive_count_result_immune, "../tables/melted_positive_count_result_immune.parquet")
# melted_positive_count_result_immune <- read_parquet("melted_positive_count_result_immune.parquet")
# 
# 




# Merge tumor_count and non_empty_counts based on sample_id
positive_count_result_immune_cells_sample_id <- merge(positive_count_result_immune_cells, tumor_count_result, by = "sample_id")
positive_count_result_immune_cells_sample_id <- merge(positive_count_result_immune_cells_sample_id, immune_count_result, by = "sample_id")
positive_count_result_immune_cells_sample_id <- merge(positive_count_result_immune_cells_sample_id, all_cell_count_immune_FOV_result, by = "sample_id")

positive_count_result_immune_cells_sample_id

tumor_func_columns_core_simple
tumor_func_columns_count <- paste0(tumor_func_columns_core_simple, "_counts")
tumor_func_columns_count_sample_id <- c(tumor_func_columns_count, "sample_id")

positive_count_result_tumor_core <- positive_count_result_tumor %>% 
  dplyr::select(all_of(tumor_func_columns_count_sample_id))
positive_count_result_tumor_core

positive_count_result_immune_cells_sample_id <- merge(positive_count_result_tumor_core, positive_count_result_immune_cells_sample_id, by = "sample_id")
colnames(positive_count_result_immune_cells_sample_id)



proportion_columns <- c(unique(cell_table_immune_thresholded$cell_meta_cluster_final))
proportion_columns <- setdiff(proportion_columns, "Tumor_cells")
proportion_columns <- paste0(proportion_columns, "_counts")
proportion_columns <- c(proportion_columns)
proportion_columns 


denominator_columns <- c(proportion_columns,tumor_func_columns_count,"all_tumor_count","all_immune_count","all_cell_count_FOV_immune" )
denominator_columns
immune_cells_prop_sample_id <- copy(positive_count_result_immune_cells_sample_id)


denominator_columns <-  c("all_tumor_count","all_immune_count","all_cell_count_FOV_immune" )
# Initialize a copy of the original data frame for simple ratios
immune_cells_prop_sample_id_ratios <- copy(positive_count_result_immune_cells_sample_id)

# Loop for simple ratios
for (num_col in proportion_columns) {
  for (denom_col in denominator_columns) {
    if (num_col == denom_col) {
      next
    }
    new_col_name = paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "prop", sep = "_")
    immune_cells_prop_sample_id_ratios[[new_col_name]] <- immune_cells_prop_sample_id[[num_col]] / immune_cells_prop_sample_id[[denom_col]]
  }
}

# Initialize a copy of the original data frame for fractions
immune_cells_prop_sample_id_fractions <- copy(positive_count_result_immune_cells_sample_id)
proportion_columns
denominator_columns_without_all <- c(tumor_func_columns_count, proportion_columns)
denominator_columns_without_all


# Assuming proportion_columns and denominator_columns_without_all are sorted or you sort them beforehand
proportion_columns <- sort(proportion_columns)
denominator_columns_without_all <- sort(denominator_columns_without_all)
denominator_columns_without_all
for (num_col in proportion_columns) {
  print(num_col)
  for (denom_col in denominator_columns_without_all) {
    print(denom_col)
    if (num_col == denom_col) {
      next
    }
    # Check if the index of num_col is less than denom_col to ensure each pair is processed once
    if (which(proportion_columns == num_col) < which(denominator_columns_without_all == denom_col)) {
      new_col_name_A = paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", num_col), "plus", sub("_counts$", "", denom_col), "prop", sep = "_")
      new_col_name_B = paste(sub("_counts$", "", denom_col), "over", sub("_counts$", "", num_col), "plus", sub("_counts$", "", denom_col), "prop", sep = "_")
      
      immune_cells_prop_sample_id_fractions[[new_col_name_A]] <- immune_cells_prop_sample_id[[num_col]] / 
        (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
      immune_cells_prop_sample_id_fractions[[new_col_name_B]] <- immune_cells_prop_sample_id[[denom_col]] / 
        (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
    }
  }
}



proportion_columns <- sort(proportion_columns)
denominator_columns_without_all <- sort(denominator_columns_without_all)

immune_cells_prop_sample_id_fractions <- list() # Initialize if not already existing



for (num_col in proportion_columns) {
  num_index <- which(proportion_columns == num_col)
  for (denom_col in denominator_columns_without_all) {
    denom_index <- which(denominator_columns_without_all == denom_col)
    
    if (num_col == denom_col || num_index >= denom_index) {
      next
    }
    
    new_col_name <- paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "plus", sub("_counts$", "", num_col), "prop", sep = "_")
    
    immune_cells_prop_sample_id_fractions[[new_col_name]] <- immune_cells_prop_sample_id[[num_col]] / 
      (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
  }
}


# Assuming immune_cells_prop_sample_id_fractions is already a defined data frame
proportion_columns <- sort(proportion_columns)
denominator_columns_without_all <- sort(denominator_columns_without_all)

for (num_col in proportion_columns) {
  num_index <- which(proportion_columns == num_col)
  for (denom_col in denominator_columns_without_all) {
    denom_index <- which(denominator_columns_without_all == denom_col)
    
    if (num_col == denom_col || num_index >= denom_index) {
      next
    }
    
    new_col_name <- paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "plus", sub("_counts$", "", num_col), "prop", sep = "_")
    
    # Append new ratio calculations directly into the data frame
    immune_cells_prop_sample_id_fractions[[new_col_name]] <- immune_cells_prop_sample_id[[num_col]] / 
      (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
  }
}


# Check and print out all generated column names
print(names(immune_cells_prop_sample_id_fractions))

# Combining the two data frames
combined_immune_cells_prop_sample_id <- cbind(immune_cells_prop_sample_id_ratios, immune_cells_prop_sample_id_fractions)

common_col_names <- intersect(names(immune_cells_prop_sample_id_ratios), names(immune_cells_prop_sample_id_fractions))
# Removing duplicated columns from the second data frame before binding
combined_immune_cells_prop_sample_id <- cbind(immune_cells_prop_sample_id_ratios, immune_cells_prop_sample_id_fractions[setdiff(names(immune_cells_prop_sample_id_fractions), names(immune_cells_prop_sample_id_ratios))])


# Identify the columns ending with '_counts' and '_prop'
immune_counts_cols <- grep("_counts$", names(combined_immune_cells_prop_sample_id), value = TRUE)
immune_prop_cols <- grep("_prop$", names(combined_immune_cells_prop_sample_id), value = TRUE)

# Add 'sample_id' column to both lists
immune_counts_cols_sample_id <- c("sample_id", immune_counts_cols)
immune_prop_cols_sample_id <- c("sample_id", immune_prop_cols)

# Create two subsets, one has the counts and other has the ratios. 
immune_cell_type_count <- combined_immune_cells_prop_sample_id[, immune_counts_cols_sample_id]
immune_cell_ratios <- combined_immune_cells_prop_sample_id[, immune_prop_cols_sample_id]

#remove tumor cell count columns as its redudent
columns_to_remove <- c(tumor_func_columns_core_simple_count, "Tumor_cells_counts")
columns_to_keep <- setdiff(names(immune_cell_type_count), columns_to_remove)
immune_cell_type_count <- immune_cell_type_count[, columns_to_keep]
#write_parquet(immune_cell_type_count, "../tables/immune_cell_type_count.parquet")

# generate long format feature table
immune_cell_ratios <- immune_cell_ratios %>%
  mutate(feature_type = "Cell_Abundance",
         feature_source = "whole_sample")

immune_cell_ratios <- combine_data_metadata(immune_cell_ratios, metadata, "sample_id")

setDT(immune_cell_ratios)


#immune_cell_ratios$who_grade <- as.character(immune_cell_ratios$who_grade)

melted_cell_table_immune_cell_ratios <- melt(data = immune_cell_ratios,
                                            #id.vars = c('sample_id','patient_id','ROI','site','final_diagnosis','cohort','feature_source','feature_type'),
                                            id.vars = c('sample_id','patient_id','feature_source','feature_type'),
                                            measure.vars = immune_prop_cols,
                                            variable.name = 'feature_variable',
                                            value.name = 'feature_value')



melted_cell_table_immune_cell_ratios <- melted_cell_table_immune_cell_ratios %>%
  mutate(feature_type = case_when(
    str_detect(feature_variable, "func") | str_detect(feature_variable, "all_tumor_count") ~ "immune_tumor_cell_features",
    str_detect(feature_variable, "all_immune_count") ~ "immune_cell_features",
    str_detect(feature_variable, "plus") & !str_detect(feature_variable, "func") ~ "immune_cell_features",
    TRUE ~ "immune_cell_ratios"
  ),
  feature_source = "whole_sample")



melted_cell_table_immune_cell_ratios <- melted_cell_table_immune_cell_ratios %>%
  mutate(bio_feature_type = case_when(
    str_detect(feature_variable, "over_all_immune_count_prop") ~ "Relative_to_all_immune_cells",
    str_detect(feature_variable, "all_cell_count_FOV_immune_prop") ~ "Relative_to_all_cells",
    str_detect(feature_variable, "all_tumor_count_prop") ~ "Relative_to_all_tumor_cells",
    str_detect(feature_variable, "func|plus") & str_detect(feature_variable, "B7H3|EGFR|VISTA|NG2|HER2|GPC2|GM2_GD2") ~ "Immune_to_Tumor_cells",
    str_detect(feature_variable, "func|plus") ~ "Immune_cells", # If only "func" or "plus" are found but not the tumor markers
    TRUE ~ "not_relevant"
  ))



melted_cell_table_immune_cell_ratios <- melted_cell_table_immune_cell_ratios %>%
  mutate(feature_type = case_when(
    str_detect(feature_variable, "all_tumor_count_prop ") ~ "Cell_Abundance",
    str_detect(feature_variable, "over_all_immune_count_prop") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_cell_count_FOV_immune_prop") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_tumor_count_prop") ~ "Cell_Abundance",
    str_detect(feature_variable, "func") ~ "Cell_Ratios",
    str_detect(feature_variable, "plus") ~ "Cell_Ratios",
    TRUE ~ "not_relevant"
  ))



# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_immune_cell_ratios$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_immune_cell_ratios$Broad_Feature_Type <- "Cells"

write_parquet(melted_cell_table_immune_cell_ratios, "../tables/melted_cell_table_immune_cell_ratios.parquet")

##### ratio of broad immune cells to other broad types and total immune, tumor and all cells ----


positive_count_result_cells_broad <- cell_table_all_merged_thresholded %>%
  filter(
    (cell_meta_cluster_final_broad == "Tumor" & source == "cell_table_tumor") |
      (cell_meta_cluster_final_broad != "Tumor" & source == "cell_table_immune")
  ) %>%
  count(sample_id, cell_meta_cluster_final_broad) %>%
  spread(key = cell_meta_cluster_final_broad, value = n, fill = 0)

colnames(positive_count_result_cells_broad)

colnames(positive_count_result_cells_broad)[-which(colnames(positive_count_result_cells_broad) == "sample_id")] <- 
  paste0(colnames(positive_count_result_cells_broad)[-which(colnames(positive_count_result_cells_broad) == "sample_id")], "_counts")
colnames(positive_count_result_cells_broad)

positive_count_result_immune_cells_columns <- colnames(positive_count_result_cells_broad)
positive_count_result_immune_cells_columns <- setdiff(colnames(positive_count_result_cells_broad), "sample_id")
positive_count_result_immune_cells_columns

write_parquet(positive_count_result_cells_broad, "../tables/positive_count_result_cells_broad.parquet")

# 
# # Assuming `positive_count_result_tumor` is already created
# melted_positive_count_result_immune <- melt(
#   positive_count_result_immune_cells,
#   id.vars = "sample_id",
#   measure.vars = positive_count_result_immune_cells_columns,
#   variable.name = "feature_variable",
#   value.name = "feature_value"
# )
# 
# # generate long format feature table
# melted_positive_count_result_immune <- melted_positive_count_result_immune %>%
#   mutate(feature_type = "immune_cell_features_counts",
#          feature_source = "whole_sample")
# 
# 
# # View the result
# head(melted_positive_count_result_immune)

# write_parquet(melted_positive_count_result_immune, "../tables/melted_positive_count_result_immune.parquet")
# melted_positive_count_result_immune <- read_parquet("melted_positive_count_result_immune.parquet")
# 
# 




# Merge tumor_count and non_empty_counts based on sample_id
positive_count_result_cells_broad_sample_id <- merge(positive_count_result_cells_broad, tumor_count_result, by = "sample_id")
positive_count_result_cells_broad_sample_id <- merge(positive_count_result_cells_broad_sample_id, immune_count_result, by = "sample_id")
positive_count_result_cells_broad_sample_id <- merge(positive_count_result_cells_broad_sample_id, all_cell_count_immune_FOV_result, by = "sample_id")
positive_count_result_cells_broad_sample_id <- merge(positive_count_result_cells_broad_sample_id, all_cell_count_tumor_FOV_result, by = "sample_id")

positive_count_result_cells_broad_sample_id

tumor_func_columns_core_simple
tumor_func_columns_count <- paste0(tumor_func_columns_core_simple, "_counts")
tumor_func_columns_count_sample_id <- c(tumor_func_columns_count, "sample_id")

positive_count_result_tumor_core <- positive_count_result_tumor %>% 
  dplyr::select(all_of(tumor_func_columns_count_sample_id))
positive_count_result_tumor_core

positive_count_result_cells_broad_sample_id <- merge(positive_count_result_tumor_core, positive_count_result_cells_broad_sample_id, by = "sample_id")
colnames(positive_count_result_cells_broad_sample_id)



proportion_columns <- c(unique(cell_table_immune_thresholded$cell_meta_cluster_final_broad))
proportion_columns <- setdiff(proportion_columns, "Tumor_cells")
proportion_columns <- paste0(proportion_columns, "_counts")
proportion_columns <- c(proportion_columns)
proportion_columns 


denominator_columns <- c(proportion_columns,tumor_func_columns_count,"all_tumor_count","all_immune_count","all_cell_count_FOV_immune","all_cell_count_tumor_FOV" )
denominator_columns
immune_cells_prop_sample_id <- copy(positive_count_result_cells_broad_sample_id)


denominator_columns <-  c("all_tumor_count","all_immune_count","all_cell_count_FOV_immune","all_cell_count_tumor_FOV" )
# Initialize a copy of the original data frame for simple ratios
immune_cells_prop_sample_id_ratios <- copy(positive_count_result_cells_broad_sample_id)

# Loop for simple ratios
for (num_col in proportion_columns) {
  for (denom_col in denominator_columns) {
    if (num_col == denom_col) {
      next
    }
    new_col_name = paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "prop", sep = "_")
    immune_cells_prop_sample_id_ratios[[new_col_name]] <- immune_cells_prop_sample_id[[num_col]] / immune_cells_prop_sample_id[[denom_col]]
  }
}

# Initialize a copy of the original data frame for fractions
immune_cells_prop_sample_id_fractions <- copy(positive_count_result_cells_broad_sample_id)
proportion_columns
denominator_columns_without_all <- c(tumor_func_columns_count, proportion_columns)
denominator_columns_without_all


# Assuming proportion_columns and denominator_columns_without_all are sorted or you sort them beforehand
proportion_columns <- sort(proportion_columns)
denominator_columns_without_all <- sort(denominator_columns_without_all)
denominator_columns_without_all
for (num_col in proportion_columns) {
  print(num_col)
  for (denom_col in denominator_columns_without_all) {
    print(denom_col)
    if (num_col == denom_col) {
      next
    }
    # Check if the index of num_col is less than denom_col to ensure each pair is processed once
    if (which(proportion_columns == num_col) < which(denominator_columns_without_all == denom_col)) {
      new_col_name_A = paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", num_col), "plus", sub("_counts$", "", denom_col), "prop", sep = "_")
      new_col_name_B = paste(sub("_counts$", "", denom_col), "over", sub("_counts$", "", num_col), "plus", sub("_counts$", "", denom_col), "prop", sep = "_")
      
      immune_cells_prop_sample_id_fractions[[new_col_name_A]] <- immune_cells_prop_sample_id[[num_col]] / 
        (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
      immune_cells_prop_sample_id_fractions[[new_col_name_B]] <- immune_cells_prop_sample_id[[denom_col]] / 
        (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
    }
  }
}



proportion_columns <- sort(proportion_columns)
denominator_columns_without_all <- sort(denominator_columns_without_all)

immune_cells_prop_sample_id_fractions <- list() # Initialize if not already existing



for (num_col in proportion_columns) {
  num_index <- which(proportion_columns == num_col)
  for (denom_col in denominator_columns_without_all) {
    denom_index <- which(denominator_columns_without_all == denom_col)
    
    if (num_col == denom_col || num_index >= denom_index) {
      next
    }
    
    new_col_name <- paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "plus", sub("_counts$", "", num_col), "prop", sep = "_")
    
    immune_cells_prop_sample_id_fractions[[new_col_name]] <- immune_cells_prop_sample_id[[num_col]] / 
      (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
  }
}


# Assuming immune_cells_prop_sample_id_fractions is already a defined data frame
proportion_columns <- sort(proportion_columns)
denominator_columns_without_all <- sort(denominator_columns_without_all)

for (num_col in proportion_columns) {
  num_index <- which(proportion_columns == num_col)
  for (denom_col in denominator_columns_without_all) {
    denom_index <- which(denominator_columns_without_all == denom_col)
    
    if (num_col == denom_col || num_index >= denom_index) {
      next
    }
    
    new_col_name <- paste(sub("_counts$", "", num_col), "over", sub("_counts$", "", denom_col), "plus", sub("_counts$", "", num_col), "prop", sep = "_")
    
    # Append new ratio calculations directly into the data frame
    immune_cells_prop_sample_id_fractions[[new_col_name]] <- immune_cells_prop_sample_id[[num_col]] / 
      (immune_cells_prop_sample_id[[num_col]] + immune_cells_prop_sample_id[[denom_col]])
  }
}


# Check and print out all generated column names
print(names(immune_cells_prop_sample_id_fractions))

# Combining the two data frames
combined_immune_cells_prop_sample_id <- cbind(immune_cells_prop_sample_id_ratios, immune_cells_prop_sample_id_fractions)

common_col_names <- intersect(names(immune_cells_prop_sample_id_ratios), names(immune_cells_prop_sample_id_fractions))
# Removing duplicated columns from the second data frame before binding
combined_immune_cells_prop_sample_id <- cbind(immune_cells_prop_sample_id_ratios, immune_cells_prop_sample_id_fractions[setdiff(names(immune_cells_prop_sample_id_fractions), names(immune_cells_prop_sample_id_ratios))])


# Identify the columns ending with '_counts' and '_prop'
immune_counts_cols <- grep("_counts$", names(combined_immune_cells_prop_sample_id), value = TRUE)
immune_prop_cols <- grep("_prop$", names(combined_immune_cells_prop_sample_id), value = TRUE)

# Add 'sample_id' column to both lists
immune_counts_cols_sample_id <- c("sample_id", immune_counts_cols)
immune_prop_cols_sample_id <- c("sample_id", immune_prop_cols)

# Create two subsets, one has the counts and other has the ratios. 
immune_cell_type_count <- combined_immune_cells_prop_sample_id[, immune_counts_cols_sample_id]
immune_cell_ratios <- combined_immune_cells_prop_sample_id[, immune_prop_cols_sample_id]

#remove tumor cell count columns as its redudent
columns_to_remove <- c(tumor_func_columns_core_simple_count, "Tumor_cells_counts")
columns_to_keep <- setdiff(names(immune_cell_type_count), columns_to_remove)
immune_cell_type_count <- immune_cell_type_count[, columns_to_keep]
#write_parquet(immune_cell_type_count, "../tables/immune_cell_type_count.parquet")

# generate long format feature table
immune_cell_ratios <- immune_cell_ratios %>%
  mutate(feature_type = "Cell_Abundance",
         feature_source = "whole_sample")

immune_cell_ratios <- combine_data_metadata(immune_cell_ratios, metadata, "sample_id")

setDT(immune_cell_ratios)


#immune_cell_ratios$who_grade <- as.character(immune_cell_ratios$who_grade)

melted_cell_table_immune_cell_ratios_broad <- melt(data = immune_cell_ratios,
                                             #id.vars = c('sample_id','patient_id','ROI','site','final_diagnosis','cohort','feature_source','feature_type'),
                                             id.vars = c('sample_id','patient_id','feature_source','feature_type'),
                                             measure.vars = immune_prop_cols,
                                             variable.name = 'feature_variable',
                                             value.name = 'feature_value')



melted_cell_table_immune_cell_ratios_broad <- melted_cell_table_immune_cell_ratios_broad %>%
  mutate(feature_type = case_when(
    str_detect(feature_variable, "func") | str_detect(feature_variable, "all_tumor_count") ~ "immune_tumor_cell_features",
    str_detect(feature_variable, "all_immune_count") ~ "immune_cell_features",
    str_detect(feature_variable, "plus") & !str_detect(feature_variable, "func") ~ "immune_cell_features",
    TRUE ~ "immune_cell_ratios"
  ),
  feature_source = "whole_sample")



melted_cell_table_immune_cell_ratios_broad <- melted_cell_table_immune_cell_ratios_broad %>%
  mutate(bio_feature_type = case_when(
    str_detect(feature_variable, "over_all_immune_count_prop") ~ "Relative_to_all_immune_cells",
    str_detect(feature_variable, "all_cell_count_FOV_immune_prop") ~ "Relative_to_all_cells",
    str_detect(feature_variable, "all_tumor_count_prop") ~ "Relative_to_all_tumor_cells",
    str_detect(feature_variable, "func|plus") & str_detect(feature_variable, "B7H3|EGFR|VISTA|NG2|HER2|GPC2|GM2_GD2") ~ "Immune_to_Tumor_cells",
    str_detect(feature_variable, "func|plus") ~ "Immune_cells", # If only "func" or "plus" are found but not the tumor markers
    TRUE ~ "not_relevant"
  ))



melted_cell_table_immune_cell_ratios_broad <- melted_cell_table_immune_cell_ratios_broad %>%
  mutate(feature_type = case_when(
    str_detect(feature_variable, "all_tumor_count_prop ") ~ "Cell_Abundance",
    str_detect(feature_variable, "over_all_immune_count_prop") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_cell_count_FOV_immune_prop") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_tumor_count_prop") ~ "Cell_Abundance",
    str_detect(feature_variable, "func") ~ "Cell_Ratios",
    str_detect(feature_variable, "plus") ~ "Cell_Ratios",
    TRUE ~ "not_relevant"
  ))



# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_immune_cell_ratios_broad$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_cell_table_immune_cell_ratios_broad$Broad_Feature_Type <- "Cells"

melted_cell_table_immune_cell_ratios_broad <- melted_cell_table_immune_cell_ratios_broad %>%
  filter(
    feature_variable != "Tumor_over_all_tumor_count_prop" &
      !(grepl("tumor_FOV_prop", feature_variable) & 
          feature_variable != "Tumor_over_all_cell_count_tumor_FOV_prop")
  )

melted_cell_table_immune_cell_ratios_broad

write_parquet(melted_cell_table_immune_cell_ratios_broad, "../tables/melted_cell_table_immune_cell_ratios_broad.parquet")



##### ratio of immune cells (based on functional marker expression/pixie cell phenotype) ###### 

# Data frame to store results
results <- data.frame(
  sample_id = character(),
  CellType1 = character(),
  CellType2 = character(),
  Marker1 = character(),
  Marker2 = character(),
  Proportion_A = numeric(),  # for A/(A+B)
  Proportion_B = numeric(),  # for B/(A+B)
  stringsAsFactors = FALSE
)




# generate a table where counts of each tumor marker combination is seen in each sample_id
positive_count_result_immune <- calculate_positive_counts(cell_table_immune_thresholded, complete_immune_func_columns, summary_by = "sample_id", summary_by_2 = "cell_meta_cluster_final")

# Merge tumor_count and non_empty_counts based on sample_id
cell_prop_immune_thresholded_sample_id <- merge(positive_count_result_immune, tumor_count_result, by = "sample_id")
cell_prop_immune_thresholded_sample_id <- merge(cell_prop_immune_thresholded_sample_id, immune_count_result, by = "sample_id")
cell_prop_immune_thresholded_sample_id <- merge(cell_prop_immune_thresholded_sample_id, all_cell_count_immune_FOV_result, by = "sample_id")

# Iterate over each sample_id
unique_sample_ids <- unique(cell_prop_immune_thresholded_sample_id$sample_id)


# Function to save results to a parquet file and clear memory
save_and_clear_results <- function(results, batch_number) {
  file_name <- sprintf("20240202_results_batch_min_10_total_%d.parquet", batch_number)
  arrow::write_parquet(results, file_name)
  return(data.frame(
    sample_id = character(),
    CellType1 = character(),
    CellType2 = character(),
    Marker1 = character(),
    Marker2 = character(),
    Proportion_A = numeric(),
    Proportion_B = numeric(),
    stringsAsFactors = FALSE
  ))
}

# Split unique_sample_ids into chunks of 20
sample_id_chunks <- split(unique_sample_ids, ceiling(seq_along(unique_sample_ids)/50))

batch_number <- 1
results <- data.frame(
  sample_id = character(),
  CellType1 = character(),
  CellType2 = character(),
  Marker1 = character(),
  Marker2 = character(),
  Proportion_A = numeric(),
  Proportion_B = numeric(),
  stringsAsFactors = FALSE
)

for (chunk in sample_id_chunks) {
  for (current_sample_id in chunk) {
    print(current_sample_id)
    print(batch_number)
    sample_id_data <- subset(cell_prop_immune_thresholded_sample_id, sample_id == current_sample_id)
    for (i in 1:length(names(func_immune_dictionary))) {
      cell_type_1 <- names(func_immune_dictionary)[i]
      markers_cell_type_1 <- unique(c(func_immune_dictionary[[cell_type_1]], universal_markers))

      for (j in i:length(names(func_immune_dictionary))) {
        cell_type_2 <- names(func_immune_dictionary)[j]
        markers_cell_type_2 <- unique(c(func_immune_dictionary[[cell_type_2]], universal_markers))

        for (marker_1 in markers_cell_type_1) {
          for (marker_2 in markers_cell_type_2) {
            if (i != j || marker_1 != marker_2) {
              proportions <- calculate_proportion_immune_func(sample_id_data, cell_type_1, marker_1, cell_type_2, marker_2)
              if (!is.na(proportions[[1]]) || !is.na(proportions[[2]])) {
                new_row_A <- data.frame(
                  sample_id = current_sample_id,
                  CellType1 = cell_type_1,
                  CellType2 = cell_type_2,
                  Marker1 = marker_1,
                  Marker2 = marker_2,
                  Proportion_A = proportions[[1]],
                  Proportion_B = NA,
                  stringsAsFactors = FALSE
                )
                new_row_B <- data.frame(
                  sample_id = current_sample_id,
                  CellType1 = cell_type_1,
                  CellType2 = cell_type_2,
                  Marker1 = marker_1,
                  Marker2 = marker_2,
                  Proportion_A = NA,
                  Proportion_B = proportions[[2]],
                  stringsAsFactors = FALSE
                )
                results <- rbind(results, new_row_A, new_row_B)
              }
            }
          }
        }
      }
    }
  }
  # Save the results for the current batch and clear memory
  results <- save_and_clear_results(results, batch_number)
  batch_number <- batch_number + 1
}

# List all Parquet files matching the pattern
files <- list.files(pattern = "20240202_results_batch_min_10_total_.*\\.parquet$", full.names = TRUE)

# Read and combine all files into a single dataframe
combined_results <- lapply(files, read_parquet) %>% bind_rows()


# Assuming 'results' is your existing dataframe
combined_results <- combined_results %>%
  group_by(sample_id, CellType1, CellType2, Marker1, Marker2) %>%
  summarise(Proportion_A = max(Proportion_A, na.rm = TRUE),
            Proportion_B = max(Proportion_B, na.rm = TRUE),
            .groups = 'drop')
combined_results

# Modify the 'combined_results' dataframe
immune_func_ratios <- combined_results %>%
  mutate(feature_source = 'cell_meta_cluster_final',
         feature_type = 'Cell_Ratios',
         feature_variable_A = paste(CellType1, Marker1, "over", CellType1, Marker1, "plus", CellType2, Marker2, sep = "_"),
         feature_variable_B = paste(CellType2, Marker2, "over", CellType1, Marker1, "plus", CellType2, Marker2, sep = "_"),
         feature_value_A = Proportion_A,
         feature_value_B = Proportion_B) %>%
  dplyr::select(sample_id, feature_source, feature_type, feature_variable_A, feature_value_A, feature_variable_B, feature_value_B)

# Melt the dataframe to long format
melted_immune_func_ratios <- immune_func_ratios %>%
  pivot_longer(cols = c(feature_variable_A, feature_variable_B, feature_value_A, feature_value_B),
               names_to = c(".value", "set"),
               names_pattern = "(feature_.*?)(_[AB])") %>%
  dplyr::select(-set) %>%
  inner_join(metadata, by = "sample_id") %>%
  dplyr::select(sample_id, patient_id ,feature_source, feature_type, feature_variable, feature_value)


# Create a function that generates a consistent ID for ratio comparisons, regardless of order
generate_consistent_id <- function(variable_name) {
  parts <- unlist(strsplit(variable_name, "_plus_|_over_"))
  # Sort the parts to create a unique ID that is order-independent
  id <- paste(sort(parts), collapse="_")
  return(id)
}

# Add a new column with the consistent ID for each entry
melted_immune_func_ratios <- melted_immune_func_ratios %>%
  mutate(consistent_id = sapply(feature_variable, generate_consistent_id))

# Now filter out the duplicates by keeping only the first occurrence of each unique consistent_id for each sample_id
melted_immune_func_ratios <- melted_immune_func_ratios %>%
  group_by(sample_id, consistent_id) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  dplyr::select(-consistent_id)  # Remove the helper column

# Now melted_immune_func_ratios should have the duplicates removed


melted_immune_func_ratios <- melted_immune_func_ratios %>%
  mutate(bio_feature_type = "Cells_and_functional_markers")




# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_immune_func_ratios$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_immune_func_ratios$Broad_Feature_Type <- "Cells"


# Write to Parquet file
write_parquet(melted_immune_func_ratios, "../tables/melted_immune_func_ratios.parquet")
#melted_immune_func_ratios<- read_parquet("melted_immune_func_ratios.parquet")

# generate counts df for other analysis

# Get all unique markers from the dictionary and universal markers
all_markers <- unique(c(unlist(func_immune_dictionary), universal_markers))

# Create a dataframe to store the counts for each marker per sample_id and cell_meta_cluster_final
summary_df <- unique(cell_prop_immune_thresholded_sample_id[, c("sample_id", "cell_meta_cluster_final")])

# Initialize columns for each marker with zeros
for(marker in all_markers) {
  summary_df[[marker]] <- integer(nrow(summary_df))
}

# Iterate over the rows of the summary_df to populate the marker counts
for(i in 1:nrow(summary_df)) {
  # Subset the full data for the current sample_id and cell_meta_cluster_final
  #print(current_sample_id)
  current_sample_id <- summary_df$sample_id[i]
  current_cluster <- summary_df$cell_meta_cluster_final[i]
  subset_data <- cell_prop_immune_thresholded_sample_id[cell_prop_immune_thresholded_sample_id$sample_id == current_sample_id & 
                                                    cell_prop_immune_thresholded_sample_id$cell_meta_cluster_final == current_cluster, ]
  
  # Retrieve the relevant markers for the current cell type from the dictionary
  relevant_markers <- c(func_immune_dictionary[[current_cluster]], universal_markers)
  
  # Sum the counts for each relevant marker
  for(marker in relevant_markers) {
    if(marker != "" && paste0(marker, "_func_counts") %in% names(subset_data)) {
      summary_df[i, marker] <- sum(subset_data[[paste0(marker, "_func_counts")]], na.rm = TRUE)
    }
  }
}
# Now summary_df has the summed counts for each marker across all sample_ids and cell_meta_cluster_final


immune_func_type_count <- copy(summary_df)

write_parquet(immune_func_type_count, "../tables/immune_func_type_count.parquet")

##### ratio of immune cells (based on functional marker expression/pixie cell phenotype) to all immune cells and tumor cell types or all cells #####

# generate a table where counts of each tumor marker combination is seen in each sample_id
positive_count_result_immune <- calculate_positive_counts(cell_table_immune_thresholded, complete_immune_func_columns, summary_by = "sample_id", summary_by_2 = "cell_meta_cluster_final")


positive_count_result_immune_func_columns <- colnames(positive_count_result_immune)
positive_count_result_immune_func_columns <- setdiff(
  colnames(positive_count_result_immune),
  c("sample_id", "cell_meta_cluster_final")
)
positive_count_result_immune_func_columns

# 
# # Assuming `positive_count_result_tumor` is already created
# melted_positive_count_result_immune_func <- melt(
#   positive_count_result_immune,
#   id.vars = c("sample_id","cell_meta_cluster_final"),
#   measure.vars = positive_count_result_immune_func_columns,
#   variable.name = "feature_variable",
#   value.name = "feature_value"
# )
# 
# # generate long format feature table
# melted_positive_count_result_immune_func <- melted_positive_count_result_immune_func %>%
#   mutate(feature_type = "immune_cell_features_counts",
#          feature_source = "cell_meta_cluster_final")
# 
# 
# # View the result
# head(melted_positive_count_result_immune_func)
# 
# write_parquet(melted_positive_count_result_immune_func, "melted_positive_count_result_immune_func.parquet")
# melted_positive_count_result_immune_func <- read_parquet("melted_positive_count_result_immune_func.parquet")




# Merge tumor_count and non_empty_counts based on sample_id
cell_prop_immune_thresholded_sample_id <- merge(positive_count_result_immune, tumor_count_result, by = "sample_id")
cell_prop_immune_thresholded_sample_id <- merge(cell_prop_immune_thresholded_sample_id, immune_count_result, by = "sample_id")
cell_prop_immune_thresholded_sample_id <- merge(cell_prop_immune_thresholded_sample_id, all_cell_count_immune_FOV_result, by = "sample_id")


cell_prop_immune_thresholded_sample_id
#columns_to_merge <- c(colnames(tumor_cell_type_count)[1], colnames(tumor_cell_type_count)[93:99])
columns_to_merge <- c("sample_id", colnames(tumor_cell_type_count)[sapply(gregexpr("func", colnames(tumor_cell_type_count)), function(x) sum(x > 0) == 1)])
columns_to_merge
# Subset 'tumor_cell_type_count' to only include the columns from 'columns_to_merge'
tumor_cell_type_count_subset <- tumor_cell_type_count[, columns_to_merge]
tumor_cell_type_count_subset
# Merge 'cell_prop_immune_thresholded_sample_id' with the subset of 'tumor_cell_type_count'
cell_prop_immune_thresholded_sample_id <- merge(cell_prop_immune_thresholded_sample_id, tumor_cell_type_count_subset, by = "sample_id")


cell_prop_immune_thresholded_sample_id
unique_sample_ids <- unique(cell_prop_immune_thresholded_sample_id$sample_id)


results_proportions <- data.frame(
  sample_id = character(),
  CellType = character(),
  Marker = character(),
  FunctionalGroup = character(),
  Value = numeric(),  # Renamed from Proportion to Value
  stringsAsFactors = FALSE
)

# Columns after 'all_tumor_count' or index 81
additional_functional_columns <- colnames(cell_prop_immune_thresholded_sample_id)[796:length(colnames(cell_prop_immune_thresholded_sample_id))]
additional_functional_columns <- colnames(cell_prop_immune_thresholded_sample_id)[(which(colnames(cell_prop_immune_thresholded_sample_id) == "all_tumor_count")):length(colnames(cell_prop_immune_thresholded_sample_id))] 
additional_functional_columns

# Iterate over each sample_id
for (current_sample_id in unique_sample_ids) {
  # Subset the data for the current sample_id
  sample_id_data <- subset(cell_prop_immune_thresholded_sample_id, sample_id == current_sample_id)
  print(current_sample_id)
  
  # Get the first row of the data for the current sample_id to use for total counts
  total_counts <- sample_id_data[1, additional_functional_columns]
  
  # Iterate through each cell type within the sample_id
  for (cell_type in names(func_immune_dictionary)) {
    # Combine the specific markers from the dictionary with the universal markers
    markers <- unique(c(func_immune_dictionary[[cell_type]], universal_markers))
    
    # Calculate the proportions for each marker
    for (marker in markers) {
      # Get the counts for the cell type and marker
      count <- sample_id_data[sample_id_data$cell_meta_cluster_final == cell_type, paste0(marker, "_func_counts")]
      
      # Iterate over each functional group to calculate the proportions
      for (functional_group in additional_functional_columns) {
        # Get the total count for the functional group
        total_functional <- total_counts[[functional_group]]
        
        # Determine the column name based on the type of calculation
        coname <- ifelse(functional_group %in% c("all_tumor_count", "all_immune_count","all_cell_count_FOV_immune"), 
                         paste(marker, "over", gsub("_func_counts", "_func", functional_group), sep = "_"),  # Ratio
                         paste(marker, "over", cell_type, marker, "plus", gsub("_func_counts", "_func", functional_group), sep = "_"))  # Fraction
        
        # Calculate the value as ratio or fraction
        value <- ifelse(functional_group %in% c("all_tumor_count", "all_immune_count","all_cell_count_FOV_immune"), 
                        ifelse(total_functional > 0, count / total_functional, NA),  # Ratio calculation
                        ifelse(total_functional > 0, count / (count + total_functional), NA))  # Fraction calculation
        
        # Add a row for the proportion/ratio if it's not NA and greater than 0
        if (!is.na(value) && value > 0) {
          results_proportions <- rbind(results_proportions, data.frame(
            sample_id = current_sample_id, 
            CellType = cell_type, 
            Marker = marker, 
            FunctionalGroup = coname,  # Use the constructed column name
            Value = value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}



results_proportions
# Order the results by sample_id, CellType, and Marker for better readability
results_proportions <- results_proportions[order(results_proportions$sample_id, results_proportions$CellType, results_proportions$Marker),]


# Check the results
head(results_proportions)


immune_func_to_all <- transform(results_proportions,
                                feature_source = 'cell_meta_cluster_final',
                                feature_type = 'immune_func_ratios_to_all',
                                feature_variable = paste(CellType, paste(FunctionalGroup, sep=""), sep = "_"),
                                feature_value = Value)



# Reorder the columns to place the new 'feature_source', 'feature_type', 'feature_variable', and 'feature_value' columns correctly
immune_func_to_all <- immune_func_to_all[, c("sample_id", "feature_source", "feature_type", "feature_variable", "feature_value")]


immune_func_to_all <- merge(immune_func_to_all, metadata[, c("sample_id", "patient_id")], by = "sample_id", all.x = TRUE)


setDT(immune_func_to_all)

# immune_func_ratios is already in the long format so just remove unneeded columns
melted_immune_func_to_all <- immune_func_to_all %>% 
  dplyr::select(sample_id,patient_id,feature_source, feature_type, feature_variable, feature_value)


melted_immune_func_to_all <- melted_immune_func_to_all %>%
  mutate(bio_feature_type = case_when(
    str_detect(feature_variable, "all_tumor_count") ~ "Relative_to_all_tumor_cells",
    str_detect(feature_variable, "all_immune_count") ~ "Relative_to_all_immune_cells",
    str_detect(feature_variable, "all_cell_count_FOV_immune") ~ "Relative_to_all_cells",
    str_detect(feature_variable, "plus") ~ "Cells_and_functional_markers",
    TRUE ~ "not_relevant"
  ))

melted_immune_func_to_all <- melted_immune_func_to_all %>%
  mutate(feature_type = case_when(
    str_detect(feature_variable, "all_tumor_count") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_immune_count") ~ "Cell_Abundance",
    str_detect(feature_variable, "all_cell_count_FOV_immune") ~ "Cell_Abundance",
    str_detect(feature_variable, "plus") ~ "Cell_Ratios",
    TRUE ~ "not_relevant"
  ))



# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_immune_func_to_all$Feature_Class <- "Non_Spatial"


# Add a new column called Feature_Class and assign 'Non_Spatial' as the default value
melted_immune_func_to_all$Broad_Feature_Type <- "Cells"


write_parquet(melted_immune_func_to_all, "../tables/melted_immune_func_to_all.parquet")

##### % positive functional marker per cell type ----

# generate a table where counts of each tumor marker combination is seen in each sample_id
positive_count_result_immune <- calculate_positive_counts(cell_table_immune_thresholded, complete_immune_func_columns, summary_by = "sample_id", summary_by_2 = "cell_meta_cluster_final")


positive_count_result_immune_func_columns <- colnames(positive_count_result_immune)
positive_count_result_immune_func_columns <- setdiff(
  colnames(positive_count_result_immune),
  c("sample_id", "cell_meta_cluster_final")
)
positive_count_result_immune_func_columns
# Step 1: Calculate total number of cells for each sample_id and cell_meta_cluster_final
total_cells <- cell_table_immune_thresholded %>%
  group_by(sample_id, cell_meta_cluster_final) %>%
  summarize(total_cells = n(), .groups = "drop") # Count rows (each row represents one cell)

# Step 2: Append total cell counts to positive_count_result_immune
positive_count_result_with_totals <- positive_count_result_immune %>%
  left_join(total_cells, by = c("sample_id", "cell_meta_cluster_final"))

# Step 3: Calculate % positive for each functional marker
percentage_positive <- positive_count_result_with_totals %>%
  mutate(across(all_of(positive_count_result_immune_func_columns), 
                ~ (. / total_cells) * 100, 
                .names = "percent_{col}")) %>%
  dplyr::select(sample_id, cell_meta_cluster_final, total_cells, starts_with("percent_"))

# View the result
percentage_positive

# Apply 0 to percentage values if total_cells is less than 5
percentage_positive <- percentage_positive %>%
  mutate(across(starts_with("percent_"), 
                ~ ifelse(total_cells < 5, 0, .)))

# View the modified result
print(percentage_positive)

library(tidyr)

# Melt the dataframe and add new columns
percentage_positive_melted <- percentage_positive %>%
  dplyr::select(sample_id, cell_meta_cluster_final, starts_with("percent_")) %>% # Ignore total_cells
  pivot_longer(
    cols = starts_with("percent_"), 
    names_to = "feature_variable", 
    values_to = "feature_value"
  ) %>%
  # Remove the "percent_" prefix from functional marker names
  mutate(feature_variable = gsub("^percent_", "", feature_variable)) %>%
  # Add new columns
  mutate(
    feature_source = 'cell_meta_cluster_final',
    feature_type = 'Cell_Ratios',
    Feature_Class = 'Non_Spatial',
    Broad_Feature_Type = 'Cells',
    source = "cell_table_immune",
    bio_feature_type = 'Functional_marker_positivity'
  )

# View the melted dataframe
percentage_positive_melted

percentage_positive_melted <- percentage_positive_melted %>%
  left_join(metadata %>% dplyr::select(sample_id, patient_id), by = "sample_id")


write_parquet(percentage_positive_melted, "../tables/percentage_positive_melted.parquet")




##### denisty (area) of tumor cells, immune cells (func marker), and immune cells (based on pixie cell phenotyping)  ######

# tumor_cell_type_count now has the sample_id column from tumor_metadata

melted_tumor_cell_type_area <- generate_density_feature(fov_area_df_tumor,tumor_cell_type_count, metadata, feature_type = "tumor_cell_features", feature_source ="whole_sample",  merge_column_metadata = "fov_tumor", merge_column_data = "fov")
melted_tumor_cell_type_area <- melted_tumor_cell_type_area %>%
  mutate(bio_feature_type = "Tumor_cells")
melted_tumor_cell_type_area$Feature_Class <- "Spatial"
melted_tumor_cell_type_area$Broad_Feature_Type <- "Cells"
melted_tumor_cell_type_area$feature_type <- "Area_Density"
write_parquet(melted_tumor_cell_type_area, "../tables/melted_tumor_cell_type_area.parquet")

melted_immune_func_cell_type_area <- generate_density_feature(fov_area_df_immune,immune_func_type_count, metadata, feature_type = "immune_cell_features", feature_source ="cell_meta_cluster_final", merge_column_metadata = "fov_immune", merge_column_data = "fov")
melted_immune_func_cell_type_area <- melted_immune_func_cell_type_area %>%
  mutate(bio_feature_type = "Functional_markers")
melted_immune_func_cell_type_area$Feature_Class <- "Spatial"
melted_immune_func_cell_type_area$Broad_Feature_Type <- "Cells"
melted_immune_func_cell_type_area$feature_type <- "Area_Density"
write_parquet(melted_immune_func_cell_type_area, "../tables/melted_immune_func_cell_type_area.parquet")

melted_immune_cell_type_area <- generate_density_feature(fov_area_df_immune,immune_cell_type_count, metadata, feature_type = "immune_cell_features", feature_source ="whole_sample",merge_column_metadata = "fov_immune", merge_column_data = "fov")
melted_immune_cell_type_area <- melted_immune_cell_type_area %>%
  mutate(bio_feature_type = "Immune_cells")
melted_immune_cell_type_area$Feature_Class <- "Spatial"
melted_immune_cell_type_area$Broad_Feature_Type <- "Cells"
melted_immune_cell_type_area$feature_type <- "Area_Density"
write_parquet(melted_immune_cell_type_area, "../tables/melted_immune_cell_type_area.parquet")

melted_all_tumor_area <- generate_density_feature(fov_area_df_tumor,tumor_count_result, metadata, feature_type = "tumor_cell_features", feature_source ="whole_sample",merge_column_metadata = "fov_tumor", merge_column_data = "fov" )
melted_all_tumor_area <- melted_all_tumor_area %>%
  mutate(bio_feature_type = "Tumor_cells")
melted_all_tumor_area$Feature_Class <- "Spatial"
melted_all_tumor_area$Broad_Feature_Type <- "Cells"
melted_all_tumor_area$feature_type <- "Area_Density"
write_parquet(melted_all_tumor_area, "../tables/melted_melted_all_tumor_area.parquet")

melted_all_immune_area <- generate_density_feature(fov_area_df_immune,immune_count_result, metadata, feature_type = "immune_cell_features", feature_source ="whole_sample",merge_column_metadata = "fov_immune", merge_column_data = "fov")
melted_all_immune_area <- melted_all_immune_area %>%
  mutate(bio_feature_type = "Immune_cells")
melted_all_immune_area$Feature_Class <- "Spatial"
melted_all_immune_area$Broad_Feature_Type <- "Cells"
melted_all_immune_area$feature_type <- "Area_Density"
write_parquet(melted_all_immune_area, "../tables/melted_all_immune_area.parquet")

##### spatial density ----


spatial_density_df <- read.csv('/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/analysis/spatial_analysis/combined/density/sample_table_summary.csv')

spatial_density_df <- spatial_density_df %>%
  dplyr::rename(sample_id = fov)

metadata <- read.csv('/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/metadata/metadata_complete.csv')



# Merge the dataframes based on 'sample_id' but only include 'patient_id' from the metadata
merged_df <- spatial_density_df %>%
  left_join(dplyr::select(metadata, sample_id, patient_id), by = "sample_id")

# Add the new columns 'feature_type' and 'bio_feature_type'
spatial_density_melted <- merged_df %>%
  mutate(feature_type = "immune_tumor_cell_features",
         bio_feature_type = "spatial_density")

# Print a preview of the merged dataframe

# Concatenate 'feature' and 'feature_value' columns with an underscore
spatial_density_melted <- spatial_density_melted %>%
  mutate(feature_variable = paste0(feature, "_", feature_value)) %>%
  dplyr::select(-feature_value) %>%  # Remove the old 'feature_value' column
  dplyr::select(-feature) %>%  # Remove the old 'feature_value' column
  dplyr::rename(feature_value = value)  # Rename 'value' column to 'feature_value'

head(spatial_density_melted)

spatial_density_melted$feature_variable <- ifelse(spatial_density_melted$feature_variable == "diversity_diversity", 
                                                  "spatial_density_diversity", 
                                                  spatial_density_melted$feature_variable)

spatial_density_melted$Feature_Class <- "Spatial"
spatial_density_melted$Broad_Feature_Type <- "Cells"
spatial_density_melted$feature_type <- "Spatial_Density"

# Save the melted dataframe as a Parquet file
write_parquet(spatial_density_melted, "../tables/spatial_density_melted.parquet")

##### immune tumor niches #####

# import spatial combined cell table
#cell_table_sptial_temp = fread("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/spatial_analysis/combined/spatial_DA_figures_updated/immunotherapy_50_nneighbors/merged_immunotherapy_50_nneighbors.csv")

cell_table_sptial_no_metadata = fread("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240809_combined_panel_all_with_quiche_niche.csv")

# # Rename the duplicated column in both dataframes
# colnames(cell_table_sptial_no_metadata)[which(colnames(cell_table_sptial_no_metadata) == "V1")[2]] <- "V2"
# colnames(cell_table_sptial_temp)[which(colnames(cell_table_sptial_temp) == "V1")[2]] <- "V2"
# 
# 
# cell_table_sptial_no_metadata <- merge(cell_table_sptial_no_metadata, cell_table_sptial_temp[, .(sample_id, combined_label, quiche_niche)], by = c("sample_id", "combined_label"), all.x = TRUE)


# Create a summary dataframe with counts of each unique niche per sample_id
cell_table_sptial_no_metadata_summary_df <- dcast(cell_table_sptial_no_metadata, sample_id ~ quiche_niche, value.var = "quiche_niche", fun.aggregate = length, fill = 0)
cell_table_sptial_no_metadata_summary_df <- cell_table_sptial_no_metadata_summary_df[, !grepl("^Var.2$", names(cell_table_sptial_no_metadata_summary_df))]

# Print the summary dataframe
print(cell_table_sptial_no_metadata_summary_df)


# Add the new columns
cell_table_sptial_no_metadata_summary_df$feature_type <- "niches"
cell_table_sptial_no_metadata_summary_df$feature_source <- "whole_sample"
cell_table_sptial_no_metadata_summary_df$bio_feature_type <- "spatial_cell_interactions"



# Merge in patient_id from tumor_metadata based on sample_id
cell_table_sptial_no_metadata_summary_df <- merge(cell_table_sptial_no_metadata_summary_df, metadata[, c("sample_id", "patient_id")], by = "sample_id", all.x = TRUE)



# Create the list of unique niche names
niche_all_names <- colnames(cell_table_sptial_no_metadata_summary_df)[!(colnames(cell_table_sptial_no_metadata_summary_df) %in% c('sample_id', 'patient_id', 'bio_feature_type', 'feature_source', 'feature_type'))]

# Melt the dataframe
melted_cell_table_sptial_no_metadata <- melt(cell_table_sptial_no_metadata_summary_df, 
                  id.vars = c('sample_id', 'patient_id', 'bio_feature_type', 'feature_source', 'feature_type'), 
                  measure.vars = niche_all_names, 
                  variable.name = 'feature_variable', 
                  value.name = 'feature_value')

# Print the melted dataframe
print(melted_cell_table_sptial_no_metadata)
write_parquet(melted_cell_table_sptial_no_metadata, "../tables/melted_cell_table_sptial_no_metadata.parquet")

##### neighborhood metrics #####

# import spatial combined cell table
# neighborhood_freqs_50 = fread("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/analysis/spatial_analysis/combined/neighborhood_mats/neighborhood_freqs-cell_meta_cluster_ML_radius50.csv")
# 
# 
# 
# neighborhood_freqs_50
# 
# 
# # Rename the fov column to sample_id
# neighborhood_freqs_50 <- neighborhood_freqs_50 %>%
#   rename(sample_id = fov)
# 
# # Summarize based on unique sample_id and cell_meta_cluster for all columns except label and patient_id
# summarized_neighborhood_freqs_50 <- neighborhood_freqs_50 %>%
#   select(-label) %>%
#   group_by(sample_id, cell_meta_cluster) %>%
#   summarize(across(everything(), sum, na.rm = TRUE), .groups = 'drop')
# 
# # Display the structure of the summarized data
# str(summarized_neighborhood_freqs_50)
# 
# 
# summarized_neighborhood_freqs_50 <- merge(summarized_neighborhood_freqs_50, 
#                                immune_metadata[, c("sample_id", "patient_id")], 
#                                by = "sample_id", 
#                                all.x = TRUE)
# 
# 
# # Add the new columns
# summarized_neighborhood_freqs_50$feature_type <- "neighborhood_freqs_50"
# summarized_neighborhood_freqs_50$feature_source <- "cell_meta_cluster"
# summarized_neighborhood_freqs_50$feature <- "spatial"
# 
# # Create the list of unique niche names
# summarized_neighborhood_freqs_50_names <- colnames(summarized_neighborhood_freqs_50)[!(colnames(summarized_neighborhood_freqs_50) %in% c('sample_id', 'patient_id', 'cell_meta_cluster','feature', 'feature_source', 'feature_type'))]
# summarized_neighborhood_freqs_50_names
# 
# # Melt the dataframe
# melted_summarized_neighborhood_freqs_50 <- melt(summarized_neighborhood_freqs_50, 
#                                                        id.vars = c('sample_id', 'patient_id', 'cell_meta_cluster','feature', 'feature_source', 'feature_type'), 
#                                                        measure.vars = summarized_neighborhood_freqs_50_names, 
#                                                        variable.name = 'feature_variable', 
#                                                        value.name = 'feature_value')






# Define a function to process each file
process_file <- function(file_path, metadata) {
  # Import the data
  neighborhood_freqs <- fread(cmd = paste0("cat '", file_path, "'"))
  
  # Check if the file was read correctly
  if (nrow(neighborhood_freqs) == 0) {
    stop(paste("Failed to read data from file:", file_path))
  }
  
  # Print the structure of the data for debugging
  print(str(neighborhood_freqs))
  
  # Rename the fov column to sample_id
  if (!"fov" %in% colnames(neighborhood_freqs)) {
    stop(paste("Column 'fov' not found in file:", file_path))
  }
  
  neighborhood_freqs <- neighborhood_freqs %>%
    dplyr::rename(sample_id = fov)
  
  # Detect and rename the 'cell_meta' column to 'cell_meta_cluster'
  cell_meta_col <- grep("^cell_meta", colnames(neighborhood_freqs), value = TRUE)
  if (length(cell_meta_col) != 1) {
    stop(paste("Expected exactly one 'cell_meta' column, but found:", length(cell_meta_col)))
  }
  neighborhood_freqs <- neighborhood_freqs %>%
    dplyr::rename(cell_meta_cluster = !!sym(cell_meta_col))
  
  # Summarize based on unique sample_id and cell_meta_cluster for all columns except label and patient_id
  summarized_neighborhood_freqs <- neighborhood_freqs %>%
    dplyr::select(-label) %>%
    group_by(sample_id, cell_meta_cluster) %>%
    summarize(across(everything(), mean, na.rm = TRUE), .groups = 'drop')
  
  # Merge with immune metadata
  summarized_neighborhood_freqs <- merge(summarized_neighborhood_freqs, 
                                         metadata[, c("sample_id", "patient_id")], 
                                         by = "sample_id", 
                                         all.x = TRUE)
  
  # Add the new columns
  summarized_neighborhood_freqs$feature_type <- "neighborhood_freqs"
  summarized_neighborhood_freqs$feature_source <- "cell_meta_cluster_final"
  summarized_neighborhood_freqs$bio_feature_type <- "spatial_cell_interactions"
  
  # Create the list of unique niche names
  niche_names <- colnames(summarized_neighborhood_freqs)[!(colnames(summarized_neighborhood_freqs) %in% c('sample_id', 'patient_id', 'cell_meta_cluster','bio_feature_type', 'feature_source', 'feature_type'))]
  
  # Melt the dataframe
  melted_summarized_neighborhood_freqs <- melt(summarized_neighborhood_freqs, 
                                               id.vars = c('sample_id', 'patient_id', 'cell_meta_cluster','bio_feature_type', 'feature_source', 'feature_type'), 
                                               measure.vars = niche_names, 
                                               variable.name = 'feature_variable', 
                                               value.name = 'feature_value')
  
  return(melted_summarized_neighborhood_freqs)
}

# List of file paths
file_paths <- c(
  "/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/analysis/spatial_analysis/combined/spatial_analysis_updated/neighborhood_mats/neighborhood_freqs-cell_meta_cluster_radius50.csv",
  "/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/analysis/spatial_analysis/combined/spatial_analysis_updated/neighborhood_mats/neighborhood_freqs-cell_meta_cluster_ML_radius50.csv",
  "/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/analysis/spatial_analysis/combined/spatial_analysis_updated/neighborhood_mats/neighborhood_freqs-cell_meta_cluster_IT_radius50.csv"
)

# Process each file
results <- lapply(file_paths, process_file, metadata = metadata)

# Assign each result to a separate dataframe
melted_neighborhood_freqs_50 <- results[[1]]
melted_neighborhood_freqs_ML_radius50 <- results[[2]]
melted_neighborhood_freqs_IT_radius50 <- results[[3]]

melted_neighborhood_freqs_50 <- melted_neighborhood_freqs_50 %>%
  dplyr::rename(cell_meta_cluster_final = cell_meta_cluster)
melted_neighborhood_freqs_ML_radius50 <- melted_neighborhood_freqs_ML_radius50 %>%
  dplyr::rename(cell_meta_cluster_final = cell_meta_cluster)
melted_neighborhood_freqs_IT_radius50 <- melted_neighborhood_freqs_IT_radius50 %>%
  dplyr::rename(cell_meta_cluster_final = cell_meta_cluster)

melted_neighborhood_freqs_50$Feature_Class <- "Spatial"
melted_neighborhood_freqs_50$Broad_Feature_Type <- "Cells"
melted_neighborhood_freqs_50$feature_type <- "Neighborhood_Frequencies"

melted_neighborhood_freqs_ML_radius50$Feature_Class <- "Spatial"
melted_neighborhood_freqs_ML_radius50$Broad_Feature_Type <- "Cells"
melted_neighborhood_freqs_ML_radius50$feature_type <- "Neighborhood_Frequencies"

melted_neighborhood_freqs_IT_radius50$Feature_Class <- "Spatial"
melted_neighborhood_freqs_IT_radius50$Broad_Feature_Type <- "Cells"
melted_neighborhood_freqs_IT_radius50$feature_type <- "Neighborhood_Frequencies"

# For melted_neighborhood_freqs_50
melted_neighborhood_freqs_50 <- melted_neighborhood_freqs_50 %>%
  dplyr::mutate(feature_variable = paste0(cell_meta_cluster_final, "-", feature_variable))

# For melted_neighborhood_freqs_ML_radius50
melted_neighborhood_freqs_ML_radius50 <- melted_neighborhood_freqs_ML_radius50 %>%
  dplyr::mutate(feature_variable = paste0(cell_meta_cluster_final, "-", feature_variable))

# For melted_neighborhood_freqs_IT_radius50
melted_neighborhood_freqs_IT_radius50 <- melted_neighborhood_freqs_IT_radius50 %>%
  dplyr::mutate(feature_variable = paste0(cell_meta_cluster_final, "-", feature_variable))


write_parquet(melted_neighborhood_freqs_50, "../tables/melted_neighborhood_freqs_50.parquet")
write_parquet(melted_neighborhood_freqs_ML_radius50, "../tables/melted_neighborhood_freqs_ML_radius50.parquet")
write_parquet(melted_neighborhood_freqs_IT_radius50, "../tables/melted_neighborhood_freqs_IT_radius50.parquet")

##### MALDI ----

MALDI_melted <- read_parquet("your_path/20241112_MALDI_BRUCE_melted.parquet") ## downloadable at www.bruce.parkerici.org

MALDI_melted$Feature_Class <- "Non_Spatial"

MALDI_melted$Broad_Feature_Type <- "Glycan"

MALDI_melted$feature_type <- "Relative_Intensity"

write_parquet(MALDI_melted,"your_path/20241112_MALDI_BRUCE_melted.parquet")

##### NS -----



metadata <- readRDS("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/My Drive/Projects/GBM/NanoString/Analysis/20221006_GBM_NS_pipeline_cleaned-sample_metadata.rds")
metadata_new <- read_csv("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/My Drive/Projects/GBM/NanoString/Analysis/GBM_Cohort_Meta_complete.csv")
data <- readRDS('/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/My Drive/Projects/GBM/NanoString/Analysis/20221006_GBM_NS_pipeline_cleaned.rds')

library(dplyr)

# Perform a left join to keep all observations from 'metadata', including duplicates
merged_metadata <- metadata %>%
  left_join(metadata_new, by = "Case")

# Check the number of observations to ensure it matches the expected count from 'metadata'
print(nrow(merged_metadata))

# This will ensure that all rows from 'metadata', including duplicates, are retained in the final merged dataset.


# This will merge rows based on matching values in the columns: Site_Sub, Site_Sub_2, Diagnosis, Diag_Subclass, ROI, Year, and Case
# 'all = TRUE' ensures that if there are rows in either dataframe that do not have a match in the other, they will still be included in the merged data frame (similar to a full join in SQL).


NS_processed_counts <- data@assayData$log_q
NS_processed_counts_raw <- copy(NS_processed_counts)

# Assuming metadata is a data frame with columns: sample.id, slide, Case, SegmentLabel
# And NS_processed_counts has column names that match the sample.id entries

# Create a new column in metadata that concatenates slide, Case, and SegmentLabel
metadata$combined_label <- paste(metadata$slide, metadata$Case, metadata$SegmentLabel,metadata$ROI, sep = "_")

# Create a named vector with the new labels, named by sample.id
new_colnames <- setNames(metadata$combined_label, metadata$sample.id)

# # Match the order of new_colnames to the order of columns in NS_processed_counts
# # and assign them as the new column names
# NS_processed_counts <- NS_processed_counts[, names(new_colnames)]  # To ensure the order matches
# colnames(NS_processed_counts) <- new_colnames[names(NS_processed_counts)]


# Ensure each column gets the correct name from new_colnames
colnames(NS_processed_counts) <- new_colnames[names(new_colnames) %in% colnames(NS_processed_counts)]






# Now NS_processed_counts should have the updated column names


# Assuming NS_processed_counts and metadata are already loaded in your R session

# Extract the current column names from NS_processed_counts
current_colnames <- colnames(NS_processed_counts)

# Use the current_colnames to match and order the combined_label from metadata
# Ensure that the sample.id matches the current column names
ordered_labels <- metadata$combined_label[match(current_colnames, metadata$combined_label)]

# Now assign the ordered combined labels as the new column names for NS_processed_counts
colnames(NS_processed_counts) <- ordered_labels

# NS_processed_counts should now have the updated column names

str(NS_processed_counts)


# Transpose the dataframe
NS_processed_counts_transposed <- t(NS_processed_counts)

# Convert the transposed matrix back to a dataframe if needed
NS_processed_counts_transposed <- as.data.frame(NS_processed_counts_transposed)

# Print the transposed dataframe
print(NS_processed_counts_transposed)



# Add the row names as a column
NS_processed_counts_transposed$original_rownames <- rownames(NS_processed_counts_transposed)

# Load the necessary library
library(dplyr)
library(tidyr)

# Split the row names into separate columns
NS_processed_counts_transposed <- NS_processed_counts_transposed %>%
  mutate(TMA = sub("_.*", "", original_rownames),
         patient_id = sub("^[^_]*_([^_]*)_.*", "\\1", original_rownames),
         region = sub("^([^_]*)_([^_]*)_([^_]*)_.*", "\\3", original_rownames)) %>%
  dplyr::select(original_rownames, TMA, patient_id, region, everything())

# Print the first few rows of the transformed dataframe
print(head(NS_processed_counts_transposed))

metadata <- read.csv('/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/metadata/metadata_complete.csv')


# Merge the dataframes based on the patient_id column using left join to avoid extra rows
merged_df <- left_join(NS_processed_counts_transposed, metadata, by = "patient_id")

# Print the first few rows of the merged dataframe
print(head(merged_df))

# Create a temporary dataframe to store the sum of columns 5 to 10000
temp_df <- merged_df %>%
  group_by(original_rownames) %>%
  mutate(sum_columns = rowSums(across(5:10000), na.rm = TRUE))

# Filter out duplicates where both conditions are met
merged_df_unique <- temp_df %>%
  group_by(original_rownames, sum_columns) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  dplyr::select(-sum_columns)

# Print the first few rows of the resulting dataframe
print(head(merged_df_unique))


# Melt the dataframe
NS_melted <- merged_df %>%
  pivot_longer(cols = 5:11166, 
               names_to = "feature_variable", 
               values_to = "feature_value") %>%
  mutate(feature_type = region,
         #feature_type = "Whole_FOV",
         bio_feature_type = "spatial_RNA") %>%
  dplyr::select(sample_id, patient_id, feature_type, feature_variable, bio_feature_type, feature_value)



# Print the first few rows of the melted dataframe
print(head(NS_melted))
NS_melted$Feature_Class <- "Spatial"

NS_melted$Broad_Feature_Type <- "RNA"

# Rename values in the feature_type column
NS_melted <- NS_melted %>%
  mutate(feature_type = recode(feature_type, 
                               "Immune" = "Immune_High", 
                               "Non" = "Immune_Low"))

# Save the melted dataframe as a Parquet file
write_parquet(NS_melted, "../tables/NS_melted.parquet")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### GENERATE MASTER FEATURE TABLE ######

# Now attempt to bind the rows together
master_feature_table <- bind_rows(
  melted_cell_table_marker_intensity,
  melted_cell_table_tumor_cell_ratios,
  melted_cell_table_immune_cell_ratios,
  #melted_tumor_antigen_coverage,
  melted_tumor_antigen_coverage_amounts,
  #melted_immune_func_ratios,
  melted_immune_func_to_all,
  melted_tumor_cell_type_area,
  melted_immune_func_cell_type_area,
  melted_immune_cell_type_area,
  melted_all_tumor_area,
  melted_all_immune_area,
  #melted_cell_table_sptial_no_metadata,
  melted_neighborhood_freqs_50,
  melted_neighborhood_freqs_ML_radius50,
  melted_neighborhood_freqs_IT_radius50,
  melted_tumor_antigen_segment_df,
  percentage_positive_melted,
  #melted_positive_count_result_immune,
  #melted_positive_count_result_tumor,
  MALDI_melted,
  NS_melted,
  spatial_density_melted,
  melted_cell_table_immune_cell_ratios_broad
)

master_feature_table$cell_meta_cluster_final <- gsub("Myeloid_CD11b_HLADR\\+", "Myeloid_CD11b_HLADRplus", master_feature_table$cell_meta_cluster_final)
master_feature_table$cell_meta_cluster_final <- gsub("Myeloid_CD11b_HLADR-", "Myeloid_CD11b_HLADRminus", master_feature_table$cell_meta_cluster_final)



saveRDS(master_feature_table, "../tables/20240810_master_feature_table.rds")


master_feature_table_filtered <- master_feature_table %>%
  filter(!is.na(feature_value) & !is.infinite(feature_value) & !is.nan(feature_value))

saveRDS(master_feature_table_filtered, "../tables/20240810_master_feature_table_na_removed.rds")


# Assuming your dataframe is called `metadata` and you want to combine `cohort` and `site`
metadata$Treatment <- paste(metadata$cohort, metadata$site, sep = "_")

# Replace specific values in the cohort_site column
metadata$Treatment <- ifelse(metadata$Treatment %in% c("pbta_all_CHOP", 
                                                           "unknown_CHOP", 
                                                           "openpbta_CHOP", 
                                                           "brain_cptac_2020_CHOP", 
                                                           "control_CoH", 
                                                           "unknown_Stanford", 
                                                           "control_UCLA", 
                                                           "pre_trial_UCSF", 
                                                           "non_trial_controls_UCSF", 
                                                           "0_UCSF", 
                                                           "lys_control_UCSF",
                                                           "pxa_group_UCSF"), 
                               "Treatment_Naive", metadata$Treatment)

# Replace specific values in the cohort_site column
metadata$Treatment <- ifelse(metadata$Treatment == "neoadjuvant_CoH", "Neoadjuvant_PD1_Trial_1", 
                               ifelse(metadata$Treatment %in% c("neoadjuvant_resp_UCLA", "neoadjuvant_nonresp_UCLA"), "Neoadjuvant_PD1_Trial_2", 
                                      ifelse(metadata$Treatment == "neoadjuvant_SPORE_CD27_UCSF", "Combinatorial_CD27_and_SPORE_Vaccine", 
                                             ifelse(metadata$Treatment == "neoadjuvant_SPORE_vaccine_UCSF", "SPORE_Vaccine", 
                                                    ifelse(metadata$Treatment == "neoadjuvant_lys_vaccine_UCSF", "Lysate_Vaccine", 
                                                           metadata$Treatment)))))


master_feature_table_filtered_metadata <- combine_data_metadata(master_feature_table_filtered, metadata,merge_by = 'sample_id', remove_duplicates = FALSE)


# Load dplyr package if not already loaded
# library(dplyr)

# Select only the specified columns from master_feature_table_filtered_metadata
master_feature_table_filtered_metadata <- master_feature_table_filtered_metadata %>%
  dplyr::select(sample_id, patient_id, tumor_region, final_diagnosis_simple, Treatment, immunotherapy, recurrence, paired, sex,
         idh_status, progression, who_grade, Broad_Feature_Type,bio_feature_type, cell_meta_cluster_final, Feature_Class, 
         feature_type, feature_variable, feature_value)


colnames(master_feature_table_filtered_metadata) <- gsub("final_diagnosis_simple", "Tumor_Diagnosis", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("tumor_region", "Tumor_Region", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("recurrence", "Recurrence", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("immunotherapy", "Immunotherapy", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("paired", "Longitudinal", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("sex", "Sex", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("idh_status", "IDH_R132H_Status", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("progression", "Progression", colnames(master_feature_table_filtered_metadata))
colnames(master_feature_table_filtered_metadata) <- gsub("who_grade", "WHO_grade", colnames(master_feature_table_filtered_metadata))

# Updated column names after your renaming step
cols_to_capitalize <- c("Tumor_Diagnosis", "Tumor_Region", "Recurrence", "Immunotherapy", "Longitudinal", "Sex", "IDH_R132H_Status", "Progression", "WHO_grade")

# Capitalize the first letter of each item in the specified columns
master_feature_table_filtered_metadata[cols_to_capitalize] <- lapply(master_feature_table_filtered_metadata[cols_to_capitalize], 
                                                                     function(x) { 
                                                                       if(is.character(x)) {
                                                                         sub("^(.)", "\\U\\1", x, perl=TRUE) 
                                                                       } else {
                                                                         x
                                                                       }
                                                                     })

# # Create a new column 'feature_class' based on the presence of 'spatial' in 'bio_feature_type'
# master_feature_table_filtered_metadata$feature_class <- ifelse(grepl("spatial", master_feature_table_filtered_metadata$bio_feature_type, ignore.case = TRUE), 
#                                                                "Spatial", "Non_Spatial")


master_feature_table_filtered_metadata <- master_feature_table_filtered_metadata %>%
  dplyr::select(sample_id, patient_id, Tumor_Diagnosis, WHO_grade, Immunotherapy, Treatment, Recurrence, Longitudinal, 
         Progression,Tumor_Region, IDH_R132H_Status, Sex, Feature_Class, Broad_Feature_Type, feature_type, bio_feature_type, 
         cell_meta_cluster_final, feature_variable, feature_value)

# View the reorganized dataframe
print(master_feature_table_filtered_metadata)

saveRDS(master_feature_table_filtered_metadata, "your_path/20240810_master_feature_table_na_removed_metadata.rds")

##### update metadata with only acquired FOVs across all modalities -----

# Define bio feature types
bio_feature_filter <- c("glycans", "spatial_RNA")
other_feature_filter <- unique(master_feature_table_filtered$bio_feature_type[!master_feature_table_filtered$bio_feature_type %in% bio_feature_filter])

# Create a presence table for each group
glycans_samples <- unique(master_feature_table_filtered$sample_id[master_feature_table_filtered$bio_feature_type == "glycans"])
spatial_RNA_samples <- unique(master_feature_table_filtered$sample_id[master_feature_table_filtered$bio_feature_type == "spatial_RNA"])
other_samples <- unique(master_feature_table_filtered$sample_id[master_feature_table_filtered$bio_feature_type %in% other_feature_filter])

# Create a data frame with each sample_id and its presence in each group
presence_table <- data.frame(
  sample_id = unique(c(glycans_samples, spatial_RNA_samples, other_samples)),
  glycans = as.integer(unique(c(glycans_samples, spatial_RNA_samples, other_samples)) %in% glycans_samples),
  spatial_RNA = as.integer(unique(c(glycans_samples, spatial_RNA_samples, other_samples)) %in% spatial_RNA_samples),
  other = as.integer(unique(c(glycans_samples, spatial_RNA_samples, other_samples)) %in% other_samples)
)

# Display the presence table
presence_table

colnames(presence_table) <- c("sample_id", "MALDI", "NS", "MIBI")

# Remove the MIBI, MALDI, and NS columns from metadata
metadata <- metadata[, !names(metadata) %in% c("MIBI", "MALDI", "NS")]

# Merge metadata with presence_table, keeping all rows in metadata
metadata <- merge(metadata, presence_table[, c("sample_id", "MALDI", "NS", "MIBI")], 
                  by = "sample_id", all.x = TRUE)

# Replace NA values with 0 in MIBI, MALDI, and NS columns after merging
metadata$MIBI[is.na(metadata$MIBI)] <- 0
metadata$MALDI[is.na(metadata$MALDI)] <- 0
metadata$NS[is.na(metadata$NS)] <- 0

# Display the resulting metadata dataframe
metadata

write_csv(metadata,"your_path/metadata_complete.csv")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### feature comparison list ######

## generate feature comparison df ----
# Define the column names for the new data frame
new_df_cols <- c("number", "filter_col1", "filter_col1_value", "filter_col2", "filter_col2_value", "filter_col3", "filter_col3_value", "summary_by", "group_feature_col", "group_1", "group_2")
# Assuming new_feature_table and new_df_cols have been defined earlier as shown in previous code.
# Check if 'new_feature_table' exists, if not, initialize it with column names
new_feature_table <- setNames(data.frame(matrix(ncol = length(new_df_cols), nrow = 0)), new_df_cols)

metadata <- read.csv("your_path/metadata_complete.csv")
master_feature_table_filtered <- combine_data_metadata(master_feature_table_filtered, metadata, merge_by = "sample_id", remove_duplicates = FALSE)

## adding each row ----
# Add row 1
new_row(
  row_number = 1,
  filter_col1 = "site",
  filter_col1_value = "Stanford",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "sample_id",
  group_feature_col = "tumor_region",
  group_1 = "tumor_infiltrating",
  group_2 = "tumor_core",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 2
new_row(
  row_number = 2,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "2",
  group_2 = "4",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 3
new_row(
  row_number = 3,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "2",
  group_2 = "3",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 4
new_row(
  row_number = 4,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "3",
  group_2 = "4",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 5
new_row(
  row_number = 5,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Oligodendroglioma",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "2",
  group_2 = "3",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 6
new_row(
  row_number = 6,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "who_grade",
  filter_col2_value = "4",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "final_diagnosis_simple",
  group_1 = "Astrocytoma",
  group_2 = "GBM",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)
# Add row 7
new_row(
  row_number = 7,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Astrocytoma_-_Oligodendroglioma_-_GBM",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "2_-_3",
  group_2 = "4",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 8
new_row(
  row_number = 8,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Astrocytoma",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "2",
  group_2 = "3",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 9
new_row(
  row_number = 9,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Astrocytoma",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "2",
  group_2 = "4",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 10
new_row(
  row_number = 10,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Astrocytoma",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "who_grade",
  group_1 = "3",
  group_2 = "4",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 11
new_row(
  row_number = 11,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Oligodendroglioma",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "recurrence",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 12
new_row(
  row_number = 12,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Astrocytoma",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "recurrence",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 13
new_row(
  row_number = 13,
  filter_col1 = "progression",
  filter_col1_value = "unknown",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "recurrence",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 14
new_row(
  row_number = 14,
  filter_col1 = "site",
  filter_col1_value = "UCLA",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "control",
  group_2 = "neoadjuvant_resp_-_neoadjuvant_nonresp",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 15
new_row(
  row_number = 15,
  filter_col1 = "site",
  filter_col1_value = "UCLA",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "control",
  group_2 = "neoadjuvant_resp",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 16
new_row(
  row_number = 16,
  filter_col1 = "site",
  filter_col1_value = "UCLA",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "control",
  group_2 = "neoadjuvant_nonresp",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 17
new_row(
  row_number = 17,
  filter_col1 = "site",
  filter_col1_value = "UCLA",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "neoadjuvant_resp",
  group_2 = "neoadjuvant_nonresp",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 18
new_row(
  row_number = 18,
  filter_col1 = "site",
  filter_col1_value = "CoH",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "control",
  group_2 = "neoadjuvant",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 19
new_row(
  row_number = 19,
  filter_col1 = "site",
  filter_col1_value = "UCSF",
  filter_col2 = "recurrence_status_paired",
  filter_col2_value = "recurrence_1",
  filter_col3 = "group",
  filter_col3_value = "A_-_C",
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 20
new_row(
  row_number = 20,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Oligodendroglioma",
  filter_col3 = "group",
  filter_col3_value = "A_-_C",
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 21
new_row(
  row_number = 21,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "Astrocytoma",
  filter_col3 = "group",
  filter_col3_value = "A_-_C",
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 22
new_row(
  row_number = 22,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "group",
  filter_col3_value = "B_-_C",
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "pre_trial_-_non_trial_controls",
  group_2 = "neoadjuvant_SPORE_vaccine",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 23
new_row(
  row_number = 23,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "group",
  filter_col3_value = "B_-_C",
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "pre_trial_-_non_trial_controls",
  group_2 = "neoadjuvant_SPORE_CD27",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 24
new_row(
  row_number = 24,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "group",
  filter_col3_value = "B",
  summary_by = "patient_id",
  group_feature_col = "cohort",
  group_1 = "neoadjuvant_SPORE_vaccine",
  group_2 = "neoadjuvant_SPORE_CD27",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 25
new_row(
  row_number = 25,
  filter_col1 = "site",
  filter_col1_value = "UCSF",
  filter_col2 = "immunotherapy",
  filter_col2_value = "no",
  filter_col3 = "group",
  filter_col3_value = "A_-_B_-_C_-_D",
  summary_by = "sample_id",
  group_feature_col = "recurrence",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 26
new_row(
  row_number = 26,
  filter_col1 = "site",
  filter_col1_value = "UCSF",
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "group",
  filter_col3_value = "A_-_B_-_C",
  summary_by = "sample_id",
  group_feature_col = "progression",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 27
new_row(
  row_number = 27,
  filter_col1 = "site",
  filter_col1_value = "UCSF",
  filter_col2 = "immunotherapy",
  filter_col2_value = "no",
  filter_col3 = "paired",
  filter_col3_value = "yes",
  summary_by = "sample_id",
  group_feature_col = "recurrence",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 28
new_row(
  row_number = 28,
  filter_col1 = "final_diagnosis_simple",
  filter_col1_value = "GBM",
  filter_col2 = "recurrence",
  filter_col2_value = "yes",
  filter_col3 = NA,
  filter_col3_value = NA,
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 29
new_row(
  row_number = 29,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "comp_1",
  group_1 = "core",
  group_2 = "recurring",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)


new_row(
  row_number = 30,
  filter_col1 = NA,
  filter_col1_value = NA,
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "comp_1",
  group_1 = "edge",
  group_2 = "recurring",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)


new_row(
  row_number = 31,
  filter_col1 = "recurrence",
  filter_col1_value = "no",
  filter_col2 = "final_diagnosis_simple",
  filter_col2_value = "GBM",
  filter_col3 = "immunotherapy",
  filter_col3_value = "no",
  summary_by = "patient_id",
  group_feature_col = "survival_status",
  group_1 = "short",
  group_2 = "long",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 32
new_row(
  row_number = 32,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "group",
  filter_col3_value = "A_-_B_-_C",
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

# Add row 33
new_row(
  row_number = 33,
  filter_col1 = "recurrence_status_paired",
  filter_col1_value = "recurrence_1",
  filter_col2 = NA,
  filter_col2_value = NA,
  filter_col3 = "group",
  filter_col3_value = "B_-_C",
  summary_by = "patient_id",
  group_feature_col = "immunotherapy",
  group_1 = "no",
  group_2 = "yes",
  data_table = master_feature_table_filtered,
  new_feature_table = new_feature_table,
  new_df_cols = new_df_cols
)

new_feature_table


## save/read feature comparison table  ----

write.csv(new_feature_table,"../tables/comparison_criteria.csv")
feature_comp_table <- read_csv("../tables/comparison_criteria.csv")


## get unique values in column in master feature table ----
selected_col_names <- colnames(master_feature_table_filtered)[1:15]

# Create a list to store unique values for each column
unique_values_list <- vector("list", length(selected_col_names))

# Populate the list with unique values for each column
for (i in seq_along(selected_col_names)) {
  col_name <- selected_col_names[i]
  unique_values_list[[i]] <- unique(master_feature_table_filtered[[col_name]])
}

# Find the maximum length to pad shorter columns
max_length <- max(sapply(unique_values_list, length))

# Pad shorter columns with NA
padded_unique_values <- lapply(unique_values_list, function(x) {
  c(x, rep(NA, max_length - length(x)))
})

# Create a dataframe from the padded list
unique_values_df <- as.data.frame(padded_unique_values, col.names = selected_col_names)



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### GLM analysis ######

## Date: Feb 1st, 2024
## Author: Meelad Amouzgar, amouzgar@stanford.edu
## runs a simple regression to identify univariate feature differences

library(tidyverse)
options(dplyr.summarise.inform = FALSE)

combine_data_metadata <- function(data, metadata, merge_by) {
  # Remove duplicates based on the merge_by column in metadata
  metadata <- metadata %>%
    distinct(!!rlang::sym(merge_by), .keep_all = TRUE)
  
  # Find duplicate columns in metadata that are also in data except for merge_by
  duplicate_columns <- setdiff(intersect(colnames(metadata), colnames(data)), merge_by)
  
  # Remove the duplicate columns from metadata
  if (length(duplicate_columns) > 0) {
    metadata <- metadata %>%
      dplyr::select(-one_of(duplicate_columns))
  }
  
  # Merge data and metadata using the column specified by merge_by
  data_merged <- merge(data, metadata, by = merge_by, all.x = TRUE)
  
  return(data_merged)
}

ft = readRDS('../tables/20240810_master_feature_table_na_removed.rds')

metadata <- read.csv("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/metadata/metadata_complete.csv")

ft <- combine_data_metadata(ft, metadata, merge_by = "sample_id")
# Exclude rows with NA in the feature_value column using na.omit
ft <- ft %>% filter(!is.na(feature_value))


cc = data.table::fread('../tables/comparison_criteria.csv',sep=',', stringsAsFactors = F)
cc
# cc <- cc %>%
#    filter(number %in% c(1))
# 
# cc <- cc %>%
#   mutate(number = ifelse(number == 1, 1, number))
# 
# cc <- cc %>%
#   mutate(V1 = ifelse(V1 == 1, 1, V1))
# 
# cc
FILTER_DATA_FXN = function(cc,...){
  DDD=apply(cc, 1, function(cc_ro){
    message(cc_ro[['number']])
    ddd=ft
    if (!is.na(cc_ro[['filter_col1']])) {
      filt_values = strsplit(cc_ro[['filter_col1_value']], split = '_-_') %>% unlist()
      ddd = ddd[ ddd[[cc_ro[['filter_col1']]]] %in%filt_values, ]
    }
    if (!is.na(cc_ro[['filter_col2']])) {
      filt_values = strsplit(cc_ro[['filter_col2_value']], split = '_-_') %>% unlist()
      ddd = ddd[ ddd[[cc_ro[['filter_col2']]]] %in% filt_values, ]
    }
    if (!is.na(cc_ro[['filter_col3']])) {
      filt_values = strsplit(cc_ro[['filter_col3_value']], split = '_-_') %>% unlist()
      ddd = ddd[ ddd[[cc_ro[['filter_col3']]]] %in% filt_values, ]
    }
    return(ddd)
  })
  names(DDD) = cc[['number']]
  return(DDD)
}

MY_GLM_FXN = function(HHH,FEATURE_NAME, cc_ro_data, ...){
  crit_num=cc_ro_data[['number']]
  group1=cc_ro_data[['group_1']]
  group2=cc_ro_data[['group_2']]
  group1_vec = strsplit(group1, split = '_-_') %>% unlist()
  group2_vec = strsplit(group2, split = '_-_') %>% unlist()
  groups=paste(group1, group2, sep = '_VS_')
  FEATURE_VALUES ='feature_value'
  GROUP_VAR= cc_ro_data[['group_feature_col']]
  SUMMARIZED_BY = cc_ro_data[['summary_by']]
  
  ## if multiple groups are being combined into 1 for a pairwise comparison, concatenate them now
  HHH[[GROUP_VAR]] = ifelse(HHH[[GROUP_VAR]]  %in% group1_vec, group1, HHH[[GROUP_VAR]] )
  HHH[[GROUP_VAR]] = ifelse(HHH[[GROUP_VAR]]  %in% group2_vec, group2, HHH[[GROUP_VAR]] )
  HHH = HHH[which(HHH[[GROUP_VAR]] %in% c(group1, group2)), ]
  
  if(cc_ro_data[['summary_by']]=='patient_id'){
    HHH=HHH %>% group_by(patient_id, !!sym(GROUP_VAR) ) %>% summarize(feature_value=mean(feature_value))
  }
  
  
  HHH[['GROUP']] = HHH[[GROUP_VAR]] 
  myFormula =  paste0(FEATURE_VALUES, '~GROUP') %>% as.formula()
  
  dif_res = data.frame(criteria_num = crit_num,
                       feature = FEATURE_NAME, 
                       group_comparison = groups,
                       group_feature_col = GROUP_VAR ,
                       model_condition = 'GLM',
                       summary_by=SUMMARIZED_BY,
                       check.rows = F,check.names = F, stringsAsFactors =F )
  
  if (nrow(HHH)<3) {
    dif_res$model_condition= 'sample_set<3'
    return(dif_res)
  }  
  if (length(na.omit(unique(HHH[['GROUP']]))) <2){
    dif_res$model_condition = 'only0or1_group_variables_in_sample_set'
    return(dif_res)
  }
  table_group_counts = table(HHH[['GROUP']])
  if (table_group_counts[[1]]<3){
    dif_res$model_condition = 'group1_LessThan3Samples'
    return(dif_res)
  }
  if (table_group_counts[[2]]<3){
    dif_res$model_condition = 'group2_LessThan3Samples'
    return(dif_res)
  }
  if (is.na(sd(HHH[['feature_value']]))){
    dif_res$model_condition = 'NA_st.dev'
    return(dif_res)
  }
  if (sd(HHH[['feature_value']])==0){
    dif_res$model_condition = '0_st.dev'
    return(dif_res)
  }
  
  ## scale feture value
  HHH[['feature_value']]  = scale(HHH[['feature_value']] , center = T, scale =T)
  
  glm.res = glm(formula = myFormula , family=gaussian(), data = HHH )
  glm.summary = summary(glm.res)
  glm.coefs = data.frame(glm.summary$coefficients, check.rows = F, check.names = F)
  
  dif_res = glm.coefs %>%
    mutate(criteria_num = crit_num,
           feature = FEATURE_NAME, 
           group_comparison = groups,
           group_feature_col = GROUP_VAR ,
           group_comparison_check = paste0(unique(HHH[['GROUP']], collapse = '_-_')),
           coef_name = rownames(glm.coefs),
           model_condition = 'GLM',
           summary_by=SUMMARIZED_BY)
  return(dif_res)
}

MY_GLM_FXN_LOOP_WRAPPER = function(GGG, cc_ro_data, feature_list){
  DIF_RESULTS = lapply(feature_list, function(FEATURE_NAME){
    HHH=GGG[GGG[['SPLIT_GROUP']] == FEATURE_NAME,]
    dif_results = MY_GLM_FXN(HHH, FEATURE_NAME, cc_ro_data)
    return(dif_results)
  }) %>% bind_rows()
  return(DIF_RESULTS)
}

LOOP_THROUGH_FEATURES = function(DDD_GROUP_SPLIT, MY_COMPARISON_LIST,...){
  lapply(MY_COMPARISON_LIST, function(cc_ro_data){ ## loop through comparisons
    print(cc_ro_data)
    message(cc_ro_data[['number']])
    FFF = DDD_GROUP_SPLIT[[cc_ro_data[['number']]]]
    
    DIF_RESULTS = lapply(names(FFF), function(FFF_name){ ## loop through celltypes (already split)
      GGG = FFF[[FFF_name]]
      feature_list = unique(GGG$SPLIT_GROUP)
      DIF_RESULTS = MY_GLM_FXN_LOOP_WRAPPER(GGG,cc_ro_data, feature_list) %>%
        mutate(SPLIT_GROUP =FFF_name)
    }) %>% bind_rows()
  }) %>% bind_rows() 
}
ft$cell_meta_cluster_final = ifelse(is.na(ft$cell_meta_cluster_final ), 'NOT_cell_meta_cluster_final',ft$cell_meta_cluster_final)
ft$SPLIT_GROUP = paste(ft$cell_meta_cluster_final,ft$feature_type, ft$feature_variable, sep = '_-_') ## PUT ALL SPLIT GROUPS HERE

start= Sys.time()
DDD = FILTER_DATA_FXN(cc) ## extracts relevant filtered data for the comparison criteria and splits into a list
DDD_GROUP_SPLIT = lapply(DDD, function(EEE){ split(EEE ,EEE$SPLIT_GROUP)})
MY_COMPARISON_LIST =split(cc, cc$number)

myDATA = LOOP_THROUGH_FEATURES(DDD_GROUP_SPLIT,MY_COMPARISON_LIST) #%>% bind_rows()
myDATA2 = myDATA %>% 
  dplyr::filter(!grepl('Intercept',coef_name) ) %>%
  mutate(coef_name = gsub('GROUP','',coef_name)) %>%
  dplyr::rename(pval = `Pr(>|t|)`,
                coef = Estimate) %>%
  dplyr::select(-c(`Std. Error`, `t value`)) %>%
  mutate(padj = p.adjust(pval, method = 'BH'))
end = Sys.time()
print(end-start)

write.table(myDATA2,'../tables/glm_feature_table_results_filtered.csv',sep =',', row.names = F, col.names = T )




##### visualize DE from GLM ----

feature_comp_table <- read_csv("../tables/comparison_criteria.csv")

glm_feature_table_results <- read_csv("../tables/glm_feature_table_results_filtered.csv")

feature_comp_table$criteria <- paste(feature_comp_table$filter_col1_value, feature_comp_table$filter_col2_value,feature_comp_table$filter_col3 ,feature_comp_table$filter_col3_value , sep = "_-_")
temp_df <- feature_comp_table[, c("number", "criteria")]
glm_feature_table_results <- merge(glm_feature_table_results, temp_df, by.x = "criteria_num", by.y = "number", all.x = TRUE)
glm_feature_table_results <- glm_feature_table_results[c("criteria", setdiff(names(glm_feature_table_results), "criteria"))]


glm_feature_table_results= glm_feature_table_results %>% mutate(SPLIT_GROUP = gsub('NA_-_','NOT_cell_meta_cluster_final_-_',SPLIT_GROUP))


glm_feature_table_results$neg_log_pval <- -log10(glm_feature_table_results$pval)

# Generate volcano plot
ggplot(glm_feature_table_results, aes(x = coef, y = neg_log_pval)) +
  geom_point(aes(color = pval < 0.05 & abs(coef) > 1)) + # Threshold for significance can be adjusted
  theme_minimal() +
  labs(x = "Coefficient (coef)",
       y = "-log10(p-value)",
       title = "Volcano Plot") +
  scale_color_manual(values = c("black", "red")) # Non-significant in blue, significant in red


plots <- glm_feature_table_results %>%
  group_by(group_comparison, group_feature_col, criteria_num) %>%
  group_split() %>%
  lapply(function(data) {
    ggplot(data, aes(x = coef, y = neg_log_pval, color = ifelse((coef > 1 | coef < -1) &pval < 0.05, "Highlighted", "Normal"))) +
      geom_point() +
      theme_minimal() +
      theme(plot.title = element_text(size = 10, lineheight = 0.8),
            legend.position = "none") +  # Remove legends
      labs(title = str_wrap(paste(data$group_comparison[1], data$group_feature_col[1], data$criteria[1], sep = " - "), width = 30),
           x = "Coefficient",
           y = "-log(p-value)") +
      scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black"))
  })


# Combine all plots
combined_plot <- wrap_plots(plots)

# Print or save the combined plot
print(combined_plot)
# ggsave("combined_volcano_plots.png", combined_plot, width = 20, height = 15)


glm_feature_table_results_sig <- glm_feature_table_results %>%
  filter(pval < 0.05, (coef > 1 | coef < -1))
#

glm_feature_table_results_sig <- glm_feature_table_results_sig %>%
  separate(SPLIT_GROUP, into = c("one", "two", "three"), sep = "_-_", remove = FALSE)


# Prepare data for plotting
glm_feature_table_results_sig$neg_log_pval <- -log10(glm_feature_table_results_sig$pval)


# get sig feature types for each comparison and plot it again total features to get % 

glm_feature_table_results <- glm_feature_table_results %>%
  separate(SPLIT_GROUP, into = c("one", "two", "three"), sep = "_-_", remove = FALSE)


# Step 1: Count occurrences in the 'two' column for the full dataset
total_counts <- glm_feature_table_results %>%
  count(two, criteria, group_comparison, group_feature_col, name = "total_n")

# Step 2: Count occurrences in the 'two' column for the significant subset
sig_counts <- glm_feature_table_results_sig %>%
  count(two, criteria, group_comparison, group_feature_col, name = "sig_n")

# Step 3: Merge the total counts with the significant counts to calculate percentages
percentage_counts <- sig_counts %>%
  left_join(total_counts, by = c("two", "criteria", "group_comparison", "group_feature_col")) %>%
  mutate(percentage = (sig_n / total_n) * 100) %>%
  arrange(criteria, group_comparison, group_feature_col, desc(percentage))

# Step 4: Create the bar plots tiled together using facet_wrap, with percentages
p <- ggplot(percentage_counts, aes(x = two, y = percentage, fill = two)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(percentage, 1)), vjust = 0, size = 3.5) +  # Adjust text position, rounding to 1 decimal
  facet_wrap(~ criteria + group_comparison + group_feature_col, scales = "free_y", ncol = 5) +  # Free y-scale for each facet
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks if desired
    legend.position = "bottom"      # Adjust legend position if needed
  ) +
  labs(x = "", y = "Percentage (%)")  # Adjust y-axis label to indicate percentage

print(p)












#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# cell - glycan correlation ----

# Base list of immune cell types
immune_cell_types_broad <- c(
  "Immune_over_all_cell_count_FOV_immune_prop", 
  "Endothelial_over_all_cell_count_FOV_immune_prop", 
  "Neurons_over_all_cell_count_FOV_immune_prop", 
  "Tumor_over_all_cell_count_tumor_FOV_prop"
)

glycan_cols <- c(
  "fucosylated",
  "sialylated",
  "highMannose",
  "hybrid",
  "paucimannose",
  "agalactosylated",
  "biantennary",
  "triantennary",
  "tetraantennary",
  "polylacnac"
)

master_feature_table_filtered <- combine_data_metadata(master_feature_table_filtered, metadata, merge_by = "sample_id", remove_duplicates = FALSE)


filtered_master_feature_table <- master_feature_table_filtered %>%
  filter(
    grepl("Stanford", sample_id) &
      (feature_variable %in% glycan_cols |feature_variable %in% immune_cell_types_broad )
  )

filtered_master_feature_table <- filtered_master_feature_table %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"))

# Pivot the filtered master feature table to have features as columns
master_feature_pivot <- filtered_master_feature_table %>%
  dplyr::select(patient_id, feature_variable, feature_value) %>%
  pivot_wider(names_from = feature_variable, values_from = feature_value, values_fn = mean) %>%
  dplyr::select(patient_id, everything())


# Identify column names from master_feature_pivot
valid_columns <- colnames(master_feature_pivot)

# Filter the combined list to keep only valid column names
filtered_combined <- immune_cell_types_broad[immune_cell_types_broad %in% valid_columns]

# View the filtered combined list
filtered_combined

# Find the elements in combined that are not in master_feature_pivot
removed_elements <- immune_cell_types_broad[!immune_cell_types_broad %in% valid_columns]

# Display removed elements
removed_elements

# Start the timer
tic("Correlation computation")

# Initialize progress bar
num_steps <- length(filtered_combined)  # Track progress per RNA column
pb <- progress_bar$new(
  format = "Computing [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = num_steps, clear = FALSE, width = 60
)

# Initialize matrices to store correlations and p-values
correlation_matrix <- matrix(NA, nrow = length(filtered_combined), ncol = length(glycan_cols))
pvalue_matrix <- matrix(NA, nrow = length(filtered_combined), ncol = length(glycan_cols))
rownames(correlation_matrix) <- filtered_combined
colnames(correlation_matrix) <- glycan_cols
rownames(pvalue_matrix) <- filtered_combined
colnames(pvalue_matrix) <- glycan_cols

# Perform correlation with progress tracking
for (i in seq_along(filtered_combined)) {
  for (j in seq_along(glycan_cols)) {
    # Perform correlation test to get both correlation and p-value
    test_result <- cor.test(
      master_feature_pivot[[filtered_combined[i]]],
      master_feature_pivot[[glycan_cols[j]]],
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    correlation_matrix[i, j] <- test_result$estimate  # Store correlation coefficient
    pvalue_matrix[i, j] <- test_result$p.value  # Store p-value
  }
  pb$tick()  # Update progress bar
}

# Stop the timer
toc()

# View the correlation and p-value matrices
print(correlation_matrix)
print(pvalue_matrix)


# Load required libraries
library(ggplot2)
library(reshape2)
library(stats)  # For hclust and dist

# Perform hierarchical clustering on rows (features) and columns (glycans)
row_order <- hclust(dist(correlation_matrix))$order
col_order <- hclust(dist(t(correlation_matrix)))$order

# Reorder the correlation and p-value matrices based on hierarchical clustering
correlation_matrix_clustered <- correlation_matrix[row_order, col_order]
pvalue_matrix_clustered <- pvalue_matrix[row_order, col_order]

# Convert the clustered matrices to long format for ggplot2
cor_long <- melt(correlation_matrix_clustered, varnames = c("Feature", "Glycan"), value.name = "Correlation")
pval_long <- melt(pvalue_matrix_clustered, varnames = c("Feature", "Glycan"), value.name = "P_value")

# Merge the two long-format data frames
plot_data <- merge(cor_long, pval_long, by = c("Feature", "Glycan"))

# Create a new column to mark non-significant correlations (p > 0.05)
plot_data$Significant <- ifelse(plot_data$P_value > 0.1, "grey", NA)

# Plot the bubble plot with clustering applied
ggplot(plot_data, aes(x = Glycan, y = Feature, size = -log10(P_value))) +
  # Plot non-significant values in grey
  geom_point(data = subset(plot_data, Significant == "grey"), color = "grey", alpha = 0.7) +
  # Plot significant correlations with a color gradient
  geom_point(data = subset(plot_data, is.na(Significant)), 
             aes(color = Correlation), alpha = 0.7) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey"
  ) +
  scale_size_continuous(range = c(2, 10), breaks = c(1, 2, 3, 5)) +  # Control size range for bubbles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Bubble Plot of Correlation Coefficients and P-values (Clustered)",
    x = "Glycan Features",
    y = "Predictor Features",
    color = "Correlation",
    size = "-log10(P-value)"
  )


# Filter for significant correlations (p-value < 0.5)
significant_data <- subset(plot_data, P_value < 0.1)


# Filter master_feature_pivot based on significant features and glycans
selected_features <- unique(significant_data$Feature)
selected_glycans <- unique(significant_data$Glycan)


# Ensure selected features and glycans are present in the column names of master_feature_pivot
valid_features <- intersect(selected_features, colnames(master_feature_pivot))
valid_features
valid_glycans <- intersect(selected_glycans, colnames(master_feature_pivot))
valid_glycans

# Filter master_feature_pivot based on valid features and glycans
filtered_data <- master_feature_pivot[, c("patient_id", valid_features, valid_glycans)]
filtered_data

# Ensure that significant_data contains valid features and glycans
print((significant_data))

# Ensure that feature and glycan are treated as character strings
significant_data$Feature <- as.character(significant_data$Feature)
significant_data$Glycan <- as.character(significant_data$Glycan)

plots <- list()  # Initialize an empty list to store plots

# Iterate through each row in significant_data to create scatter plots
for (i in 1:nrow(significant_data)) {
  feature <- significant_data$Feature[i]
  glycan <- significant_data$Glycan[i]
  
  # Debugging: Print each feature and glycan combination
  print(paste("Plotting:", feature, "vs", glycan))
  
  # Create scatter plot using ggplot
  p <- ggplot(filtered_data, aes(x = .data[[feature]], y = .data[[glycan]])) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = paste(feature, "vs", glycan),
         x = feature, 
         y = glycan) +
    theme(legend.position = "none")
  
  # Store the plot in the list
  plots[[i]] <- p
}

# Combine all plots using patchwork
if (length(plots) > 0) {
  combined_plot <- wrap_plots(plots, ncol = 4)  # Adjust layout if needed
  print(combined_plot)
} else {
  print("No valid plots were generated.")
}

saveRDS(significant_data, file = "your_path/glycan_cell_significant_data_fig_5.rds")

significant_data <- readRDS("your_path/glycan_cell_significant_data_fig_5.rds")


# RNA - glycan correlations----

RNA_cols <- unique(master_feature_table_filtered$feature_variable[master_feature_table_filtered$bio_feature_type == "spatial_RNA"])

RNA_cols


glycan_cols <- c(
  "fucosylated",
  "sialylated",
  "highMannose",
  "hybrid",
  "paucimannose",
  "agalactosylated",
  "biantennary",
  "triantennary",
  "tetraantennary",
  "polylacnac"
)

filtered_master_feature_table <- master_feature_table_filtered %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"))

filtered_master_feature_table <- filtered_master_feature_table %>%
  filter(
    feature_variable %in% glycan_cols |
      (feature_variable %in% RNA_cols & bio_feature_type == "spatial_RNA")
  )
# Pivot the filtered master feature table to have features as columns
master_feature_pivot <- filtered_master_feature_table %>%
  dplyr::select(patient_id, feature_variable, feature_value) %>%
  pivot_wider(names_from = feature_variable, values_from = feature_value, values_fn = mean) %>%
  dplyr::select(patient_id, everything())


# Start the timer
tic("Correlation computation")

# Initialize progress bar
num_steps <- length(RNA_cols) # Track progress per RNA column
pb <- progress_bar$new(
  format = "Computing [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = num_steps, clear = FALSE, width = 60
)

# Initialize an empty matrix to store correlations
correlation_matrix <- matrix(NA, nrow = length(RNA_cols), ncol = length(glycan_cols))
rownames(correlation_matrix) <- RNA_cols
colnames(correlation_matrix) <- glycan_cols

# Perform correlation with progress tracking
for (i in seq_along(RNA_cols)) {
  for (j in seq_along(glycan_cols)) {
    correlation_matrix[i, j] <- cor(
      master_feature_pivot[[RNA_cols[i]]],
      master_feature_pivot[[glycan_cols[j]]],
      use = "pairwise.complete.obs",
      method = "pearson"
    )
  }
  pb$tick()  # Update progress bar
}

# Stop the timer
toc()

# View the correlation matrix
print(correlation_matrix)

saveRDS(correlation_matrix, file = "your_path/correlation_matrix_RNA_glycan.rds")




# Convert the correlation matrix to a long format
correlation_long <- melt(correlation_matrix, varnames = c("RNA", "Glycan"), value.name = "Correlation")

# Extract gene lists with correlations > 0.5 for each glycan
glycan_gene_lists <- correlation_long %>%
  filter(Correlation > 0.4) %>%
  group_by(Glycan) %>%
  summarise(Gene_List = list(RNA))

# Print the gene lists for inspection
print(glycan_gene_lists)


saveRDS(glycan_gene_lists, file = "your_path/glycan_gene_lists_fig_5.rds")

glycan_gene_lists <- readRDS("your_path/glycan_gene_lists_fig_5.rds")

# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)

# Initialize progress bar
go_pb <- progress_bar$new(
  format = "GO Analysis [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = nrow(glycan_gene_lists), clear = FALSE, width = 60
)

# GO analysis for each glycan's gene list with progress tracking
go_results <- list()  # Initialize empty list to store results

for (i in seq_len(nrow(glycan_gene_lists))) {
  genes <- unlist(glycan_gene_lists$Gene_List[[i]])
  glycan <- glycan_gene_lists$Glycan[i]
  
  # Perform GO analysis
  result <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",  # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # Add glycan name as an attribute to the result
  if (!is.null(result)) {
    attr(result, "glycan_name") <- glycan
  }
  
  # Store the result in the list with the glycan name as key
  go_results[[glycan]] <- result
  
  # Update progress bar
  go_pb$tick()
}

# Iterate over go_results and print glycan names
cat("Glycan Names in GO Results:\n")

for (i in seq_along(go_results)) {
  result <- go_results[[i]]
  
  if (!is.null(result)) {
    glycan_name <- attr(result, "glycan_name")  # Extract glycan name
    cat("Index", i, ":", glycan_name, "\n")
  } else {
    cat("Index", i, ": No glycan name (NULL result)\n")
  }
}


saveRDS(go_results, file = "your_path/go_results_fig_5.rds")

go_results <- readRDS("your_path/go_results_fig_5.rds")


# Load the progress package
library(progress)

# Initialize lists to store descriptions and gene lists from significant pathways
all_descriptions <- c()
all_gene_lists <- list()

# Collect descriptions and gene lists from significant pathways in all GO results
for (i in seq_along(go_results)) {
  result <- go_results[[i]]
  
  if (!is.null(result)) {
    # Try to get the glycan name as a slot or attribute
    glycan_name <- if (!is.null(result@glycan_name)) {
      as.character(result@glycan_name)
    } else if (!is.null(attr(result, "glycan_name"))) {
      as.character(attr(result, "glycan_name"))
    } else {
      "Unknown_Glycan"
    }
    
    # Filter only significant pathways (p.adjust < 0.05)
    significant_results <- result@result[result@result$p.adjust < 0.05, ]
    
    if (nrow(significant_results) > 0) {
      # Append glycan name to each description
      descriptions <- paste0(significant_results$Description, "_", glycan_name)
      gene_lists <- strsplit(significant_results$geneID, "/")
      
      # Store them in the cumulative lists
      all_descriptions <- c(all_descriptions, descriptions)
      all_gene_lists <- c(all_gene_lists, gene_lists)
    }
  }
}

# Ensure that all descriptions and gene lists have the same length
n <- length(all_descriptions)

if (n != length(all_gene_lists)) {
  stop("Mismatch between descriptions and gene lists length.")
}

# Initialize the overlap matrix
overlap_matrix <- matrix(0, nrow = n, ncol = n, 
                         dimnames = list(all_descriptions, all_descriptions))

# Initialize the progress bar
pb <- progress_bar$new(
  format = "  Calculating overlap [:bar] :percent in :elapsed",
  total = (n * (n - 1)) / 2,  # Total iterations for upper triangle matrix excluding diagonal
  clear = FALSE, width = 60
)

# Fill the overlap matrix with percentage overlaps
for (i in 1:(n - 1)) {  # Avoid the last row, as it has no upper triangle elements
  for (j in (i + 1):n) {  # Only calculate the upper triangle
    # Calculate intersection and union of unique genes for descriptions i and j
    intersection <- length(intersect(all_gene_lists[[i]], all_gene_lists[[j]]))
    union <- length(unique(c(all_gene_lists[[i]], all_gene_lists[[j]])))
    
    # Calculate percentage overlap
    overlap_matrix[i, j] <- (intersection / union) * 100
    overlap_matrix[j, i] <- overlap_matrix[i, j]  # Symmetric matrix
    
    # Update the progress bar
    pb$tick()
  }
}
# View the overlap matrix
print(overlap_matrix)

# Check for any NA values
any_na <- any(is.na(overlap_matrix))
cat("\nAny NA values in overlap matrix:", any_na, "\n")
saveRDS(overlap_matrix, file = "your_path/overlap_matrix_fig_5.rds")
overlap_matrix <- readRDS("your_path/overlap_matrix_fig_5.rds")

# Load necessary libraries
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(readr)

# List of glycans to include
included_glycans <- c("fucosylated", "sialylated", "triantennary", "agalactosylated")

# Initialize a list to store the top GO terms for all glycans
all_top_go <- list()

cat("\nExtracting top GO terms...\n")

# Iterate through the GO results and extract top 10 GO terms for each glycan
for (i in seq_along(go_results)) {
  cat("\nProcessing GO result:", i, "\n")  # Debug message to track progress
  result <- go_results[[i]]
  
  if (!is.null(result) && nrow(result@result) > 0) {
    tryCatch({
      # Extract glycan name from the result's attribute
      glycan_name <- attr(result, "glycan_name")
      if (is.null(glycan_name)) glycan_name <- paste("Unknown Glycan", i)
      
      if (!(glycan_name %in% included_glycans)) {
        cat("Skipping glycan:", glycan_name, "\n")
        next
      }
      
      # Extract the original GeneRatio from the result
      cat("Extracting GeneRatio...\n")
      top_terms <- result@result %>%
        arrange(pvalue) %>%
        slice_head(n = 15) %>%
        mutate(
          glycan = glycan_name,
          gene_count = lengths(strsplit(geneID, "/")),  # Calculate gene counts
          GeneRatio = GeneRatio, # Use existing GeneRatio from the data
          Description = paste0(Description, "_", glycan_name)  # Append glycan name to Description
        )
      
      # Store the top terms for further merging
      all_top_go[[glycan_name]] <- top_terms
      cat("Top GO terms extracted for glycan:", glycan_name, "\n")
      
    }, error = function(e) {
      cat("Error processing glycan:", glycan_name, "-", e$message, "\n")
    })
  } else {
    cat("Skipping GO result:", i, "as it is NULL or has no enriched terms.\n")
  }
}

# Combine all top GO terms into a single dataframe
cat("\nCombining all top GO terms into a single dataframe...\n")
all_top_go_df <- bind_rows(all_top_go)
all_top_go_df

all_top_go_df <- all_top_go_df %>%
  mutate(
    # Extract the denominator from GeneRatio and convert it to a numeric column
    total_genes_in_pathway = as.numeric(sub(".*/", "", GeneRatio)),
    
    # Calculate Gratio by dividing Count by total_genes_in_pathway
    Gratio = Count / total_genes_in_pathway
  )



# Get the list of Description values from all_top_go_df
descriptions <- all_top_go_df$Description

# Filter overlap_matrix to include only rows and columns that match the Description values in all_top_go_df
filtered_overlap_matrix <- overlap_matrix[rownames(overlap_matrix) %in% descriptions, colnames(overlap_matrix) %in% descriptions]

print(filtered_overlap_matrix)

saveRDS(filtered_overlap_matrix, file = "your_path/filtered_overlap_matrix_fig_5.rds")

filtered_overlap_matrix <- readRDS("your_pathg/filtered_overlap_matrix_fig_5.rds")

# RNA - cell correlations----

RNA_cols <- unique(master_feature_table_filtered$feature_variable[master_feature_table_filtered$bio_feature_type == "spatial_RNA"])

RNA_cols


immune_cells <- c("Myeloid_CD14", 
                  "APC", "Tcell_CD8", 
                  "Tcell_FoxP3", "Tcell_CD4", 
                  "Microglia", "Macrophage_CD68", 
                  "Macrophage_CD206", "Myeloid_CD141", "Myeloid_CD11b_HLADRminus", 
                  "Myeloid_CD11b_HLADRplus", "Neutrophils", 
                  "DC_Mac_CD209", "Bcells", "Mast_cells", 
                  "Myeloid_CD14_CD163","Macrophage_CD68_CD163", 
                  "Microglia_CD163")

other_cells <- c("Tumor_over_all_cell_count_tumor_FOV_prop", "Endothelial_over_all_cell_count_FOV_immune_prop","Neurons_over_all_cell_count_FOV_immune_prop")


# Append "_over_all_immune_count_prop" to each element in immune_cells
immune_cell_types <- paste0(immune_cells, "_over_all_cell_count_FOV_immune_prop")

combined_list <- c(immune_cell_types, other_cells)

filtered_master_feature_table <- master_feature_table_filtered %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"))

filtered_master_feature_table <- filtered_master_feature_table %>%
  filter((feature_variable %in% combined_list |feature_variable %in% RNA_cols )
  )

# Pivot the filtered master feature table to have features as columns
master_feature_pivot <- filtered_master_feature_table %>%
  dplyr::select(patient_id, feature_variable, feature_value) %>%
  pivot_wider(names_from = feature_variable, values_from = feature_value, values_fn = mean) %>%
  dplyr::select(patient_id, everything())



# Start the timer
tic("Correlation computation")

# Initialize progress bar
num_steps <- length(RNA_cols) # Track progress per RNA column
pb <- progress_bar$new(
  format = "Computing [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = num_steps, clear = FALSE, width = 60
)

# Initialize an empty matrix to store correlations
correlation_matrix <- matrix(NA, nrow = length(RNA_cols), ncol = length(combined_list))
rownames(correlation_matrix) <- RNA_cols
colnames(correlation_matrix) <- combined_list

# Perform correlation with progress tracking
for (i in seq_along(RNA_cols)) {
  for (j in seq_along(combined_list)) {
    correlation_matrix[i, j] <- cor(
      master_feature_pivot[[RNA_cols[i]]],
      master_feature_pivot[[combined_list[j]]],
      use = "pairwise.complete.obs",
      method = "pearson"
    )
  }
  pb$tick()  # Update progress bar
}

# Stop the timer
toc()

# View the correlation matrix
print(correlation_matrix)

# prep for cell - glycan - RNA correlation analysis ----


# Initialize an empty dataframe to store results for all glycans
merged_all_top_go_df <- data.frame()

# Identify groups of correlated descriptions across all glycans using >80% threshold
adj_matrix <- (filtered_overlap_matrix > 70) * 1  # Convert to binary matrix
graph <- graph_from_adjacency_matrix(as.matrix(adj_matrix), mode = "undirected", diag = FALSE)
clusters <- clusters(graph)$membership


# List of glycans to include
included_glycans <- c("fucosylated", "sialylated", "triantennary", "agalactosylated")




# Loop through each glycan and merge descriptions based on correlation clusters
for (glycan in included_glycans) {
  # Filter `all_top_go_df` for the current glycan
  glycan_df <- all_top_go_df %>%
    filter(endsWith(Description, paste0("_", glycan)))
  
  # Loop over clusters to create merged descriptions
  glycan_merged <- list()
  for (cluster_id in unique(clusters)) {
    cluster_descriptions <- names(clusters)[clusters == cluster_id]
    merged_name <- paste(cluster_descriptions, collapse = " _ ")  # Merged name including all glycans
    
    # Filter glycan-specific entries in `cluster_descriptions`
    glycan_specific_descriptions <- cluster_descriptions[endsWith(cluster_descriptions, paste0("_", glycan))]
    
    # Calculate mean values for glycan-specific entries in this cluster
    if (length(glycan_specific_descriptions) > 0) {
      merged_row <- glycan_df %>%
        filter(Description %in% glycan_specific_descriptions) %>%
        summarise(
          glycan = glycan,
          pvalue = mean(pvalue, na.rm = TRUE),
          p.adjust = mean(p.adjust, na.rm = TRUE),
          qvalue = mean(qvalue, na.rm = TRUE),
          Count = mean(Count, na.rm = TRUE),
          gene_count = mean(gene_count, na.rm = TRUE),
          Gratio = mean(Gratio, na.rm = TRUE),
          Original_Descriptions = paste(glycan_specific_descriptions, collapse = " / ")
        ) %>%
        mutate(Description = merged_name)  # Use the combined description name
      
      glycan_merged[[length(glycan_merged) + 1]] <- merged_row
    }
  }
  
  # Combine all merged results for this glycan and add to the main dataframe
  merged_all_top_go_df <- bind_rows(merged_all_top_go_df, bind_rows(glycan_merged))
}

# Display the final merged dataframe with results from all glycans
print(merged_all_top_go_df)

cat("\nAll glycans processed and merged with cross-glycan correlations.\n")

# Remove duplicates from the merged results if any
merged_all_top_go_df <- merged_all_top_go_df %>%
  distinct(Description, .keep_all = TRUE)

# Display the final merged dataframe with duplicates removed
print(merged_all_top_go_df)

cat("\nDuplicates removed. Merging complete.\n")


# Save the dataframe as an .RData file
saveRDS(merged_all_top_go_df, file = "your_path/merged_all_top_go_df_fig_5.rds")

