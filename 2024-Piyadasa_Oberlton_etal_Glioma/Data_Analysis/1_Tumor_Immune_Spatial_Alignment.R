library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(tidyr)
library(data.table)
library(arrow)
# change to be base notebook
data_Immune <- fread("/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/immune_panel/20240919_immune_no_controls_cell_table_size_norm_metacluser_ID_and_cell_labels_and_func_with_remaining_final.csv")
data_Tumor <- fread("/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/tumor_panel/20240503_cell_table_size_norm_metacluser_ID_and_cell_labels_post_clustering_with_rem.csv")

# Data Type Alignment -----
shared_columns <- intersect(names(data_Tumor), names(data_Immune))

# Loop through the shared columns and change data types
for (col in shared_columns) {
  if (class(data_Tumor[[col]]) != class(data_Immune[[col]])) {
    # Detect the class of the column in data_Stanford_Tumor
    target_class <- class(data_Tumor[[col]])
    
    # Convert the column in data_Stanford_Immune to match
    if (target_class == "factor") {
      data_Immune[[col]] <- as.factor(data_Immune[[col]])
    } else if (target_class == "numeric") {
      data_Immune[[col]] <- as.numeric(as.character(data_Immune[[col]]))
    } else if (target_class == "integer") {
      data_Immune[[col]] <- as.integer(as.character(data_Immune[[col]]))
    } else if (target_class == "character") {
      data_Immune[[col]] <- as.character(data_Immune[[col]])
    }
  }
}
# Create centroid-1_test and centroid-2_test -----------
# Centroid-1_test and centroid-2_test will act as editable copies 

# Immune
data_Immune[, `:=` (
  `centroid-0_test` = `centroid-0`,
  `centroid-0_nuclear_test` = `centroid-0_nuclear`,
  `centroid-1_test` = `centroid-1`,
  `centroid-1_nuclear_test` = `centroid-1_nuclear`,
  Panel = "Immune"
)]

# Tumor
data_Tumor[, `:=` (
  `centroid-0_test` = `centroid-0`,
  `centroid-0_nuclear_test` = `centroid-0_nuclear`,
  `centroid-1_test` = `centroid-1`,
  `centroid-1_nuclear_test` = `centroid-1_nuclear`,
  Panel = "Tumor"
)]

# Create_and_plot_data----
create_and_plot_data <- function(working_FOV, image_size, centroid_0_change, centroid_1_change, 
                                 Immune_meta_cluster, Immune_cell_type, 
                                 Tumor_meta_cluster, Tumor_cell_type) {
  
  # move image
  # Apply changes to centroid values
  rows_to_change_Immune <- data_Immune$fov == working_FOV
  data_Immune[rows_to_change_Immune, `:=` (
    `centroid-0_test` = `centroid-0` + centroid_0_change,
    `centroid-0_nuclear_test` = `centroid-0_nuclear` + centroid_0_change,
    `centroid-1_test` = `centroid-1` + centroid_1_change,
    `centroid-1_nuclear_test` = `centroid-1_nuclear` + centroid_1_change
  )]
  
  # Filter data_Immune and data_Tumor for fov == 'CHOP_907_R11C1'
  filtered_Immune <- data_Immune[fov == working_FOV]
  filtered_Tumor <- data_Tumor[fov == working_FOV]
  
  # Further subset if meta cluster and cell type are specified
  if (!is.null(Immune_meta_cluster) && !is.null(Immune_cell_type)) {
    filtered_Immune <- filtered_Immune[get(Immune_meta_cluster) == Immune_cell_type]
  }
  
  if (!is.null(Tumor_meta_cluster) && !is.null(Tumor_cell_type)) {
    filtered_Tumor <- filtered_Tumor[get(Tumor_meta_cluster) == Tumor_cell_type]
  }
  
  # Combine the filtered data
  combined_data <- rbindlist(list(filtered_Immune, filtered_Tumor), fill = TRUE)
  
  # Add a new column to identify the source
  combined_data[, source := ifelse(.I <= nrow(filtered_Immune), 'Immune', 'Tumor')]
  
  # Plotting
  
  ### INDIVIDUAL PLOTS ###
  immune_plot <- ggplot(filtered_Immune, aes(x = `centroid-1`, y = `centroid-0`)) +
    geom_point(color = "blue") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = paste(working_FOV, "Immune"))
  
  # Plot for filtered_Tumor
  tumor_plot <- ggplot(filtered_Tumor, aes(x = `centroid-1`, y = `centroid-0`)) +
    geom_point(color = "red") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = paste(working_FOV, "Tumor"))
  
  # Display the plots
  print(immune_plot)
  print(tumor_plot)
  
  ### INDIVIDUAL PLOTS END ###
  
  # Plot with white background - Original points
  pre_change_plot <- ggplot(combined_data, aes(x = `centroid-1`, y = `centroid-0`, color = source)) +
    geom_point() +
    scale_color_manual(values = c("Immune" = "blue", "Tumor" = "red")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = working_FOV)
  
  # Plot with white background - Adjusted points
  post_change_plot <- ggplot(combined_data, aes(x = `centroid-1_test`, y = `centroid-0_test`, color = source)) +
    geom_point() +
    scale_color_manual(values = c("Immune" = "blue", "Tumor" = "red")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = working_FOV)
  
  
  print(pre_change_plot)
  print(post_change_plot) 
  
  assign("data_Immune", data_Immune, envir = .GlobalEnv)
  assign("data_Tumor", data_Tumor, envir = .GlobalEnv)
  
  return(list(pre_change_plot = pre_change_plot, post_change_plot = post_change_plot))
}


# Create_and plot_data_separate_FOVs----
create_and_plot_data_separate_FOVs <- function(working_FOV_Immune, working_FOV_Tumor, image_size, centroid_0_change, centroid_1_change, 
                                 Immune_meta_cluster, Immune_cell_type, 
                                 Tumor_meta_cluster, Tumor_cell_type) {
  
  # move image
  # Apply changes to centroid values
  rows_to_change_Immune <- data_Immune$fov == working_FOV_Immune
  data_Immune[rows_to_change_Immune, `:=` (
    `centroid-0_test` = `centroid-0` + centroid_0_change,
    `centroid-0_nuclear_test` = `centroid-0_nuclear` + centroid_0_change,
    `centroid-1_test` = `centroid-1` + centroid_1_change,
    `centroid-1_nuclear_test` = `centroid-1_nuclear` + centroid_1_change
  )]
  
  # Filter data_Immune and data_Tumor for fov 
  filtered_Immune <- data_Immune[fov == working_FOV_Immune]
  filtered_Tumor <- data_Tumor[fov == working_FOV_Tumor]
  
  # Further subset if meta cluster and cell type are specified
  if (!is.null(Immune_meta_cluster) && !is.null(Immune_cell_type)) {
    filtered_Immune <- filtered_Immune[get(Immune_meta_cluster) == Immune_cell_type]
  }
  
  if (!is.null(Tumor_meta_cluster) && !is.null(Tumor_cell_type)) {
    filtered_Tumor <- filtered_Tumor[get(Tumor_meta_cluster) == Tumor_cell_type]
  }
  
  # Combine the filtered data
  combined_data <- rbindlist(list(filtered_Immune, filtered_Tumor), fill = TRUE)
  
  # Add a new column to identify the source
  combined_data[, source := ifelse(.I <= nrow(filtered_Immune), 'Immune', 'Tumor')]
  
  # Plotting
  
  ### INDIVIDUAL PLOTS ###
  immune_plot <- ggplot(filtered_Immune, aes(x = `centroid-1`, y = `centroid-0`)) +
    geom_point(color = "blue") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = paste(working_FOV_Immune, "Immune"))
  
  # Plot for filtered_Tumor
  tumor_plot <- ggplot(filtered_Tumor, aes(x = `centroid-1`, y = `centroid-0`)) +
    geom_point(color = "red") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = paste(working_FOV_Immune, "Tumor"))
  
  # Display the plots
  print(immune_plot)
  print(tumor_plot)
  
  ### INDIVIDUAL PLOTS END ###
  
  # Plot with white background - Original points
  pre_change_plot <- ggplot(combined_data, aes(x = `centroid-1`, y = `centroid-0`, color = source)) +
    geom_point() +
    scale_color_manual(values = c("Immune" = "blue", "Tumor" = "red")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = working_FOV_Immune)
  
  # Plot with white background - Adjusted points
  post_change_plot <- ggplot(combined_data, aes(x = `centroid-1_test`, y = `centroid-0_test`, color = source)) +
    geom_point() +
    scale_color_manual(values = c("Immune" = "blue", "Tumor" = "red")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = working_FOV_Immune)
  
  
  print(pre_change_plot)
  print(post_change_plot) 
  
  assign("data_Immune", data_Immune, envir = .GlobalEnv)
  assign("data_Tumor", data_Tumor, envir = .GlobalEnv)
  
  
  return(list(pre_change_plot = pre_change_plot, post_change_plot = post_change_plot))
}
# Apply_and_save_changes----
apply_changes_and_save <- function(working_FOV, image_size, centroid_0_change, centroid_1_change, 
                                   Immune_meta_cluster, Immune_cell_type, 
                                   Tumor_meta_cluster, Tumor_cell_type, 
                                   pre_change_plot, post_change_plot) {
  

  
  ### SAVE PLOTS ###
  # Define the directory path for the new folder
  new_folder_path <- file.path("/Users/boberlto/Desktop/Research/Bendall Lab/analysis pipeline/final_analysis/Alignment_Images", working_FOV)
  
  # Create the new folder if it doesn't exist
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  # Define file paths for the plots
  pre_change_plot_path_png <- file.path(new_folder_path, "pre_change_plot.png")
  pre_change_plot_path_tiff <- file.path(new_folder_path, "pre_change_plot.tiff")
  post_change_plot_path_png <- file.path(new_folder_path, "post_change_plot.png")
  post_change_plot_path_tiff <- file.path(new_folder_path, "post_change_plot.tiff")
  
  # Save the plots
  ggsave(pre_change_plot_path_png, pre_change_plot, width = 10, height = 8)
  ggsave(pre_change_plot_path_tiff, pre_change_plot, width = 10, height = 8, dpi = 100)
  ggsave(post_change_plot_path_png, post_change_plot, width = 10, height = 8)
  ggsave(post_change_plot_path_tiff, post_change_plot, width = 10, height = 8, dpi = 100)
  ### SAVE PLOTS END ###
  
  ###
  
  ### APPLY CHANGES ###\
  rows_to_change_Immune <- data_Immune$fov == working_FOV
  # Apply changes to centroid values
  data_Immune[rows_to_change_Immune, `:=` (
    `centroid-0` = `centroid-0_test`,
    `centroid-0_nuclear` = `centroid-0_nuclear_test`,
    `centroid-1` = `centroid-1_test`,
    `centroid-1_nuclear` = `centroid-1_nuclear_test`
  )]
  ### APPLY CHANGES END###
  
  ###
  
  ### REMOVE POINTS OUTSIDE OF BOUNDS###
  
  # Function to filter and identify rows to delete
  min_centroid_0 = max(0, 0 + centroid_0_change)
  min_centroid_1 = max(0, 0 + centroid_1_change)
  max_centroid_0 = min(image_size, image_size + centroid_0_change)
  max_centroid_1 = min(image_size, image_size + centroid_1_change)
  
  delete_rows <- function(df) {
    # Identify rows where centroid_0 or centroid_1 are outside the specified ranges
    rows_to_delete <- df$fov == working_FOV & (
      df$`centroid-0` < min_centroid_0 | 
        df$`centroid-0` > max_centroid_0 |
        df$`centroid-1` < min_centroid_1 | 
        df$`centroid-1` > max_centroid_1
    )
    
    # Delete rows
    df <- df[!rows_to_delete, ]
    
    return(df)
  }

  # delete rows outside of constraints
  data_Immune <- delete_rows(data_Immune)
  data_Tumor <- delete_rows(data_Tumor)
  assign("data_Immune", data_Immune, envir = .GlobalEnv)
  assign("data_Tumor", data_Tumor, envir = .GlobalEnv)
  
  ### REMOVE POINTS OUTSIDE OF BOUNDS END###
  
  ###
  
  ### FINAL PLOT ###
  
  # Filter data_Immune and data_Tumor for fov 
  filtered_Immune <- data_Immune[fov == working_FOV]
  filtered_Tumor <- data_Tumor[fov == working_FOV]
  
  # Further subset if meta cluster and cell type are specified
  if (!is.null(Immune_meta_cluster) && !is.null(Immune_cell_type)) {
    filtered_Immune <- filtered_Immune[get(Immune_meta_cluster) == Immune_cell_type]
  }
  
  if (!is.null(Tumor_meta_cluster) && !is.null(Tumor_cell_type)) {
    filtered_Tumor <- filtered_Tumor[get(Tumor_meta_cluster) == Tumor_cell_type]
  }
  
  # Combine the filtered data
  combined_data <- rbindlist(list(filtered_Immune, filtered_Tumor), fill = TRUE)
  
  # Add a new column to identify the source
  combined_data[, source := ifelse(.I <= nrow(filtered_Immune), 'Immune', 'Tumor')]
  
  # Plotting
  
  # Plot with white background - Original points
  final_plot <- ggplot(combined_data, aes(x = `centroid-1`, y = `centroid-0`, color = source)) +
    geom_point() +
    scale_color_manual(values = c("Immune" = "blue", "Tumor" = "red")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = working_FOV)
  
  
  print(final_plot)
  ### FINAL PLOT END ###
  
  ###
  
  ### SAVE FINAL PLOT ###
  
  # Create the new folder if it doesn't exist
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  # Define file paths for the plots
  final_plot_path_png <- file.path(new_folder_path, "final_plot.png")
  final_plot_path_tiff <- file.path(new_folder_path, "final_plot.tiff")
  
  # Save the plots
  ggsave(final_plot_path_png, final_plot, width = 10, height = 8)
  ggsave(final_plot_path_tiff, final_plot, width = 10, height = 8, dpi = 100)
}

# Apply_and_save_changes_separate_FOVs----
apply_changes_and_save_separate_FOVs <- function(working_FOV_Immune, working_FOV_Tumor, image_size, centroid_0_change, centroid_1_change, 
                                   Immune_meta_cluster, Immune_cell_type, 
                                   Tumor_meta_cluster, Tumor_cell_type, 
                                   pre_change_plot, post_change_plot) {
  
  
  
  ### SAVE PLOTS ###
  # Define the directory path for the new folder
  new_folder_path <- file.path("/Users/boberlto/Desktop/Research/Bendall Lab/analysis pipeline/final_analysis/Alignment_Images", working_FOV_Immune)
  
  # Create the new folder if it doesn't exist
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  # Define file paths for the plots
  pre_change_plot_path_png <- file.path(new_folder_path, "pre_change_plot.png")
  pre_change_plot_path_tiff <- file.path(new_folder_path, "pre_change_plot.tiff")
  post_change_plot_path_png <- file.path(new_folder_path, "post_change_plot.png")
  post_change_plot_path_tiff <- file.path(new_folder_path, "post_change_plot.tiff")
  
  # Save the plots
  ggsave(pre_change_plot_path_png, pre_change_plot, width = 10, height = 8)
  ggsave(pre_change_plot_path_tiff, pre_change_plot, width = 10, height = 8, dpi = 100)
  ggsave(post_change_plot_path_png, post_change_plot, width = 10, height = 8)
  ggsave(post_change_plot_path_tiff, post_change_plot, width = 10, height = 8, dpi = 100)
  ### SAVE PLOTS END ###
  
  ###
  
  ### APPLY CHANGES ###\
  rows_to_change_Immune <- data_Immune$fov == working_FOV_Immune
  # Apply changes to centroid values
  data_Immune[rows_to_change_Immune, `:=` (
    `centroid-0` = `centroid-0_test`,
    `centroid-0_nuclear` = `centroid-0_nuclear_test`,
    `centroid-1` = `centroid-1_test`,
    `centroid-1_nuclear` = `centroid-1_nuclear_test`
  )]
  ### APPLY CHANGES END###
  
  ###
  
  ### REMOVE POINTS OUTSIDE OF BOUNDS###
  
  # Function to filter and identify rows to delete
  min_centroid_0 = max(0, 0 + centroid_0_change)
  min_centroid_1 = max(0, 0 + centroid_1_change)
  max_centroid_0 = min(image_size, image_size + centroid_0_change)
  max_centroid_1 = min(image_size, image_size + centroid_1_change)
  
  delete_rows <- function(df, unique_FOV) {
    # Identify rows where centroid_0 or centroid_1 are outside the specified ranges
    rows_to_delete <- df$fov == unique_FOV & (
      df$`centroid-0` < min_centroid_0 | 
        df$`centroid-0` > max_centroid_0 |
        df$`centroid-1` < min_centroid_1 | 
        df$`centroid-1` > max_centroid_1
    )
    
    # Delete rows
    df <- df[!rows_to_delete, ]
    
    return(df)
  }
  
  # delete rows outside of constraints
  data_Immune <- delete_rows(data_Immune, working_FOV_Immune)
  data_Tumor <- delete_rows(data_Tumor, working_FOV_Tumor)
  assign("data_Immune", data_Immune, envir = .GlobalEnv)
  assign("data_Tumor", data_Tumor, envir = .GlobalEnv)
  
  ### REMOVE POINTS OUTSIDE OF BOUNDS END###
  
  ###
  
  ### FINAL PLOT ###
  
  # Filter data_Immune and data_Tumor for fov 
  filtered_Immune <- data_Immune[fov == working_FOV_Immune]
  filtered_Tumor <- data_Tumor[fov == working_FOV_Tumor]
  
  # Further subset if meta cluster and cell type are specified
  if (!is.null(Immune_meta_cluster) && !is.null(Immune_cell_type)) {
    filtered_Immune <- filtered_Immune[get(Immune_meta_cluster) == Immune_cell_type]
  }
  
  if (!is.null(Tumor_meta_cluster) && !is.null(Tumor_cell_type)) {
    filtered_Tumor <- filtered_Tumor[get(Tumor_meta_cluster) == Tumor_cell_type]
  }
  
  # Combine the filtered data
  combined_data <- rbindlist(list(filtered_Immune, filtered_Tumor), fill = TRUE)
  
  # Add a new column to identify the source
  combined_data[, source := ifelse(.I <= nrow(filtered_Immune), 'Immune', 'Tumor')]
  
  # Plotting
  
  # Plot with white background - Original points
  final_plot <- ggplot(combined_data, aes(x = `centroid-1`, y = `centroid-0`, color = source)) +
    geom_point() +
    scale_color_manual(values = c("Immune" = "blue", "Tumor" = "red")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_y_reverse() +
    labs(x = "Centroid-1", y = "Centroid-0", title = working_FOV_Immune)
  
  
  print(final_plot)
  ### FINAL PLOT END ###
  
  ###
  
  ### SAVE FINAL PLOT ###
  
  # Create the new folder if it doesn't exist
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  # Define file paths for the plots
  final_plot_path_png <- file.path(new_folder_path, "final_plot.png")
  final_plot_path_tiff <- file.path(new_folder_path, "final_plot.tiff")
  
  # Save the plots
  ggsave(final_plot_path_png, final_plot, width = 10, height = 8)
  ggsave(final_plot_path_tiff, final_plot, width = 10, height = 8, dpi = 100)
}

# 'CHOP_907_R11C1'----
plots <- create_and_plot_data("CHOP_907_R11C1", 1024, -100, 20, NULL, NULL, NULL, NULL)
apply_changes_and_save("CHOP_907_R11C1", 1024, -100, 20, NULL, NULL, NULL, NULL, 
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R13C14'----
plots <- create_and_plot_data("CHOP_907_R13C14", 1024, -100, 0,"cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune")
apply_changes_and_save("CHOP_907_R13C14", 1024, -100, 0, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune",
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R13C7'----
plots <- create_and_plot_data("CHOP_907_R13C7", 1024, -120, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save("CHOP_907_R13C7", 1024, -120, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R1C11'----
plots <- create_and_plot_data("CHOP_907_R1C11", 1024, 110, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save("CHOP_907_R1C11", 1024, 110, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R3C3'----
plots <- create_and_plot_data('CHOP_907_R3C3', 1024, -50, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R3C3', 1024, -50, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R4C7'----
plots <- create_and_plot_data('CHOP_907_R4C7', 1024, -110, 30, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R4C7', 1024, -110, 30, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R5C11'----
plots <- create_and_plot_data('CHOP_907_R5C11', 1024, -30, 200, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R5C11', 1024, -30, 200, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R5C4'----
plots <- create_and_plot_data('CHOP_907_R5C4', 1024, -50, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R5C4', 1024, -50, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'CHOP_907_R6C2'----
plots <- create_and_plot_data('CHOP_907_R6C2', 1024, -125, -20, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R6C2', 1024, -125, -20, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'CHOP_907_R6C4'----
plots <- create_and_plot_data('CHOP_907_R6C4', 1024, -100, -60, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R6C4', 1024, -100, -60, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'CHOP_907_R6C5'----
plots <- create_and_plot_data('CHOP_907_R6C5', 1024, -120, -20, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_907_R6C5', 1024, -120, -20, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'CHOP_941_R4C1'----
plots <- create_and_plot_data('CHOP_941_R4C1', 1024, 0, 50, NULL, NULL, NULL, NULL)
apply_changes_and_save('CHOP_941_R4C1', 1024, 0, 50, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'CoH_C001_R2C1'----
plots <- create_and_plot_data('CoH_C001_R2C1', 2048, -200, -80, NULL, NULL, NULL, NULL)
apply_changes_and_save('CoH_C001_R2C1', 2048, -200, -80, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R11C2'----
plots <- create_and_plot_data('Stanford_TA550_R11C2', 2048, 200, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R11C2', 2048, 200, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R11C3'----
plots <- create_and_plot_data('Stanford_TA550_R11C3', 2048, 150, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R11C3', 2048, 150, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R1C3'----
plots <- create_and_plot_data('Stanford_TA550_R1C3', 2048, -80, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R1C3', 2048, -80, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R2C2' ----
plots <- create_and_plot_data('Stanford_TA550_R2C2', 2048, -500, -200, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R2C2', 2048, -500, -200, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'Stanford_TA550_R2C4' ----
plots <- create_and_plot_data('Stanford_TA550_R2C4', 1024, 0, -100, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R2C4', 1024, 0, -100, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'Stanford_TA550_R3C1'----
plots <- create_and_plot_data('Stanford_TA550_R3C1', 2048, 130, 320, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R3C1', 2048, 130, 320, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R3C3'----
plots <- create_and_plot_data('Stanford_TA550_R3C3', 2048, -140, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R3C3', 2048, -140, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R6C3'----
plots <- create_and_plot_data('Stanford_TA550_R6C3', 2048, -150, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R6C3', 2048, -150, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R7C1'----
plots <- create_and_plot_data('Stanford_TA550_R7C1', 2048, 80, 150, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R7C1', 2048, 80, 150, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R7C2'----
plots <- create_and_plot_data('Stanford_TA550_R7C2', 2048, 0, 100, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R7C2', 2048, 0, 100, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA550_R7C3'----
plots <- create_and_plot_data('Stanford_TA550_R7C3', 2048, -50, 50, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA550_R7C3', 2048, -50, 50, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA551_R2C1'----
plots <- create_and_plot_data('Stanford_TA551_R2C1', 2048, -180, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA551_R2C1', 2048, -180, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA551_R3C1'----
plots <- create_and_plot_data('Stanford_TA551_R3C1', 1024, 70, 5, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA551_R3C1', 1024, 70, 5, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'Stanford_TA551_R3C2'----
plots <- create_and_plot_data('Stanford_TA551_R3C2', 2048, 110, 70, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA551_R3C2', 2048, 110, 70, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA551_R5C2'----
plots <- create_and_plot_data('Stanford_TA551_R5C2', 1024, 90, 50, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA551_R5C2', 1024, 90, 50, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R10C1'----
plots <- create_and_plot_data('Stanford_TA552_R10C1', 2048, 800, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R10C1', 2048, 800, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R10C4'----
plots <- create_and_plot_data('Stanford_TA552_R10C4', 2048, 550, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R10C4', 2048, 550, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R3C1'----
plots <- create_and_plot_data('Stanford_TA552_R3C1', 2048, 450, 15, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R3C1', 2048, 450, 15, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R3C2'----
plots <- create_and_plot_data('Stanford_TA552_R3C2', 2048, 520, -10, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R3C2', 2048, 520, -10, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R3C4'----
plots <- create_and_plot_data('Stanford_TA552_R3C4', 2048, 510, -30, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R3C4', 2048, 510, -30, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R4C1'----
plots <- create_and_plot_data('Stanford_TA552_R4C1', 2048, 570, 80, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R4C1', 2048, 570, 80, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R4C4'----
plots <- create_and_plot_data('Stanford_TA552_R4C4', 2048, 530, -10, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R4C4', 2048, 530, -10, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R5C1'----
plots <- create_and_plot_data('Stanford_TA552_R5C1', 2048, 530, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R5C1', 2048, 530, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R6C4'----
plots <- create_and_plot_data('Stanford_TA552_R6C4', 2048, 500, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R6C4', 2048, 500, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R7C1'----
plots <- create_and_plot_data('Stanford_TA552_R7C1', 2048, 500, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R7C1', 2048, 500, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)


# 'Stanford_TA552_R7C2'----
plots <- create_and_plot_data('Stanford_TA552_R7C2', 2048, 500, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R7C2', 2048, 500, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R7C4'----
plots <- create_and_plot_data('Stanford_TA552_R7C4', 2048, 520, -40, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R7C4', 2048, 520, -40, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R7C5'----
plots <- create_and_plot_data('Stanford_TA552_R7C5', 2048, 480, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R7C5', 2048, 480, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R8C1'----
plots <- create_and_plot_data('Stanford_TA552_R8C1', 2048, 480, 60, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R8C1', 2048, 480, 60, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R8C2'----
plots <- create_and_plot_data('Stanford_TA552_R8C2', 2048, 480, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R8C2', 2048, 480, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R8C4'----
plots <- create_and_plot_data('Stanford_TA552_R8C4', 2048, 570, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R8C4', 2048, 570, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R9C2'----
plots <- create_and_plot_data('Stanford_TA552_R9C2', 2048, 500, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R9C2', 2048, 500, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R9C4'----
plots <- create_and_plot_data('Stanford_TA552_R9C4', 2048, 500, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R9C4', 2048, 500, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'Stanford_TA552_R9C5'----
plots <- create_and_plot_data('Stanford_TA552_R9C5', 2048, 515, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('Stanford_TA552_R9C5', 2048, 515, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCLA_R1C1'----
plots <- create_and_plot_data('UCLA_R1C1', 2048, 230, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCLA_R1C1', 2048, 230, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCLA_R1C5'----
plots <- create_and_plot_data('UCLA_R1C5', 2048, 100, 35, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCLA_R1C5', 2048, 100, 35, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCLA_R6C3'----
plots <- create_and_plot_data('UCLA_R6C3', 2048, 100, 50, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCLA_R6C3', 2048, 100, 50, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCLA_R6C4'----
plots <- create_and_plot_data('UCLA_R6C4', 2048, 100, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCLA_R6C4', 2048, 100, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCLA_R6C5'----
plots <- create_and_plot_data('UCLA_R6C5', 1024, 100, 150, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCLA_R6C5', 1024, 100, 150, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCLA_R7C5'----
plots <- create_and_plot_data('UCLA_R7C5', 2048, 100, 50, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCLA_R7C5', 2048, 100, 50, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)
# 'UCLA_R98C3' (different) ----
plots <- create_and_plot_data_separate_FOVs('UCLA_R98C3','UCLA_R99C3', 2048, -50, -175, NULL, NULL, NULL, NULL)
apply_changes_and_save_separate_FOVs('UCLA_R98C3','UCLA_R99C3', 2048, -50, -175, NULL, NULL, NULL, NULL,
                                     plots$pre_change_plot, plots$post_change_plot)

# 'UCSF_TA2_R2C5'----
plots <- create_and_plot_data('UCSF_TA2_R2C5', 2048, 250, -70, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune")
apply_changes_and_save('UCSF_TA2_R2C5', 2048, 250, -70, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune",
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA2_R3C5'----
plots <- create_and_plot_data('UCSF_TA2_R3C5', 2048, -50, 240, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA2_R3C5', 2048, -50, 240, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA2_R4C5'----
plots <- create_and_plot_data('UCSF_TA2_R4C5', 2048, -140, 0, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA2_R4C5', 2048, -140, 0, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA3_R4C4'----
plots <- create_and_plot_data('UCSF_TA3_R4C4', 2048, 70, -100, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA3_R4C4', 2048, 70, -100, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA3_R6C5'----
plots <- create_and_plot_data('UCSF_TA3_R6C5', 2048, 100, -80, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA3_R6C5', 2048, 100, -80, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA4_R2C5'----
plots <- create_and_plot_data('UCSF_TA4_R2C5', 2048, 0, -150, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA4_R2C5', 2048, 0, -150, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA4_R6C1'----
plots <- create_and_plot_data('UCSF_TA4_R6C1', 2048, 120, -160, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA4_R6C1', 2048, 120, -160, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA5_R4C2'----
plots <- create_and_plot_data('UCSF_TA5_R4C2', 2048, 120, -190, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA5_R4C2', 2048, 120, -190, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA6_R1C4' (same)----
plots <- create_and_plot_data('UCSF_TA6_R1C4', 2048, 50, -80, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA6_R1C4', 2048, 50, -80, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA6_R1C5' (same)----
plots <- create_and_plot_data('UCSF_TA6_R1C5', 2048, 200, -90, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA6_R1C5', 2048, 200, -90, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA6_R2C1' (same)----
plots <- create_and_plot_data('UCSF_TA6_R2C1', 2048, 0, 100, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA6_R2C1', 2048, 0, 100, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA6_R2C2' (same)----
plots <- create_and_plot_data('UCSF_TA6_R2C2', 2048, 100, 100, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA6_R2C2', 2048, 100, 100, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA6_R2C3' (same)----
plots <- create_and_plot_data('UCSF_TA6_R2C3', 2048, -70, 30, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA6_R2C3', 2048, -70, 30, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'UCSF_TA6_R2C6' (different)----
plots <- create_and_plot_data_separate_FOVs('UCSF_TA6_R2C6','UCSF_TA6_R1C6', 2048, 200, -70, NULL, NULL, NULL, NULL)
apply_changes_and_save_separate_FOVs('UCSF_TA6_R2C6','UCSF_TA6_R1C6', 2048, 200, -70, NULL, NULL, NULL, NULL,
                                     plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA6_R7C6' (same)----
plots <- create_and_plot_data('UCSF_TA6_R7C6', 2048, -100, -100, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA6_R7C6', 2048, -100, -100, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)

# 'UCSF_TA7_R5C5'----
plots <- create_and_plot_data('UCSF_TA7_R5C5', 2048, -130, -20, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune")
apply_changes_and_save('UCSF_TA7_R5C5', 2048, -130, -20, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune",
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA8_R2C6'----
plots <- create_and_plot_data('UCSF_TA8_R2C6', 2048, -130, -100, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune")
apply_changes_and_save('UCSF_TA8_R2C6', 2048, -130, -100, "cell_meta_cluster_2", "Immune",  "cell_meta_cluster_2", "Immune",
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA8_R4C1'----
plots <- create_and_plot_data('UCSF_TA8_R4C1', 2048, -130, -30, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA8_R4C1', 2048, -130, -30, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA8_R4C2'----
plots <- create_and_plot_data('UCSF_TA8_R4C2', 2048, -100, -80, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA8_R4C2', 2048, -100, -80, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA8_R6C3'----
plots <- create_and_plot_data('UCSF_TA8_R6C3', 2048, -100, 80, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA8_R6C3', 2048, -100, 80, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# 'UCSF_TA8_R7C6'----
plots <- create_and_plot_data('UCSF_TA8_R7C6', 2048, -120, 60, NULL, NULL, NULL, NULL)
apply_changes_and_save('UCSF_TA8_R7C6', 2048, -120, 60, NULL, NULL, NULL, NULL,
                       plots$pre_change_plot, plots$post_change_plot)



# Save data tables----


# Define the file path for the new CSV file
file_path_Immune <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/immune_panel/alignment_modified/immune_alignment_modified.csv"

# Save data_Immune as a CSV file
fwrite(data_Immune, file_path_Immune, row.names = FALSE)

# Define the file path for the new CSV file
file_path_Tumor <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/tumor_panel/alignment_modified/tumor_alignment_modified.csv"

# Save data_Immune as a CSV file
fwrite(data_Tumor, file_path_Tumor, row.names = FALSE)


### combine data tables into a single data table ----

### Load Aligned data ----
file_path_Immune <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/immune_panel/alignment_modified/immune_alignment_modified.csv"
file_path_Tumor <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/tumor_panel/alignment_modified/tumor_alignment_modified.csv"
file_path_metadata <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/metadata/20240529_metadata_complete.csv"
data_Immune <- fread(file_path_Immune)
data_Tumor <- fread(file_path_Tumor)
metadata <- fread(file_path_metadata)
### Remove unnamed columns  ----
# Find columns that do NOT contain "Unnamed: " in their title
named_columns <- grep("Unnamed: ", names(data_Tumor), value = TRUE, invert = TRUE)

# Subset the data.table to keep only the named columns
data_Tumor <- data_Tumor[, ..named_columns]

# Find columns that do NOT contain "Unnamed: " in their title
named_columns <- grep("Unnamed: ", names(data_Immune), value = TRUE, invert = TRUE)

# Subset the data.table to keep only the named columns
data_Immune <- data_Immune[, ..named_columns]

### Align column names ----
data_Tumor <- data_Tumor %>% 
  rename(
    CD3 = CD3e,
    CD3_nuclear = CD3e_nuclear,
    CD8 = CD8a,
    CD8_nuclear = CD8a_nuclear,
    FoxP3 = FOXP3,
    FOXP3_nuclear = FOXP3_nuclear
  )
### Identify unique columns to each data table-----
# Find columns that exist in data_Immune but not in data_Tumor
unique_to_Immune <- setdiff(names(data_Immune), names(data_Tumor))

# Find columns that exist in data_Tumor but not in data_Immune
unique_to_Tumor <- setdiff(names(data_Tumor), names(data_Immune))

# Get data types of the unique columns in data_Immune
data_types_Immune <- sapply(data_Immune[, ..unique_to_Immune], class)

# Get data types of the unique columns in data_Tumor
data_types_Tumor <- sapply(data_Tumor[, ..unique_to_Tumor], class)

# Print the results
if (length(unique_to_Immune) > 0) {
  print("Columns unique to data_Immune and their data types:")
  print(data_types_Immune)
} else {
  print("There are no columns unique to data_Immune.")
}

if (length(unique_to_Tumor) > 0) {
  print("Columns unique to data_Tumor and their data types:")
  print(data_types_Tumor)
} else {
  print("There are no columns unique to data_Tumor.")
}
### Create missing unique columns to each data table ----
# Filter columns that do not contain "_func" in their title
unique_to_Immune_filtered <- unique_to_Immune[!grepl("_func", unique_to_Immune)]
unique_to_Tumor_filtered <- unique_to_Tumor[!grepl("_func", unique_to_Tumor)]

# Function to add missing columns with default values based on type
add_missing_columns <- function(to, from, data_types) {
  for (col in names(data_types)) {
    if (!col %in% names(to)) {
      # Determine default value based on data type
      default_value <- ifelse(data_types[[col]] == "numeric", 0, "0")
      # Add column to 'to' data.table
      set(to, j = col, value = rep(default_value, nrow(to)))
    }
  }
}

# Add missing columns to data_Tumor and data_Immune
add_missing_columns(data_Tumor, data_Immune, data_types_Immune[unique_to_Immune_filtered])
add_missing_columns(data_Immune, data_Tumor, data_types_Tumor[unique_to_Tumor_filtered])
### Assign cell_table of origin ----
data_Immune[, table_origin := 'Immune']
data_Tumor[, table_origin := 'Tumor']


### Assign sample_id and alignment status ----
data_Immune[metadata, on = .(fov = fov_immune),
            `:=` (alignment_status = i.alignment_status, sample_id = i.sample_id)]
data_Tumor[metadata, on = .(fov = fov_tumor),
          `:=` (alignment_status = i.alignment_status, sample_id = i.sample_id)]

### Reassign cell_meta clusters ----

  data_Immune_updated <- data_Immune %>%
    mutate(cell_meta_cluster_final = case_when(
      str_detect(cell_meta_cluster, "^Tumor_") ~ "Tumor_cells", # Assign "Tumor_cells" to strings starting with "Tumor_"
      cell_meta_cluster == "Astrocytes" ~ "Tumor_cells", # Assign "Tumor_cells" to "Astrocytes"
      cell_meta_cluster == "Cancer_cells" ~ "Tumor_cells", # Assign "Tumor_cells" to "Cancer_cells"
      cell_meta_cluster == "Tumor" ~ "Tumor_cells", # Assign "Tumor_cells" to "Tumor"
      cell_meta_cluster == "Bcells_Unknown" ~ "Bcells", # Assign "Bcells" to "Bcells_Unknown"
      cell_meta_cluster == "Oligodendrocytes" ~ "Tumor_cells", # Assign "Tumor_cells" to "Oligodendrocytes"
      cell_meta_cluster == "CD163_Myeloid" ~ "Macrophage_CD163", # Assign "Macrophage_CD163" to "CD163_Myeloid"
      cell_meta_cluster == "Mono_mac_CD163" ~ "Macrophage_CD163", # Assign "Macrophage_CD163" to "Mono_mac_CD163"
      cell_meta_cluster == "Mono_Mac_CD14" ~ "Myeloid_CD14", # Assign "Myeloid_CD14" to "Mono_Mac_CD14"
      cell_meta_cluster == "CD4_Tcells" ~ "Tcell_CD4", # Assign "Tcell_CD4" to "CD4_Tcells"
      cell_meta_cluster == "CD4_cells" ~ "Tcell_CD4", # Assign "Tcell_CD4" to "CD4_cells"
      cell_meta_cluster == "CD8_Tcells" ~ "Tcell_CD8", # Assign "Tcell_CD8" to "CD8_Tcells"
      cell_meta_cluster == "FoxP3_Tcells" ~ "Tcell_FoxP3", # Assign "Tcell_FoxP3" to "FoxP3_Tcells"
      cell_meta_cluster == "FoxP3_Tcell" ~ "Tcell_FoxP3", # Assign "Tcell_FoxP3" to "FoxP3_Tcell"
      cell_meta_cluster == "Treg" ~ "Tcell_FoxP3", # Assign "Tcell_FoxP3" to "Treg"
      cell_meta_cluster == "CD14_Myeloid" ~ "Myeloid_CD14", # Assign "Myeloid_CD14" to "CD14_Myeloid"
      cell_meta_cluster == "CD11b_Myeloid" ~ "Myeloid_CD11b", # Assign "Myeloid_CD11b" to "CD11b_Myeloid"
      cell_meta_cluster == "CD123_DC" ~ "DC_CD123", # Assign "DC_CD123" to "CD123_DC"
      cell_meta_cluster == "CD141_Myeloid" ~ "Myeloid_CD141", # Assign "Myeloid_CD141" to "CD141_Myeloid"
      cell_meta_cluster == "CD206_Macrophages" ~ "Macrophage_CD206", # Assign "Macrophage_CD206" to "CD206_Macrophages"
      cell_meta_cluster == "CD208_DC" ~ "DC_CD208", # Correctly assign "DC_CD208" to "CD208_DC"
      cell_meta_cluster == "CD209_DC_Mac" ~ "DC_Mac_CD209", # Correctly assign "DC_Mac_CD209" to "CD209_DC_Mac"
      cell_meta_cluster == "CD68_Myeloid" ~ "Macrophage_CD68", # Assign "Macrophage_CD68" to "CD68_Myeloid"
      cell_meta_cluster == 'CD20_TMEM119_Unknown' ~ "Microglia", # Assign "Microglia" to 'CD20_TMEM119_Unknown'
      cell_meta_cluster == 'Mast_Cells' ~ "Mast_cells", # Assign "Microglia" to 'CD20_TMEM119_Unknown'
      cell_meta_cluster == "CD45_Only" ~ "Immune_unassigned", # Assign "Immune_unassigned" to "CD45_Only"
      cell_meta_cluster == "CD14_CD163_Myeloid" ~ "Myeloid_CD14_CD163",
      cell_meta_cluster == "CD11b_HLADR_Myeloid" ~ "Myeloid_CD11b_HLADR+",
      cell_meta_cluster == "CD68_CD163_Myeloid" ~ "Macrophage_CD68_CD163",
      cell_meta_cluster == "CD20_TMEM119_CD163_Unknown" ~ "Microglia_CD163",
      cell_meta_cluster == "Immune_other" ~ "Immune_unassigned",
      cell_meta_cluster == "Tcell_Unknown" ~ "Immune_unassigned", # Assign "Immune_unassigned" to "Tcell_Unknown"
      cell_meta_cluster == "Astrocytes_Unknown" ~ "Unassigned", # Assign "Unassigned" to "Astrocytes_Unknown"
      cell_meta_cluster == "Unknown" ~ "Unassigned", # Assign "Unassigned" to "Unknown"
      cell_meta_cluster == "Undetermined" ~ "Unassigned", # Assign "Unassigned" to "Undetermined"
      
      TRUE ~ cell_meta_cluster # Keep other values as they are
    ))

data_Tumor_updated <- data_Tumor %>%
  mutate(cell_meta_cluster_final = case_when(
    str_detect(cell_meta_cluster, "^Tumor_") ~ "Tumor_cells", # Assign "Tumor_cells" to strings starting with "Tumor_"
    cell_meta_cluster == "Astrocytes" ~ "Tumor_cells", # Assign "Tumor_cells" to "Astrocytes"
    cell_meta_cluster == "Cancer_cells" ~ "Tumor_cells", # Assign "Tumor_cells" to "Cancer_cells"
    cell_meta_cluster == "Tumor" ~ "Tumor_cells", # Assign "Tumor_cells" to "Tumor"
    cell_meta_cluster == "Bcells_Unknown" ~ "Bcells", # Assign "Bcells" to "Bcells_Unknown"
    cell_meta_cluster == "Oligodendrocytes" ~ "Tumor_cells", # Assign "Tumor_cells" to "Oligodendrocytes"
    cell_meta_cluster == "CD163_Myeloid" ~ "Macrophage_CD163", # Assign "Macrophage_CD163" to "CD163_Myeloid"
    cell_meta_cluster == "Mono_mac_CD163" ~ "Macrophage_CD163", # Assign "Macrophage_CD163" to "Mono_mac_CD163"
    cell_meta_cluster == "Mono_Mac_CD14" ~ "Myeloid_CD14", # Assign "Myeloid_CD14" to "Mono_Mac_CD14"
    cell_meta_cluster == "CD4_Tcells" ~ "Tcell_CD4", # Assign "Tcell_CD4" to "CD4_Tcells"
    cell_meta_cluster == "CD4_cells" ~ "Tcell_CD4", # Assign "Tcell_CD4" to "CD4_cells"
    cell_meta_cluster == "CD8_Tcells" ~ "Tcell_CD8", # Assign "Tcell_CD8" to "CD8_Tcells"
    cell_meta_cluster == "FoxP3_Tcells" ~ "Tcell_FoxP3", # Assign "Tcell_FoxP3" to "FoxP3_Tcells"
    cell_meta_cluster == "FoxP3_Tcell" ~ "Tcell_FoxP3", # Assign "Tcell_FoxP3" to "FoxP3_Tcell"
    cell_meta_cluster == "Treg" ~ "Tcell_FoxP3", # Assign "Tcell_FoxP3" to "Treg"
    cell_meta_cluster == "CD14_Myeloid" ~ "Myeloid_CD14", # Assign "Myeloid_CD14" to "CD14_Myeloid"
    cell_meta_cluster == "CD11b_Myeloid" ~ "Myeloid_CD11b", # Assign "Myeloid_CD11b" to "CD11b_Myeloid"
    cell_meta_cluster == "CD123_DC" ~ "DC_CD123", # Assign "DC_CD123" to "CD123_DC"
    cell_meta_cluster == "CD141_Myeloid" ~ "Myeloid_CD141", # Assign "Myeloid_CD141" to "CD141_Myeloid"
    cell_meta_cluster == "CD206_Macrophages" ~ "Macrophage_CD206", # Assign "Macrophage_CD206" to "CD206_Macrophages"
    cell_meta_cluster == "CD208_DC" ~ "DC_CD208", # Correctly assign "DC_CD208" to "CD208_DC"
    cell_meta_cluster == "CD209_DC_Mac" ~ "DC_Mac_CD209", # Correctly assign "DC_Mac_CD209" to "CD209_DC_Mac"
    cell_meta_cluster == "CD68_Myeloid" ~ "Macrophage_CD68", # Assign "Macrophage_CD68" to "CD68_Myeloid"
    cell_meta_cluster == 'CD20_TMEM119_Unknown' ~ "Microglia", # Assign "Microglia" to 'CD20_TMEM119_Unknown'
    cell_meta_cluster == 'Mast_Cells' ~ "Mast_cells", # Assign "Microglia" to 'CD20_TMEM119_Unknown'
    cell_meta_cluster == "CD45_Only" ~ "Immune_unassigned", # Assign "Immune_unassigned" to "CD45_Only"
    cell_meta_cluster == "Tcell_Unknown" ~ "Immune_unassigned", # Assign "Immune_unassigned" to "Tcell_Unknown"
    cell_meta_cluster == "Astrocytes_Unknown" ~ "Unassigned", # Assign "Unassigned" to "Astrocytes_Unknown"
    cell_meta_cluster == "Unknown" ~ "Unassigned", # Assign "Unassigned" to "Unknown"
    cell_meta_cluster == "Undetermined" ~ "Tumor_cells", # Assign "Unassigned" to "Undetermined"
    
    TRUE ~ cell_meta_cluster # Keep other values as they are
  ))

### Choose selected cells from each data_table ----
data_Immune_updated_filtered <- data_Immune_updated %>%
  filter(!(cell_meta_cluster_final %in% c("Tumor_cells", "Unassigned"))) %>%
  filter(!(fov == 'CHOP_1101_R1C3')) %>% 
  filter(alignment_status == 1) %>% 
  dplyr::select(-cell_meta_cluster, -cell_meta_cluster_2, -cell_meta_cluster_3, -cell_meta_cluster_4, -cell_meta_cluster_5, -cell_meta_cluster_all)

# UCLA_R99C3 seemed to be missing from cell_meta_cluster_3 so used cell_meta_cluster_final with modified Undetermined cells
data_Tumor_updated_filtered <- data_Tumor_updated %>%
  filter(cell_meta_cluster_final == "Tumor_cells") %>%
  filter(alignment_status == 1) %>% 
  dplyr::select(-cell_meta_cluster, -cell_meta_cluster_2, -cell_meta_cluster_3, -cell_meta_cluster_4, -cell_meta_cluster_5, -cell_meta_cluster_all) 


sum(data_Tumor_updated$cell_meta_cluster_final == "Tumor_cells")
sum(data_Tumor_updated$cell_meta_cluster_final == "Unnassigned")
sum(data_Tumor_updated$cell_meta_cluster_3 == "Tumor")
# Unique sample_id values in data_Immune_updated_filtered not in data_Tumor_updated_filtered
unique_to_immune <- setdiff(unique(data_Immune_updated_filtered$sample_id), unique(data_Tumor_updated_filtered$sample_id))

# Unique sample_id values in data_Tumor_updated_filtered not in data_Immune_updated_filtered
unique_to_tumor <- setdiff(unique(data_Tumor_updated_filtered$sample_id), unique(data_Immune_updated_filtered$sample_id))

# Print the results
unique_to_immune
unique_to_tumor
### Combine both data tables  ----
# Combine the tables, filling missing columns with NA
combined_data <- rbindlist(list(data_Immune_updated_filtered, data_Tumor_updated_filtered), fill = TRUE)

# convert func columns to binary
func_comb_data <- combined_data %>%
  mutate(across(ends_with("_func"), ~ as.integer(!is.na(.) & . != ""), .names = "{.col}")) 

### Check binary data matches final data ----
# Step 1: Count NA and "" in original data
original_counts <- sapply(names(combined_data), function(col_name) {
  if(grepl("_func$", col_name)) {
    col <- combined_data[[col_name]]
    sum(is.na(col) | col == "", na.rm = TRUE)
  } else {
    NA # Indicate non-applicable for columns not ending with "_func"
  }
}, USE.NAMES = TRUE)

# Step 2: Count 0s in converted data
converted_counts <- sapply(names(func_comb_data), function(col_name) {
  if(grepl("_func$", col_name)) {
    col <- func_comb_data[[col_name]]
    sum(col == 0, na.rm = TRUE)
  } else {
    NA # Indicate non-applicable for columns not ending with "_func"
  }
}, USE.NAMES = TRUE)

# Step 3: Compare the counts
comparison <- data.frame(
  Original_Counts = original_counts,
  Converted_Counts = converted_counts
)

# Display the comparison for review
comparison <- na.omit(comparison) # Removing NAs to focus on relevant "_func" columns
print(comparison)

colnames(func_comb_data)
### save combined dataframe ----
# Define the file path for the new CSV file
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"

# Save data_Immune as a CSV file
fwrite(func_comb_data, file_path_combined, row.names = FALSE)


### Load combined data ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"
combined_data <- fread(file_path_combined)

### update additional meta clusters ----
combined_data_updated <- combined_data %>%
  mutate(cell_meta_cluster_ML = case_when(
    cell_meta_cluster_final == "Macrophage_CD163" ~ "Myeloid", 
    cell_meta_cluster_final == "Myeloid_CD14" ~ "Myeloid", 
    cell_meta_cluster_final == "Myeloid_CD14_CD163" ~ "Myeloid",
    cell_meta_cluster_final == "Neutrophils" ~ "Myeloid", 
    cell_meta_cluster_final == "Myeloid_CD11b" ~ "Myeloid", 
    cell_meta_cluster_final == "Myeloid_CD11b_HLADR+" ~ "Myeloid",
    cell_meta_cluster_final == "Macrophage_CD68" ~ "Myeloid", 
    cell_meta_cluster_final == "Macrophage_CD68_CD163" ~ "Myeloid",
    cell_meta_cluster_final == "Microglia"  ~ "Myeloid", 
    cell_meta_cluster_final == "Microglia_CD163"  ~ "Myeloid",
    cell_meta_cluster_final == "DC_CD123" ~ "Myeloid", 
    cell_meta_cluster_final == "Endothelial_cells" ~ "Endothelial", 
    cell_meta_cluster_final == "DC_Mac_CD209" ~ "Endothelial", 
    cell_meta_cluster_final == "Neurons" ~ "Neurons",
    cell_meta_cluster_final == "Tcell_CD8" ~ "Lymphoid", 
    cell_meta_cluster_final == "Myeloid_CD141" ~ "Myeloid",
    cell_meta_cluster_final == "APC" ~ "Myeloid",
    cell_meta_cluster_final == "Macrophage_CD206" ~ "Myeloid",
    cell_meta_cluster_final == "Tcell_CD4" ~ "Lymphoid", 
    cell_meta_cluster_final == "Tcell_FoxP3" ~ "Lymphoid",  
    cell_meta_cluster_final == "DC_CD208" ~ "Myeloid", 
    cell_meta_cluster_final == "DC_Mac_CD209" ~ "Myeloid", 
    cell_meta_cluster_final == "Bcells" ~ "Lymphoid", 
    cell_meta_cluster_final == "Mast_cells" ~ "Myeloid", 
    
    TRUE ~ cell_meta_cluster_final # Keep other values as they are
  ))

combined_data_updated <- combined_data_updated %>%
  mutate(cell_meta_cluster_IT = case_when(
    cell_meta_cluster_ML == "Myeloid" ~ "Immune",
    cell_meta_cluster_ML == "Lymphoid" ~ "Immune",
    
    TRUE ~ cell_meta_cluster_ML # Keep other values as they are
  ))


### save new metaclusters ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"

# Save data_Immune as a CSV file
fwrite(combined_data_updated, file_path_combined, row.names = FALSE)

### Load combined data ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"
combined_data <- fread(file_path_combined)
### update cell labels ----
combined_data <- combined_data %>%
  group_by(sample_id) %>%
  mutate(combined_label = row_number()) %>%
  ungroup() 
### save new labels ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"

# Save data_Immune as a CSV file
fwrite(combined_data, file_path_combined, row.names = FALSE)

### Load combined data ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"
combined_data <- fread(file_path_combined)


### Convert charcter to numeric with NA being replaced with "0" ----
library(data.table)

# Assuming combined_data is your data.table
nuclear_columns <- names(combined_data)[grepl("_nuclear$", names(combined_data))]

for (col_name in nuclear_columns) {
  original_char_col <- as.character(combined_data[[col_name]]) # Convert to character first
  # Replace NA character representations and empty/whitespace-only strings with "0" before conversion
  original_char_col[original_char_col == "" | grepl("^\\s*$", original_char_col)] <- "0"
  numeric_col <- as.numeric(original_char_col) # Then attempt to convert to numeric
  
  # After conversion, all values should be numeric, so additional checks for NAs should not find any.
  # However, we include this step for completeness and future-proofing.
  na_indices <- which(is.na(numeric_col))
  
  if (length(na_indices) > 0) {
    unique_na_values <- unique(original_char_col[na_indices])
    cat("Column:", col_name, "\n")
    cat("Unexpected non-numeric values converted to NA:", toString(unique_na_values), "\n")
    cat("Total unexpected NAs introduced by non-numeric values:", length(na_indices), "\n\n")
  }
  
  # Assign the converted column back, any NA values at this point would be unexpected
  combined_data[[col_name]] <- numeric_col
}

### save new alignment of columns ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"

# Save data_Immune as a CSV file
fwrite(combined_data, file_path_combined, row.names = FALSE)
### Load combined data ----
file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel.csv"
combined_data <- fread(file_path_combined)

# ### incorporate metadata ----
# metadata_immune <- fread("/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/metadata/20240529_metadata_complete.csv")
# metadata_immune <- metadata_immune %>% 
#   dplyr::select(-c(X, V1, alignment_status, combined_label))
# 
# # Perform the merge
# merged_data <- merge(combined_data, metadata_immune, by = "sample_id", all.x = TRUE)
# ### Save combined data with metadata ----
# file_path_combined <- "/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/cell_tables/combined_panel/20240922_combined_panel_with_metadata.csv"
# 
# # Save data_Immune as a CSV file
# fwrite(merged_data, file_path_combined, row.names = FALSE)





