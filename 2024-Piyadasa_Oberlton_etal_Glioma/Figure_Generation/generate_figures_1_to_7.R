stop("Stopping the script execution.")

# 20240529
# Authors - Hadeesha Piyadasa, Benjamin Oberlton

# Enter your main directory here after downloading all the required files from https://bruce.parkerici.org/pages/raw-data-access.html
base_dir <- "your_path/Files"
output_dir <- "your_path/Output"


###### Figure generation notebook #######
###### Load required libraries #################################################
library(ComplexHeatmap)
library(circlize)
library(svglite)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)
library(stringr)
library(ggplot2)
library(gridExtra)
library(arrow)
library(scales)
library(reshape2)
library(forcats)
library(patchwork)
library(rlang)
library(FactoMineR)
library(factoextra)
library(MASS)
library(corrplot)
library(caret)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(webr)
library(ggforce)
library(plotly)
library(tidyverse)
library(ggdendro)
library(htmlwidgets)
library(webshot)
library(ggalluvial)
library(progress)
library(tictoc)
library(extrafont)
library(igraph)
library(ggraph)

###### Functions #####################################################

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

generate_marker_pair_annotation <- function(column_names) {
  # Get all unique combinations of 2 columns
  combinations <- combn(column_names, 2)
  
  # Create a data frame for the combinations
  comb_df <- data.frame(
    Combination = apply(combinations, 2, function(x) paste(x, collapse = "_")),
    stringsAsFactors = FALSE
  )
  
  # Initialize a matrix for the heatmap
  heatmap_matrix <- matrix(NA, nrow = ncol(combinations), ncol = length(column_names))
  colnames(heatmap_matrix) <- column_names
  
  # Fill the heatmap matrix with combination data
  for (i in 1:ncol(combinations)) {
    comb <- combinations[, i]
    heatmap_matrix[i, comb[1]] <- 1
    heatmap_matrix[i, comb[2]] <- 1
  }
  
  # Convert the matrix to a long format data frame
  heatmap_df <- as.data.frame(heatmap_matrix)
  heatmap_df$Row <- 1:nrow(heatmap_df)
  heatmap_df_melt <- melt(heatmap_df, id.vars = "Row", variable.name = "Column", value.name = "Value")
  
  
  # Invert the order of rows
  heatmap_df_melt$Row <- factor(heatmap_df_melt$Row, levels = rev(unique(heatmap_df_melt$Row)))
  
  print(heatmap_df_melt)
  
  # Plot the heatmap-like matrix with modifications
  p <- ggplot(heatmap_df_melt, aes(x = Row, y = Column, fill = as.factor(Value))) +
    geom_tile(color = "black", size = 0.1) +
    scale_fill_manual(values = c("1" = "#85B7CE", "NA" = "white"), na.value = "white", guide = "none") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1,family = "Helvetica", color = "black",margin = margin(t = -2)),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank()) +
    scale_y_discrete(position = "left") +
    coord_flip()
  
  print(p)
  # Export the plot as PNG with 1200 dpi
  ggsave("/Users/piyadasa/Library/CloudStorage/GoogleDrive-piyadasa@stanford.edu/Shared drives/BRUCE_data/resource_paper_drafts/Figures/Figure_2/final/heatmap_matrix_antigen_combinations.png", plot = p, width = 25, height = 65, units = "mm", dpi = 1200)
}

generate_abundance_plots <- function(data_file, metadata_file, column_of_interest, desired_order, nrow_wrap = 1, name, height, width) {
  # Reading the data
  top_seven_tumor_antigen_df <- (data_file)
  
  # Summing marker_3_prop to marker_7_prop and creating new columns
  top_seven_tumor_antigen_df$marker_1_to_2_prop <- rowSums(top_seven_tumor_antigen_df[, c("marker_1_prop", "marker_2_prop")])
  top_seven_tumor_antigen_df$marker_1_to_7_prop <- rowSums(top_seven_tumor_antigen_df[, c("marker_1_prop", "marker_2_prop","marker_3_prop", "marker_4_prop", "marker_5_prop", "marker_6_prop", "marker_7_prop")])
  
  # Filtering out rows where column_of_interest is NA
  filtered_df <- top_seven_tumor_antigen_df[!is.na(top_seven_tumor_antigen_df[[column_of_interest]]), ]
  
  # Melting the dataframe to long format for plotting
  melted_df <- melt(filtered_df, id.vars = c("sample_id","patient_id", column_of_interest), measure.vars = c("marker_1_prop", "marker_1_to_2_prop", "marker_1_to_7_prop"))
  # Recode variable names with symbols, using \n for new lines
  melted_df <- melted_df %>%
    mutate(variable = recode(variable,
                             "marker_1_prop" = "1",  # Single antigen
                             "marker_1_to_2_prop" = "2",  # Two antigens
                             "marker_1_to_7_prop" = "7"))  # Seven antigens
  melted_df <- melted_df %>%
    mutate(across(where(is.numeric), ~ . * 100))
  
  print(melted_df)
  # Creating individual plots with lines connecting points for each unique value in column_of_interest
  unique_values <- unique(melted_df[[column_of_interest]])
  plot_list <- list()
  
  # Grouping by patient_id, recurrence, and variable, then summarizing the value column using the mean
  summarized_df <- melted_df %>%
    group_by(patient_id, recurrence = !!sym(column_of_interest), variable) %>%
    summarize(value = mean(value, na.rm = TRUE)) %>%
    ungroup()
  
  # Creating individual plots with lines connecting points for each unique value in column_of_interest
  unique_values <- unique(summarized_df$recurrence)
  plot_list <- list()
  
  for (value in unique_values) {
    df <- summarized_df[summarized_df$recurrence == value,]
    
    medians <- df %>% group_by(variable) %>% summarize(median_value = median(value))
    
    # Save df as a global variable
    global_df_name <- paste0("df_", value)
    assign(global_df_name, df, envir = .GlobalEnv)
    
    df <- get(global_df_name, envir = .GlobalEnv)
    
    # Adding jitter to the point positions
    set.seed(123)  # For reproducibility
    df_dodged <- df %>%
      mutate(dmeasure = as.numeric(factor(variable)) +
               ifelse(variable == "1", 1 + runif(n(), -0.4, 0.4),
                      ifelse(variable == "2", 2 + runif(n(), -0.4, 0.4),
                             4 + runif(n(), -0.4, 0.4))))
    
    medians_dodged <- medians %>%
      mutate(dmeasure = as.numeric(factor(variable)) +
               ifelse(variable == "1", 1,
                      ifelse(variable == "2", 2,
                             4)))
    
    # Adjust plotting code
    plot <- ggplot() +
      geom_col(data = medians_dodged, aes(x = dmeasure, y = 100, fill = "Antigen Negative"), 
               color = "black", alpha = 0.3, width = 1.5, size = 0.1) +
      geom_col(data = medians_dodged, aes(x = dmeasure, y = median_value, fill = "Antigen Positive"), 
               alpha = 0.7, color = "black", width = 1.5, size = 0.1) +
      geom_point(data = df_dodged, aes(x = dmeasure, y = value), size = 0.1, color = "black") +  
      geom_line(data = df_dodged, aes(group = patient_id, x = dmeasure, y = value), color = "gray", alpha = 0.5, size = 0.2) +
      geom_line(data = medians_dodged, aes(x = dmeasure, y = median_value), color = "#93003a", size = 0.5) +
      scale_x_continuous(breaks = c(2, 4, 7), labels = c("1", "2", "7")) +  # Correctly label and display the x-axis
      ylim(0, 100) +
      labs(x = "", y = "% of tumor cells") +
      scale_fill_manual(name = "Antigen Status", 
                        values = c("Antigen Positive" = "#85B7CE", "Antigen Negative" = "#cccccc")) + 
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 1, color = "black", lineheight = 0.5, family = "Helvetica"),
        axis.text.y = element_text(size = 5, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size = 6, color = "black", family = "Helvetica"),
        plot.title = element_text(size = 5.5, hjust = 0.5, color = "black", family = "Helvetica"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(color = "black", size = 0.25)
      )
    
    plot_list[[value]] <- plot
  }
  
  final_plot <- wrap_plots(plot_list[desired_order], nrow = nrow_wrap)
  
  # Displaying the final plot
  print(final_plot)
  ggsave(name, plot = final_plot, dpi = 1200, width = width, height = height, units = "mm", device = "png")
}

generate_cumulative_line_plot <- function(data, column_of_interest, col_order, col_colors = NULL, output_filename = NULL, width = 60, height = 100) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)  # For custom transformation
  library(grid)    # For unit function
  
  # Filter out rows with NA in the column of interest
  data <- data %>% filter(!is.na(.data[[column_of_interest]]))
  
  # Convert column_of_interest to a factor with levels specified by col_order
  data[[column_of_interest]] <- factor(data[[column_of_interest]], levels = col_order)
  
  # Filter out rows where column_of_interest is NA after factor conversion
  data <- data %>% filter(!is.na(.data[[column_of_interest]]))
  
  # Calculate cumulative sums for markers
  data <- data %>%
    mutate(
      marker_1_prop = marker_1_prop,
      marker_1_to_2_prop = marker_1_prop + marker_2_prop,
      marker_1_to_3_prop = marker_1_to_2_prop + marker_3_prop,
      marker_1_to_4_prop = marker_1_to_3_prop + marker_4_prop,
      marker_1_to_5_prop = marker_1_to_4_prop + marker_5_prop,
      marker_1_to_6_prop = marker_1_to_5_prop + marker_6_prop,
      marker_1_to_7_prop = marker_1_to_6_prop + marker_7_prop
    )
  
  # Pivot data to long format
  melted_df <- data %>%
    select(all_of(column_of_interest), matches("^marker_.*_prop$")) %>%
    pivot_longer(
      cols = matches("^marker_.*_prop$"),
      names_to = "marker_range",
      values_to = "value"
    )
  
  # Calculate mean, SEM, and n for each group and marker_range
  summary_df <- melted_df %>%
    group_by(.data[[column_of_interest]], marker_range) %>%
    summarize(
      n = sum(!is.na(value)),
      mean_value = mean(value, na.rm = TRUE),
      sem_value = sd(value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  
  # Map marker_range to numeric labels starting from 1
  marker_labels <- c(
    "marker_1_prop" = 1,
    "marker_1_to_2_prop" = 2,
    "marker_1_to_3_prop" = 3,
    "marker_1_to_4_prop" = 4,
    "marker_1_to_5_prop" = 5,
    "marker_1_to_6_prop" = 6,
    "marker_1_to_7_prop" = 7
  )
  summary_df$marker_numeric <- marker_labels[summary_df$marker_range]
  
  # Print out n for SEM values
  print(summary_df[, c(column_of_interest, "marker_range", "n")])
  
  # Custom formatter function for y-axis labels
  percent_no_sign <- function(x) {
    round(x * 100)
  }
  
  ### **Compute Difference Between marker_1_to_7_prop and marker_1_prop**
  
  # Filter for marker_1_prop and marker_1_to_7_prop
  diff_df <- summary_df %>%
    filter(marker_range %in% c("marker_1_prop", "marker_1_to_7_prop")) %>%
    select(.data[[column_of_interest]], marker_range, mean_value) %>%
    pivot_wider(names_from = marker_range, values_from = mean_value) %>%
    mutate(
      difference = marker_1_to_7_prop - marker_1_prop
    ) %>%
    select(.data[[column_of_interest]], difference)
  
  # Print out the difference data frame
  print(diff_df)
  
  # Assign the difference data frame as a global variable
  assign("diff_df", diff_df, envir = .GlobalEnv)
  
  ### **Define Custom Transformation Functions**
  
  # Transformation function: Compresses 0-0.4 to 0-0.1 (10% of space)
  trans_break_y <- function(y) {
    ifelse(y <= 0.4, y * 0.25, 0.1 + (y - 0.4) * (0.9 / 0.6))
  }
  
  # Inverse transformation: Maps back to original y-values
  inv_trans_break_y <- function(y_prime) {
    ifelse(y_prime <= 0.1, y_prime / 0.25, 0.4 + (y_prime - 0.1) * (0.6 / 0.9))
  }
  
  # Create the custom transformation object
  break_axis_trans <- trans_new(
    name = "break_axis",
    transform = trans_break_y,
    inverse = inv_trans_break_y
  )
  
  ### **Create the Plot with Custom Y-Axis Transformation**
  
  p <- ggplot(summary_df, aes(
    x = marker_numeric,
    y = mean_value,
    color = .data[[column_of_interest]],
    group = .data[[column_of_interest]]
  )) +
    geom_point(size = 0.25) +     # Set point size to 0.25
    geom_line(size = 0.25) +      # Set line size to 0.25
    geom_errorbar(aes(
      ymin = mean_value,
      ymax = mean_value + sem_value
    ), width = 0.2, size = 0.1) +  # Set error bar line size to 0.1
    scale_x_continuous(
      breaks = 1:7,
      labels = 1:7,
      expand = expansion(add = c(1.2, 0))  # Increase space before first tick
    ) +
    scale_y_continuous(
      trans = break_axis_trans,
      limits = c(0, 1.03),  # Set y-axis limits from 0 to 1
      breaks = c(seq(0, 0.4, by = 0.4), seq(0.5, 1, by = 0.1)),
      labels = percent_no_sign
    ) +
    labs(
      x = "",
      y = "",
      color = NULL  # Remove legend title
    ) +
    theme_minimal() +
    theme(
      # Set base text properties
      text = element_text(family = "Helvetica", size = 8, face = "plain", colour = "black"),
      # Ensure axis titles and text inherit from base text
      axis.title.x = element_text(),
      axis.title.y = element_text(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      # Other theme adjustments
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),  # Set axis line size to 0.25
      axis.ticks = element_line(size = 0.1, color = "black"),  # Tick mark thickness
      axis.ticks.length = unit(0.5, "mm"),                     # Tick mark length
      legend.position = "none"  # Remove legend
    )
  
  # Apply custom colors if provided
  if (!is.null(col_colors)) {
    p <- p + scale_color_manual(values = col_colors)
  }
  
  # Save the plot if output_filename is provided
  if (!is.null(output_filename)) {
    ggsave(filename = output_filename, plot = p, width = width, height = height, units = "mm", dpi = 1200)
  }
  
  return(p)
}

generate_violin_plot <- function(data, column_of_interest, col_order, col_colors = NULL, output_filename = NULL, width = 60, height = 100) {
  library(dplyr)
  library(ggplot2)
  library(grid)    # For unit function
  
  # Filter out rows with NA in the column of interest or true_remaining_prop
  data <- data %>%
    filter(
      !is.na(.data[[column_of_interest]]),
      !is.na(true_remaining_prop)
    )
  
  # Convert column_of_interest to a factor with levels specified by col_order
  data[[column_of_interest]] <- factor(data[[column_of_interest]], levels = col_order)
  
  # Filter out rows where column_of_interest is NA after factor conversion
  data <- data %>% filter(!is.na(.data[[column_of_interest]]))
  
  data <- data %>%
    mutate(true_remaining_prop = true_remaining_prop * 100)
  
  # Print data counts
  print(table(data[[column_of_interest]]))
  
  # Create the violin plot with individual data points
  p <- ggplot(data, aes(
    x = .data[[column_of_interest]],
    y = true_remaining_prop,
    color = .data[[column_of_interest]]
  )) +
    geom_violin(trim = FALSE, fill = NA, size = 0.25) +  # Violin plot outline
    geom_point(position = position_jitter(width = 0.2), size = 0.25) +  # Individual data points
    scale_y_continuous(
      limits = c(0, 102),
      expand = c(0, 0)
    ) +
    labs(
      x = "",
      y = "",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      # Set base text properties
      text = element_text(family = "Helvetica", size = 8, face = "plain", colour = "black"),
      # Ensure axis titles and text inherit from base text
      axis.title.x = element_text(),
      axis.title.y = element_text(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      # Other theme adjustments
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.ticks.y = element_line(size = 0.1, color = "black"),
      axis.ticks.length = unit(0.5, "mm"),
      legend.position = "none"
    )
  
  # Apply custom colors if provided
  if (!is.null(col_colors)) {
    p <- p + scale_color_manual(values = col_colors)
  }
  
  # Save the plot if output_filename is provided
  if (!is.null(output_filename)) {
    ggsave(filename = output_filename, plot = p, width = width, height = height, units = "mm", dpi = 1200)
  }
  
  return(p)
}

generate_marker_pair_summary_heatmap_all <- function(data_file, metadata_file, column_names, column_of_interest, col_order = NULL, color_gradient = colorRamp2(c(0, 1, 50), c("white", "white", "#b81b4a"))) {
  # Reading the data
  top_seven_tumor_antigen_df <- data_file
  
  # Function to remove '_func' from each element
  remove_func_suffix <- function(x) {
    str_replace_all(x, "_func", "")
  }
  
  # Get all unique combinations of 2 columns
  combinations <- combn(column_names, 2)
  
  # Create a data frame for the combinations
  comb_df <- data.frame(
    marker_pair = apply(combinations, 2, function(x) paste(sort(x), collapse = ",")),
    stringsAsFactors = FALSE
  )
  
  # Standardize the order of combinations in top_combinations
  top_seven_tumor_antigen_df$marker_pair <- apply(top_seven_tumor_antigen_df, 1, function(row) {
    markers <- c(row["marker1"], row["marker2"])
    sorted_combination <- paste(sort(markers), collapse = ",")
    return(sorted_combination)
  })
  
  # Create a consistent marker pair by ordering marker names alphabetically
  top_seven_tumor_antigen_df <- top_seven_tumor_antigen_df %>%
    mutate(
      marker_pair = ifelse(marker_1_name < marker_2_name,
                           paste(marker_1_name, marker_2_name, sep = ","),
                           paste(marker_2_name, marker_1_name, sep = ","))
    )
  
  # Summarize data by marker_pair and column_of_interest
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_df %>%
    group_by(!!sym(column_of_interest), marker_pair) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(!!sym(column_of_interest)) %>%
    mutate(percentage = count / sum(count) * 100)
  
  top_seven_tumor_antigen_summary_df$marker_pair <- remove_func_suffix(top_seven_tumor_antigen_summary_df$marker_pair)
  
  # Filter out NA from column_of_interest
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    filter(!is.na(!!sym(column_of_interest)))
  
  # Ungroup the dataframe
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    ungroup()
  
  # Add empty rows for combinations present in comb_df but not in patient_combinations
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    complete(!!sym(column_of_interest), marker_pair = comb_df$marker_pair, fill = list(percentage = 0))
  
  # Exclude rows not found in comb_df$marker_pair
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    semi_join(comb_df, by = "marker_pair")
  
  # Ensure the combinations are factors with the desired order
  combination_levels <- rev(unique(top_seven_tumor_antigen_summary_df$marker_pair))  # Reverse the order of unique combinations
  top_seven_tumor_antigen_summary_df$marker_pair <- factor(top_seven_tumor_antigen_summary_df$marker_pair, levels = combination_levels)
  
  assign("top_seven_tumor_antigen_summary_df", top_seven_tumor_antigen_summary_df, envir = .GlobalEnv)
  
  # Convert to wide format for ComplexHeatmap
  heatmap_data <- top_seven_tumor_antigen_summary_df %>%
    dplyr::select(!!sym(column_of_interest), marker_pair, percentage) %>%
    spread(key = !!sym(column_of_interest), value = percentage, fill = 0)
  
  
  # Convert to data frame to set row names
  heatmap_data <- heatmap_data %>%
    column_to_rownames(var = "marker_pair")
  
  # Flip the order of the rows
  heatmap_data <- heatmap_data[rev(rownames(heatmap_data)), ]
  heatmap_data <- as.matrix(heatmap_data)
  heatmap_data <- heatmap_data[, col_order]
  print(heatmap_data)
  heatmap <- Heatmap(heatmap_data, 
                     name = "% patients", 
                     col = color_gradient, 
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE, 
                     show_row_names = FALSE, 
                     show_column_names = TRUE, 
                     column_names_rot = 90, 
                     border = FALSE, 
                     column_title = NULL, 
                     rect_gp = gpar(col = "black", lwd = 0.25),  # Thinner grid lines
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if (heatmap_data[i, j] != 0) {
                         grid.text(sprintf("%.0f%%", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 4, fontfamily = "Helvetica", col = "black"))
                       }
                     },
                     column_order = col_order,
                     show_heatmap_legend = FALSE,  # Remove the legend
                     row_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"),
                     column_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
  
  draw(heatmap)
  
  # Export the heatmap to a PNG file with correct dimensions
  png(file.path(output_dir, "coverage_heatmap.png"), width = 31, height = 85, units = "mm", res = 1200)
  draw(heatmap)
  dev.off()
  
}

generate_marker_pair_summary_heatmap <- function(data_file, metadata_file, column_names, column_of_interest, col_order = NULL, color_gradient = colorRamp2(c(0, 1, 50), c("white", "white", "#b81b4a"))) {
  # Reading the data
  top_seven_tumor_antigen_df <- data_file
  
  # Function to remove '_func' from each element
  remove_func_suffix <- function(x) {
    str_replace_all(x, "_func", "")
  }
  
  # Get all unique combinations of 2 columns
  combinations <- combn(column_names, 2)
  
  # Create a data frame for the combinations
  comb_df <- data.frame(
    marker_pair = apply(combinations, 2, function(x) paste(sort(x), collapse = ",")),
    stringsAsFactors = FALSE
  )
  
  # Standardize the order of combinations in top_combinations
  top_seven_tumor_antigen_df$marker_pair <- apply(top_seven_tumor_antigen_df, 1, function(row) {
    markers <- c(row["marker1"], row["marker2"])
    sorted_combination <- paste(sort(markers), collapse = ",")
    return(sorted_combination)
  })
  
  # Create a consistent marker pair by ordering marker names alphabetically
  top_seven_tumor_antigen_df <- top_seven_tumor_antigen_df %>%
    mutate(
      marker_pair = ifelse(marker_1_name < marker_2_name,
                           paste(marker_1_name, marker_2_name, sep = ","),
                           paste(marker_2_name, marker_1_name, sep = ","))
    )
  
  # Summarize data by marker_pair and column_of_interest
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_df %>%
    group_by(!!sym(column_of_interest), marker_pair) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(!!sym(column_of_interest)) %>%
    mutate(percentage = count / sum(count) * 100)
  
  top_seven_tumor_antigen_summary_df$marker_pair <- remove_func_suffix(top_seven_tumor_antigen_summary_df$marker_pair)
  
  # Filter out NA from column_of_interest
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    filter(!is.na(!!sym(column_of_interest)))
  
  # Ungroup the dataframe
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    ungroup()
  
  # Add empty rows for combinations present in comb_df but not in patient_combinations
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    complete(!!sym(column_of_interest), marker_pair = comb_df$marker_pair, fill = list(percentage = 0))
  
  # Exclude rows not found in comb_df$marker_pair
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    semi_join(comb_df, by = "marker_pair")
  
  # Ensure the combinations are factors with the desired order
  combination_levels <- rev(unique(top_seven_tumor_antigen_summary_df$marker_pair))  # Reverse the order of unique combinations
  top_seven_tumor_antigen_summary_df$marker_pair <- factor(top_seven_tumor_antigen_summary_df$marker_pair, levels = combination_levels)
  
  top_seven_tumor_antigen_summary_df <- top_seven_tumor_antigen_summary_df %>%
    filter(
      final_diagnosis_simple %in% c("GBM", "Astrocytoma", "Oligodendroglioma") &
        (
          str_detect(marker_pair, "B7H3") |
            marker_pair %in% c("EGFR,HER2", "GM2_GD2,HER2")
        )
    )
  
  assign("top_seven_tumor_antigen_summary_df", top_seven_tumor_antigen_summary_df, envir = .GlobalEnv)
  
  # Convert to wide format for ComplexHeatmap
  heatmap_data <- top_seven_tumor_antigen_summary_df %>%
    dplyr::select(!!sym(column_of_interest), marker_pair, percentage) %>%
    spread(key = !!sym(column_of_interest), value = percentage, fill = 0)
  
  
  # Convert to data frame to set row names
  heatmap_data <- heatmap_data %>%
    column_to_rownames(var = "marker_pair")
  
  # Flip the order of the rows
  heatmap_data <- heatmap_data[rev(rownames(heatmap_data)), ]
  heatmap_data <- as.matrix(heatmap_data)
  heatmap_data <- heatmap_data[, col_order]
  print(heatmap_data)
  heatmap <- Heatmap(heatmap_data, 
                     name = "% patients", 
                     col = color_gradient, 
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE, 
                     show_row_names = FALSE, 
                     show_column_names = FALSE, 
                     column_names_rot = 90, 
                     border = FALSE, 
                     column_title = NULL, 
                     rect_gp = gpar(col = "black", lwd = 0.25),  # Thinner grid lines
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if (heatmap_data[i, j] != 0) {
                         grid.text(sprintf("%.0f%%", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 5, fontfamily = "Helvetica", col = "black"))
                       }
                     },
                     column_order = col_order,
                     show_heatmap_legend = FALSE,  # Remove the legend
                     row_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"),
                     column_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
  
  draw(heatmap)
  
  # Export the heatmap to a PNG file with correct dimensions
  png(file.path(output_dir, "coabundance_heatmap_output_updated.png"), width = 30, height = 30, units = "mm", res = 1200)
  draw(heatmap)
  dev.off()
  
}

generate_coexpression_heatmap <- function(top_combinations_file, tumor_cells_file, metadata_file, column_names, column_of_interest, col_order = NULL, color_gradient = colorRamp2(c(0, 1, 50), c("white", "white", "#b81b4a"))) {
  # Read the data
  top_combinations <- top_combinations_file
  tumor_cells <- tumor_cells_file
  
  # Function to remove '_func' from each element
  remove_func_suffix <- function(x) {
    str_replace_all(x, "_func", "")
  }
  
  # Apply the function to marker1, marker2, and combination columns
  top_combinations$marker1 <- remove_func_suffix(top_combinations$marker1)
  top_combinations$marker2 <- remove_func_suffix(top_combinations$marker2)
  top_combinations$combination <- remove_func_suffix(top_combinations$combination)
  
  # Standardize the order of combinations in top_combinations
  top_combinations$combination <- apply(top_combinations, 1, function(row) {
    markers <- c(row["marker1"], row["marker2"])
    sorted_combination <- paste(sort(markers), collapse = ",")
    return(sorted_combination)
  })
  
  # Get all unique combinations of 2 columns
  combinations <- combn(column_names, 2)
  
  # Create a data frame for the combinations
  comb_df <- data.frame(
    Combination = apply(combinations, 2, function(x) paste(sort(x), collapse = ",")),
    stringsAsFactors = FALSE
  )
  
  # Aggregate by column_of_interest
  col_of_interest <- ensym(column_of_interest)
  patient_combinations <- top_combinations %>%
    left_join(tumor_cells %>% dplyr::select(patient_id, !!col_of_interest) %>% distinct(), by = "patient_id") %>%
    group_by(!!col_of_interest, combination) %>%
    summarise(patient_count = n(), .groups = 'drop')
  
  total_patients <- patient_combinations %>%
    group_by(!!col_of_interest) %>%
    summarise(total_patients = sum(patient_count), .groups = 'drop')
  
  patient_combinations <- patient_combinations %>%
    left_join(total_patients, by = as_name(col_of_interest)) %>%
    mutate(percentage = (patient_count / total_patients) * 100)
  
  # Filter out NA from column_of_interest
  patient_combinations <- patient_combinations %>%
    filter(!is.na(!!col_of_interest))
  
  # Add empty rows for combinations present in comb_df but not in patient_combinations
  patient_combinations <- patient_combinations %>%
    complete(!!col_of_interest, combination = comb_df$Combination, fill = list(percentage = 0))
  
  # Ensure the combinations are factors with the desired order
  combination_levels <- rev(unique(patient_combinations$combination))  # Reverse the order of unique combinations
  patient_combinations$combination <- factor(patient_combinations$combination, levels = combination_levels)
  
  # Convert to wide format for ComplexHeatmap
  heatmap_data <- patient_combinations %>%
    dplyr::select(!!col_of_interest, combination, percentage) %>%
    spread(key = as_name(col_of_interest), value = percentage, fill = 0)
  
  # Convert to data frame to set row names
  heatmap_data <- heatmap_data %>%
    column_to_rownames(var = "combination")
  rownames(heatmap_data) <- gsub("_", "/", rownames(heatmap_data))
  
  # Flip the order of the rows
  heatmap_data <- heatmap_data[rev(rownames(heatmap_data)), ]
  heatmap_data <- as.matrix(heatmap_data)
  print(heatmap_data)
  heatmap <- Heatmap(heatmap_data, 
                     name = "% patients", 
                     col = color_gradient, 
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE, 
                     show_row_names = FALSE, 
                     show_column_names = TRUE, 
                     column_names_rot = 90, 
                     border = FALSE, 
                     column_title = NULL, 
                     rect_gp = gpar(col = "black", lwd = 0.25),  # Thinner grid lines
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if (heatmap_data[i, j] != 0) {
                         grid.text(sprintf("%.0f%%", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 4, fontfamily = "Helvetica", col = "black"))
                       }
                     },
                     column_order = col_order,
                     show_heatmap_legend = FALSE,  # Remove the legend
                     row_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"),
                     column_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
  
  
  print(heatmap)
  
  # Export the heatmap to a PNG file with correct dimensions
  png(file.path(output_dir, "heatmap_output_co_expression.png"), width = 31, height = 85, units = "mm", res = 1200)
  draw(heatmap)
  dev.off()
}

generate_non_cumulative_bar_plot <- function(data, column_of_interest, col_order) {
  # Filter out rows with NA in the column of interest
  data <- data %>% filter(!is.na(!!sym(column_of_interest)))
  
  # Select only the columns we need
  marker_columns <- c("marker_1_prop", "marker_2_prop", "marker_3_prop", 
                      "marker_4_prop", "marker_5_prop", "marker_6_prop", "marker_7_prop")
  selected_columns <- c(column_of_interest, marker_columns)
  
  # Summarize the data by calculating the median for each unique value in column_of_interest
  summary_df <- data %>%
    group_by(!!sym(column_of_interest)) %>%
    summarize(across(all_of(marker_columns), median, na.rm = TRUE), .groups = 'drop')
  
  # Convert column_of_interest to a factor with levels specified by col_order
  summary_df[[column_of_interest]] <- factor(summary_df[[column_of_interest]], levels = col_order)
  
  # Melt the summarized dataframe to long format for plotting
  melted_df <- summary_df %>%
    pivot_longer(cols = -all_of(column_of_interest), names_to = "marker_range", values_to = "value")
  
  # Save the resulting data frame to the global environment
  assign("melted_df", melted_df, envir = .GlobalEnv)
  
  # Creating the bar plot with points for each unique patient_id
  ggplot(melted_df, aes(x = marker_range, y = value, fill = !!sym(column_of_interest))) +
    geom_col(position = "dodge", color = "black") +  # Bar plot using geom_col
    scale_x_discrete(labels = c("\u25A0", "\u25A0\u25A0", "\u25A0\u25A0\u25A0", "\u25A0\u25A0\u25A0\u25A0", 
                                "\u25A0\u25A0\u25A0\u25A0\u25A0", "\u25A0\u25A0\u25A0\u25A0\u25A0\u25A0", "\u25A0\u25A0\u25A0\u25A0\u25A0\u25A0\u25A0")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1), labels = function(x) round(x * 100)) +  # Set y-axis breaks, limits, and custom labels
    labs(x = "Cumulative Tumor Antigens", y = "Tumor Cells (%)", fill = column_of_interest, color = column_of_interest) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5, color = "red", face = "bold", lineheight = 0.5), 
          legend.position = "right")  # Remove legend and make x-axis text red and bold
}

generate_coexp_percentage_heatmap <- function(tumor_cells_file, metadata, column_of_interest, col_order, comparison_type = "greater_than_or_equal_to") {
  # Define the column names
  column_names <- c("B7H3", "EGFR", "HER2", "NG2", "GM2_GD2", "GPC2", "VISTA")
  column_names_func <- paste0(column_names, "_func")
  
  # Calculate the percentage of cells based on comparison_type
  summary_df <- tumor_cells_file %>%
    group_by(patient_id) %>%
    summarise(
      total_cells = n(),
      at_least_1 = sum(rowSums(across(all_of(column_names_func))) >= 1) / total_cells * 100,
      at_least_2 = sum(rowSums(across(all_of(column_names_func))) >= 2) / total_cells * 100,
      at_least_3 = sum(rowSums(across(all_of(column_names_func))) >= 3) / total_cells * 100,
      at_least_4 = sum(rowSums(across(all_of(column_names_func))) >= 4) / total_cells * 100,
      at_least_5 = sum(rowSums(across(all_of(column_names_func))) >= 5) / total_cells * 100,
      at_least_6 = sum(rowSums(across(all_of(column_names_func))) >= 6) / total_cells * 100,
      at_least_7 = sum(rowSums(across(all_of(column_names_func))) >= 7) / total_cells * 100,
      exactly_1 = sum(rowSums(across(all_of(column_names_func))) == 1) / total_cells * 100,
      exactly_2 = sum(rowSums(across(all_of(column_names_func))) == 2) / total_cells * 100,
      exactly_3 = sum(rowSums(across(all_of(column_names_func))) == 3) / total_cells * 100,
      exactly_4 = sum(rowSums(across(all_of(column_names_func))) == 4) / total_cells * 100,
      exactly_5 = sum(rowSums(across(all_of(column_names_func))) == 5) / total_cells * 100,
      exactly_6 = sum(rowSums(across(all_of(column_names_func))) == 6) / total_cells * 100,
      exactly_7 = sum(rowSums(across(all_of(column_names_func))) == 7) / total_cells * 100
    )
  
  # Select the appropriate columns based on comparison_type
  if (comparison_type == "greater_than_or_equal_to") {
    summary_df <- summary_df %>%
      select(patient_id, total_cells, starts_with("at_least"))
  } else if (comparison_type == "equal_to") {
    summary_df <- summary_df %>%
      select(patient_id, total_cells, starts_with("exactly"))
  }
  
  # Remove duplicates based on patient_id in metadata
  metadata <- metadata %>%
    distinct(patient_id, .keep_all = TRUE)
  
  # Merge with metadata to get final_diagnosis_simple
  summary_df <- left_join(summary_df, metadata, by = "patient_id")
  
  # Remove rows where final_diagnosis_simple is NA
  summary_df <- summary_df %>%
    filter(!is.na(final_diagnosis_simple))
  
  # Summarize by final_diagnosis_simple
  summary_df <- summary_df %>%
    group_by(!!sym(column_of_interest)) %>%
    summarise(across(starts_with(if (comparison_type == "greater_than_or_equal_to") "at_least" else "exactly"), mean))
  
  # Convert column_of_interest to a factor with levels specified by col_order
  summary_df[[column_of_interest]] <- factor(summary_df[[column_of_interest]], levels = col_order)
  
  # Melt the dataframe to long format for plotting
  melted_df <- summary_df %>%
    pivot_longer(cols = starts_with(if (comparison_type == "greater_than_or_equal_to") "at_least" else "exactly"), names_to = "marker_count", values_to = "percent_positive") %>%
    mutate(percent_positive = ifelse(percent_positive <= 1, 0, percent_positive))
  
  # Prepare the matrix for the heatmap
  heatmap_matrix <- melted_df %>%
    pivot_wider(names_from = marker_count, values_from = percent_positive) %>%
    column_to_rownames(var = column_of_interest) %>%
    as.matrix()
  
  # Rename column names to 1 to 7
  colnames(heatmap_matrix) <- 1:7
  
  # Define the color function
  col_fun <- colorRamp2(c(0, 1, 1, 100), c("gray", "white", "white", "red"))
  
  # Rename column names to include >= if comparison_type is "greater_than_or_equal_to"
  if (comparison_type == "greater_than_or_equal_to") {
    colnames(heatmap_matrix) <- paste0("≥ ", 1:7)
  } else {
    colnames(heatmap_matrix) <- 1:7
  }
  
  
  # Plotting the heatmap
  Heatmap(heatmap_matrix, name = "% of Tumor Cells", col = col_fun,
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_dend = FALSE, show_column_dend = FALSE,
          column_title = "Number of Co-expressed Antigens", row_title = NULL,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", heatmap_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
          },
          rect_gp = gpar(col = "black", lwd = 1),
          row_order = col_order,
          row_names_side = "left",  # Move row names to the left
          column_names_rot = 0)  # Rotate x-axis labels 90 degrees
}

run_analysis <- function(data, group_var, columns, analysis_type = "pca_lda", use_pca = TRUE, use_corr = FALSE, corr_threshold = 0.9, color_var = NULL) {
  
  # Function to handle NA values
  handle_na <- function(data) {
    data[is.na(data)] <- 0
    return(data)
  }
  
  # Function to scale data
  scale_data <- function(data) {
    return(scale(data))
  }
  
  # Function to plot histograms
  plot_histograms <- function(data, scaled_data) {
    # Extract the first column
    first_column <- data[[1]]
    scaled_first_column <- scaled_data[, 1]
    
    # Create a combined dataframe for plotting
    combined_df <- data.frame(
      value = c(first_column, scaled_first_column),
      type = rep(c("Before Scaling", "After Scaling"), each = length(first_column))
    )
    
    # Plot the histograms
    hist_plot <- ggplot(combined_df, aes(x = value, fill = type)) +
      geom_histogram(binwidth = 1, alpha = 0.5, position = "identity") +
      facet_wrap(~ type, scales = "free") +
      theme_minimal() +
      labs(title = "Histogram of the First Column Before and After Scaling", x = "Value", y = "Frequency")
    
    print(hist_plot)
  }
  
  # Select the columns for analysis
  data_selected <- data %>%
    dplyr::select(all_of(columns))
  
  data_selected <- handle_na(data_selected)
  scaled_data_selected <- as.data.frame(scale_data(data_selected))
  
  # Plot histograms before and after scaling
  plot_histograms(data_selected, scaled_data_selected)
  
  if (is.null(color_var)) {
    color_var <- group_var
  }
  
  if (analysis_type == "pca") {
    pca_result <- PCA(scaled_data_selected, graph = FALSE)
    
    # Export PC1 and PC2 components for the variables and samples
    pca_var_components <- pca_result$var$coord[, 1:2]
    pca_sample_components <- pca_result$ind$coord[, 1:2]
    rownames(pca_sample_components) <- rownames(data)
    
    assign("PCA_var_components", pca_var_components, envir = .GlobalEnv)
    assign("PCA_sample_components", pca_sample_components, envir = .GlobalEnv)
    
    print("PCA_var_components and PCA_sample_components saved to global environment.")
    
    # Plot PCA
    pca_plot <- fviz_pca_ind(pca_result, 
                             geom.ind = "point", 
                             col.ind = data[[color_var]], 
                             palette = "jco",
                             addEllipses = TRUE, ellipse.level = 0.95,
                             ellipse.type ="confidence",
                             legend.title = color_var) +
      labs(title = paste("PCA of Stanford MALDI Data by", color_var))
    print(pca_plot)
  } else if (analysis_type == "pca_lda") {
    if (use_pca) {
      pca_result <- PCA(scaled_data_selected, graph = FALSE)
      
      # Export PC1 and PC2 components for the variables and samples
      PCA_lda_var_components <- pca_result$var$coord[, 1:2]
      PCA_lda_sample_components <- pca_result$ind$coord[, 1:2]
      rownames(PCA_lda_sample_components) <- rownames(data)
      
      assign("PCA_lda_var_components", PCA_lda_var_components, envir = .GlobalEnv)
      assign("PCA_lda_sample_components", PCA_lda_sample_components, envir = .GlobalEnv)
      
      print("PCA_lda_var_components and PCA_sample_components saved to global environment.")
      
      # Plot PCA
      pca_plot <- fviz_pca_ind(pca_result, 
                               geom.ind = "point", 
                               col.ind = data[[color_var]], 
                               palette = "jco",
                               addEllipses = TRUE, ellipse.level = 0.95,
                               ellipse.type ="confidence", 
                               legend.title = color_var) +
        labs(title = paste("PCA of Stanford MALDI Data by", color_var))
      print(pca_plot)
      
      # Use PCA values for LDA
      pca_values <- pca_result$ind$coord
      pca_df <- data.frame(pca_values, group = data[[group_var]])
      names(pca_df)[ncol(pca_df)] <- group_var  # Rename the group column
      lda_result <- lda(as.formula(paste(group_var, "~ .")), data = pca_df)
    } else if (use_corr) {
      # Calculate the correlation matrix
      corr_matrix <- cor(scaled_data_selected, use = "pairwise.complete.obs")
      
      # Plot the correlation matrix
      corrplot(corr_matrix, method = "color", tl.cex = 0.7, number.cex = 0.7)
      
      # Identify highly correlated variables
      high_corr <- findCorrelation(corr_matrix, cutoff = corr_threshold, names = TRUE)
      
      # Print highly correlated variables
      print(high_corr)
      
      # Remove highly correlated variables
      data_reduced <- scaled_data_selected %>%
        dplyr::select(-one_of(high_corr))
      
      data_reduced <- handle_na(data_reduced)
      
      # Add group variable back to the reduced data
      data_reduced <- cbind(data[[group_var]], data_reduced)
      colnames(data_reduced)[1] <- group_var  # Ensure the first column name matches group_var
      
      # Perform LDA on reduced data
      lda_result <- lda(as.formula(paste(group_var, "~ .")), data = data_reduced)
    } else {
      # Add group variable back to the selected data
      scaled_data_selected <- cbind(data[[group_var]], scaled_data_selected)
      colnames(scaled_data_selected)[1] <- group_var  # Ensure the first column name matches group_var
      
      # Perform LDA directly on the selected data
      lda_result <- lda(as.formula(paste(group_var, "~ .")), data = scaled_data_selected)
    }
    
    # Export LD1 and LD2 coefficients for the variables and samples
    lda_values <- predict(lda_result)$x
    lda_coefficients <- lda_result$scaling[, 1:2]
    lda_sample_components <- lda_values[, 1:2]
    rownames(lda_sample_components) <- rownames(data)
    
    assign("LDA_var_coefficients", lda_coefficients, envir = .GlobalEnv)
    assign("LDA_sample_components", lda_sample_components, envir = .GlobalEnv)
    
    # Create a data frame with LDA results
    lda_values <- predict(lda_result)$x
    lda_df <- data.frame(lda_values, group = data[[color_var]])
    names(lda_df)[ncol(lda_df)] <- color_var  # Rename the group column
    
    # Plot LDA results with decision boundaries
    lda_plot <- ggplot(lda_df, aes_string(x = "LD1", y = "LD2", color = color_var, fill = color_var)) +
      geom_point() +
      stat_ellipse(geom = "polygon", alpha = 0.1) +
      scale_color_jco() +  # Use the jco palette
      scale_fill_jco() +   # Ensure fill uses the same palette
      labs(color = "WHO Grade", fill = "WHO Grade", x = "LD1", y = "LD2") +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(color = "black", fill = NA),  # Add black border
        plot.title = element_blank()  # Remove title
      )
    print(lda_plot)
  } else {
    stop("Invalid analysis type. Choose 'pca' or 'pca_lda'.")
  }
}

###### Load files required #####################################################

##### all files listed here can be downloaded through - https://bruce.parkerici.org/pages/raw-data-access.html #####

metadata <- read.csv(file.path(base_dir, "metadata_complete.csv"))
master_feature_table_filtered <- readRDS(file.path(base_dir, "20240810_master_feature_table_na_removed.rds"))
MALDI_data <- read.csv(file.path(base_dir, "20241112_MALDI_BRUCE_annotated.csv"))
classes <- read.csv(file.path(base_dir, "glycan_classes.csv"))
tumor_antigen_segment_df <- read_parquet(file.path(base_dir, "tumor_antigen_segment_df.parquet"))
cell_table_all_merged_thresholded <- read_parquet(file.path(base_dir, "cell_table_all_merged_thresholded.parquet"))
data_file <- read_parquet(file.path(base_dir, "top_seven_df.parquet"))
tumor_count_result <- read_parquet(file.path(base_dir, "tumor_count_result.parquet"))
all_cell_count_tumor_FOV_result <- read_parquet(file.path(base_dir, "all_cell_count_tumor_FOV_result.parquet"))
antigen_intensity <- read_parquet(file.path(base_dir, "melted_cell_table_marker_intensity_median.parquet"))
auc_df_WHO <- read_csv(file.path(base_dir, "auc_data_who_grade.csv"))
importance_df_filtered <- read_csv(file.path(base_dir, "rf_gini_importance_df_who_grade.csv")) %>%
  arrange(desc(Importance))
importance_df_filtered_survival <- read_csv(file.path(base_dir, "rf_gini_importance_df_survival_status.csv")) %>%
  arrange(desc(Importance))
auc_df_GBM <- read_csv(file.path(base_dir, "auc_data_survival_status.csv"))
PCA <- read_csv(file.path(base_dir, "pca_subset.csv"))
fold_change_dt <- arrow::read_parquet(file.path(base_dir, "fold_change_dt.parquet"))
merged_all_top_go_df <- readRDS(file.path(base_dir, "merged_all_top_go_df_fig_5.rds"))
significant_data <- readRDS(file.path(base_dir, "glycan_cell_significant_data_fig_5.rds"))

###### lists ###################################################################

cell_type_immune <- c("Myeloid_CD14", "Immune_unassigned", 
                      "APC", "Tcell_CD8", 
                      "Endothelial_cells", "Tcell_FoxP3", "Tcell_CD4", 
                      "Microglia", "Macrophage_CD68", 
                      "Macrophage_CD206", "Myeloid_CD141", "Myeloid_CD11b_HLADRminus", 
                      "Myeloid_CD11b_HLADRplus", "Neutrophils", 
                      "DC_Mac_CD209", "Bcells", "Mast_cells", 
                      "Myeloid_CD14_CD163","Macrophage_CD68_CD163", 
                      "Microglia_CD163", "Unassigned")



tumor_func_columns_extended <- c(
  "B7H3_func_over_all_tumor_count_prop", "EGFR_func_over_all_tumor_count_prop",
  "GM2_GD2_func_over_all_tumor_count_prop", "GPC2_func_over_all_tumor_count_prop",
  "HER2_func_over_all_tumor_count_prop", "NG2_func_over_all_tumor_count_prop",
  "VISTA_func_over_all_tumor_count_prop"
)

tumor_columns_core <- c("B7H3", 
                        "EGFR",  
                        "GM2_GD2", "GPC2",  
                        "HER2", 
                        "NG2", "VISTA")

# Define the desired order for the X-axis labels for who_diagnosis figures
desired_order_who_diagnosis <- c(
  "4 - GBM",
  "4 - Astrocytoma",
  "3 - Astrocytoma",
  "2 - Astrocytoma",
  "3 - Oligodendroglioma",
  "2 - Oligodendroglioma",
  "3 - PXA",
  "2 - PXA",
  "Pediatric DIPG",
  "Pediatric HGG (other)"
)
###### Prep ###################################################################

master_feature_table_filtered <- combine_data_metadata(master_feature_table_filtered, metadata, merge_by = "sample_id", remove_duplicates = FALSE)
data_file <- combine_data_metadata(data_file, metadata, merge_by = "patient_id")


###### Figure 1 ################################################################
# Panel a Cohort graphic ----

# Filter rows where at least one of the columns MIBI, NS, or MALDI has a value of 1
metadata <- metadata[rowSums(metadata[, c("MIBI", "NS", "MALDI")] == 1) > 0, ]

# Create tables
diagnosis_table <- table(metadata$final_diagnosis_simple, useNA = "ifany")
who_grade_table <- table(metadata$who_grade)
tumor_region_table <- table(metadata$tumor_region)
survival_table <- table(metadata$survival_status)
idh_status_table <- table(metadata$idh_status)
recurrence_table <- table(metadata$recurrence)
immunotherapy_yes_table <- table(metadata$cohort[metadata$immunotherapy == "yes"])
immunotherapy_no_count <- sum(metadata$immunotherapy == "no")
paired_yes_data <- metadata[metadata$paired == "yes" & !is.na(metadata$paired), ]
paired_yes_table <- table(paired_yes_data$final_diagnosis_simple, useNA = "ifany")
immunotherapy_paired_yes_count <- sum(paired_yes_data$immunotherapy == "yes", na.rm = TRUE)

# Calculate metrics
total_count <- sum(diagnosis_table)
filtered_table <- diagnosis_table[!names(diagnosis_table) %in% c("GBM_other", NA)]
sum_excluding <- sum(filtered_table)
difference <- total_count - sum_excluding

# Print the results
cat(
  "Unique patients:", length(unique(metadata$patient_id)), "\n",
  "Unique samples:", length(unique(metadata$sample_id)), "\n"
)

cat(
  "Diagnosis Breakdown:\n",
  "Astrocytoma:", filtered_table["Astrocytoma"], "\n",
  "GBM:", filtered_table["GBM"], "\n",
  "Oligodendroglioma:", filtered_table["Oligodendroglioma"], "\n",
  "Pediatric DIPG:", filtered_table["Pediatric DIPG"], "\n",
  "Pediatric HGG (other):", filtered_table["Pediatric HGG (other)"], "\n",
  "PXA:", filtered_table["PXA"], "\n",
  "Other:", difference, "\n"
)

cat(
  "WHO Grade Breakdown:\n",
  "Grade 2:", who_grade_table["2"], "\n",
  "Grade 3:", who_grade_table["3"], "\n",
  "Grade 4:", who_grade_table["4"], "\n"
)

cat(
  "Tumor Region Breakdown:\n",
  paste(names(tumor_region_table), tumor_region_table, sep = ": ", collapse = "\n"), "\n"
)

cat(
  "Survival Breakdown:\n",
  paste(names(survival_table), survival_table, sep = ": ", collapse = "\n"), "\n"
)

cat(
  "IDH Status Breakdown:\n",
  paste(names(idh_status_table), idh_status_table, sep = ": ", collapse = "\n"), "\n"
)

cat(
  "Recurrence Breakdown:\n",
  paste(names(recurrence_table), recurrence_table, sep = ": ", collapse = "\n"), "\n"
)

# Combine and rename groups
cohort_combined <- c(
  Antibody = sum(immunotherapy_yes_table[c("neoadjuvant", "neoadjuvant_resp", "neoadjuvant_nonresp")], na.rm = TRUE),
  Vaccine = sum(immunotherapy_yes_table[c("neoadjuvant_lys_vaccine", "neoadjuvant_SPORE_vaccine")], na.rm = TRUE),
  Combinatorial = sum(immunotherapy_yes_table[c("neoadjuvant_SPORE_CD27")], na.rm = TRUE)
)

cat(
  "Combined Cohort Breakdown for Immunotherapy (Yes):\n",
  paste(names(cohort_combined), cohort_combined, sep = ": ", collapse = "\n"), "\n"
)

cat(
  "Count of Immunotherapy (No):", immunotherapy_no_count, "\n"
)

cat(
  "Final Diagnosis Breakdown for Longitudinal (Yes):\n",
  paste(names(paired_yes_table), paired_yes_table, sep = ": ", collapse = "\n"), "\n"
)

cat(
  "Final Diagnosis Breakdown for Paired (Yes):\n",
  paste(names(paired_yes_table), paired_yes_table, sep = ": ", collapse = "\n"), "\n",
  "Number of Samples with Immunotherapy (Yes) and Paired (Yes):", immunotherapy_paired_yes_count, "\n"
)


# Filter for NS == 1
metadata_ns <- metadata %>% filter(NS == 1)
unique_sample_id_count_ns <- length(unique(metadata_ns$sample_id))
unique_patient_id_count_ns <- length(unique(metadata_ns$patient_id))

# Filter for MIBI == 1
metadata_mibi <- metadata %>% filter(MIBI == 1)
unique_sample_id_count_mibi <- length(unique(metadata_mibi$sample_id))
unique_patient_id_count_mibi <- length(unique(metadata_mibi$patient_id))

# Filter for MALDI == 1
metadata_maldi <- metadata %>% filter(MALDI == 1)
unique_sample_id_count_maldi <- length(unique(metadata_maldi$sample_id))
unique_patient_id_count_maldi <- length(unique(metadata_maldi$patient_id))

# Print results
cat(
  "Unique counts for NS == 1:\n",
  "Sample ID:", unique_sample_id_count_ns, "\n",
  "Patient ID:", unique_patient_id_count_ns, "\n\n",
  
  "Unique counts for MIBI == 1:\n",
  "Sample ID:", unique_sample_id_count_mibi, "\n",
  "Patient ID:", unique_patient_id_count_mibi, "\n\n",
  
  "Unique counts for MALDI == 1:\n",
  "Sample ID:", unique_sample_id_count_maldi, "\n",
  "Patient ID:", unique_patient_id_count_maldi, "\n"
)


###### Figure 2 ################################################################
# Panel b pie and stacked bar ----


# Assuming your dataframe is called filtered_data
filtered_data <- cell_table_all_merged_thresholded %>%
  filter((cell_meta_cluster_final_broad == "Tumor" & source == "cell_table_tumor") |
           (cell_meta_cluster_final_broad == "Immune" & source == "cell_table_immune") |
           (cell_meta_cluster_final_broad == "Neurons" & source == "cell_table_tumor") |
           (cell_meta_cluster_final_broad == "Endothelial" & source == "cell_table_immune")) %>%
  filter(cell_meta_cluster_final != "Immune_unassigned")

# Assuming metadata is already loaded in your environment
filtered_data <- filtered_data %>%
  left_join(metadata %>% select(sample_id, final_diagnosis_simple), by = "sample_id")

filtered_data <- filtered_data %>%
  filter(!is.na(final_diagnosis_simple)) %>%
  filter(!final_diagnosis_simple %in% c("GBM_other"))

# Filter the data where source == "immune"
filtered_data <- filtered_data %>%
  mutate(tumor_ag = if_else(
    cell_meta_cluster_final_broad == "Tumor" & 
      (EGFR_func == 1 | NG2_func == 1 | GM2_GD2_func == 1 | HER2_func == 1 | GPC2_func == 1 | B7H3_func == 1 | VISTA_func == 1),
    "Antigen_positive",
    "Antigen_negative"
  ))%>%
  mutate(cell_meta_cluster_final_broad = case_when(
    cell_meta_cluster_final == "Tumor_cells" & tumor_ag == "Antigen_positive" ~ "Antigen_positive",
    cell_meta_cluster_final == "Tumor_cells" & tumor_ag == "Antigen_negative" ~ "Antigen_negative",
    TRUE ~ cell_meta_cluster_final_broad  # Keep other values as they are
  ))

filtered_data <- filtered_data %>%
  filter(!cell_meta_cluster_final_broad %in% c('Undetermined'))

filtered_data <- filtered_data %>%
  mutate(cell_meta_cluster_ML = case_when(
    cell_meta_cluster_final %in% c("DC_Mac_CD209",
                                   "Macrophage_CD206",
                                   "Macrophage_CD68",
                                   "Macrophage_CD68_CD163",
                                   "Microglia_CD163",
                                   "Myeloid_CD11b_HLADRplus",
                                   "Myeloid_CD11b_HLADRminus",
                                   "Myeloid_CD14",
                                   "Myeloid_CD141",
                                   "Mast_cells",
                                   "APC",
                                   "Microglia",
                                   "Neutrophils",
                                   "Myeloid_CD14_CD163") ~ "myeloid",
    cell_meta_cluster_final %in% c("Tcell_CD4", 
                                   "Tcell_CD8", 
                                   "Tcell_FoxP3", 
                                   "Bcells") ~ "lymphoid",
    TRUE ~ "other"
  ))

filtered_data <- filtered_data %>%
  mutate(cell_meta_cluster_final_broad = case_when(
    cell_meta_cluster_ML == "myeloid" ~ "Myeloid",
    cell_meta_cluster_ML == "lymphoid" ~ "Lymphoid",
    TRUE ~ cell_meta_cluster_final_broad  # Keep other values as they are
  ))


# Rename cell_meta_cluster_final values
filtered_data <- filtered_data %>%
  mutate(cell_meta_cluster_final = case_when(
    cell_meta_cluster_final == "DC_Mac_CD209" ~ "DC/Mac CD209<sup>+</sup>",
    cell_meta_cluster_final == "Macrophage_CD206" ~ "Mac CD206<sup>+</sup>",
    cell_meta_cluster_final == "Macrophage_CD68" ~ "Mac CD68<sup>+</sup>",
    cell_meta_cluster_final == "Macrophage_CD68_CD163" ~ "Mac CD68<sup>+</sup>CD163<sup>+</sup>",
    cell_meta_cluster_final == "Microglia_CD163" ~ "Microglia CD163<sup>+</sup>",
    cell_meta_cluster_final == "Myeloid_CD11b_HLADRplus" ~ "Myeloid CD11b<sup>+</sup>HLADR<sup>+</sup>",
    cell_meta_cluster_final == "Myeloid_CD11b_HLADRminus" ~ "Myeloid CD11b<sup>+</sup>HLADR<sup>-</sup>",
    cell_meta_cluster_final == "Myeloid_CD14" ~ "Myeloid CD14<sup>+</sup>",
    cell_meta_cluster_final == "Myeloid_CD141" ~ "Myeloid CD141<sup>+</sup>",
    cell_meta_cluster_final == "Myeloid_CD14_CD163" ~ "Myeloid CD14<sup>+</sup>CD163<sup>+</sup>",
    cell_meta_cluster_final == "Tcell_CD4" ~ "T CD4<sup>+</sup>",
    cell_meta_cluster_final == "Tcell_CD8" ~ "T CD8<sup>+</sup>",
    cell_meta_cluster_final == "Tcell_FoxP3" ~ "T<sub>reg</sub>",
    cell_meta_cluster_final == "Mast_cells" ~ "Mast",
    cell_meta_cluster_final == "Bcells" ~ "B",
    cell_meta_cluster_final == "Endothelial_cells" ~ "Endothelial",
    #cell_meta_cluster_final == "Tumor_cells" ~ "Tumor",
    TRUE ~ cell_meta_cluster_final
  ))


# Summarize the data for broad and final categories
summary_data_broad <- filtered_data %>%
  count(cell_meta_cluster_final_broad) %>%
  rename(count_broad = n)

# Calculate the angles for the broad categories
summary_data_broad <- summary_data_broad %>%
  mutate(start_angle = lag(cumsum(count_broad) / sum(count_broad) * 360, default = 0),
         end_angle = cumsum(count_broad) / sum(count_broad) * 360)


# Define the desired order for the broad categories
desired_order_broad <- c(
  "Neurons",
  "Endothelial", "Lymphoid","Myeloid","Antigen_negative","Antigen_positive"
)

# Convert cell_meta_cluster_final_broad to a factor with the desired order and calculate percentage
summary_data_broad <- summary_data_broad %>%
  mutate(
    cell_meta_cluster_final_broad = factor(cell_meta_cluster_final_broad, levels = desired_order_broad),
    percentage = (count_broad / sum(count_broad)) * 100
  ) %>%
  arrange(cell_meta_cluster_final_broad)

# View the updated dataframe
summary_data_broad


color_broad <- c("#d2691e", "#c2b280", "#93003a",  "#85b7ce", "#696969","#d3d3d3")


pull_values_broad <- ifelse(summary_data_broad$cell_meta_cluster_final_broad %in% c("Lymphoid"), 
                            0.0, 
                            ifelse(summary_data_broad$cell_meta_cluster_final_broad %in% c("Myeloid"), 
                                   0.0, 
                                   0))

# Create the nested pie chart
pie <- plot_ly()

# Add inner pie for cell_meta_cluster_final_broad
pie <- pie %>%
  add_pie(data = summary_data_broad,
          labels = ~paste(cell_meta_cluster_final_broad, round(count_broad / sum(count_broad) * 100, 2), "%"),
          values = ~count_broad,
          hole = 0.5,
          textinfo = 'none',
          textposition = 'inside',
          sort = FALSE,
          pull = pull_values_broad,
          marker = list(colors = color_broad, line = list(color = 'white', width = 1)),
          domain = list(x = c(0, 1), y = c(0, 1)),
          showlegend = FALSE)

pie

# Save fig2 with transparent background
htmlwidgets::saveWidget(as_widget(pie), "fig2.html", selfcontained = FALSE)
webshot("fig2.html", file = file.path(output_dir, "broad_pie.png"), zoom = 10, vwidth = 600, vheight = 600, selector = "div.plot-container")

# Filter data for Immune broad category
immune_data <- filtered_data %>%
  filter(cell_meta_cluster_final_broad == "Myeloid") %>%
  count(cell_meta_cluster_final) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(n)) %>%
  mutate(cell_meta_cluster_final = factor(cell_meta_cluster_final, levels = cell_meta_cluster_final))


# Define the color palette
color_final <- c(
  "#b1dfdb",  # Myeloid CD14+
  "#a3c1d9",  # Mac CD68+
  "#1f4e79",  # Mac CD206+
  "#4a90d9",  # Myeloid CD11b+HLADR-
  "#6aaed6",  # Myeloid CD11b+HLADR+
  "#1f78b4",  # Mac CD68+CD163+
  "#4e78b5",  # APC
  "#66cdaa",  # Microglia
  "#3cb371",  # Microglia CD163+
  "#85b7ce",  # Myeloid CD14+CD163+
  "#6694c1",  # Myeloid CD141+
  "#00429d",  # DC/Mac CD209+
  "#32cd32",  # Neutrophils
  "#228b22"   # Mast
)

plot <- ggplot(immune_data, aes(x = "", y = n, fill = cell_meta_cluster_final)) +
  geom_bar(stat = "identity", position = "stack", width = 1, color = "black", size = 0.1) +  # Add borders with specified line width
  scale_fill_manual(values = color_final) +
  labs(x = "Fine", y = NULL, fill = NULL, title = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),  # Keep grid minimal
    plot.title = element_blank()
  )

plot

# Export the plot with adjusted dimensions to reduce empty space
ggsave(
  filename = file.path(output_dir, "Lymphoid_stacked.png"),
  plot = plot,
  width = 13,        # Reduce width
  height = 53,       # Reduce height
  units = "mm",
  dpi = 1200         # Keep high resolution
)


# Filter data for Immune broad category
immune_data <- filtered_data %>%
  filter(cell_meta_cluster_final_broad == "Lymphoid") %>%
  count(cell_meta_cluster_final) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(n)) %>%
  mutate(cell_meta_cluster_final = factor(cell_meta_cluster_final, levels = cell_meta_cluster_final))


color_final <- c(
  "#ffb3a7",  # T CD8+
  "#b81b4a",  # T CD4+
  "#93003a",  # Treg
  "#d5405e"   # B
)

plot <- ggplot(immune_data, aes(x = "", y = n, fill = cell_meta_cluster_final)) +
  geom_bar(stat = "identity", position = "stack", width = 1, color = "black", size = 0.1) +  # Add borders with specified line width
  scale_fill_manual(values = color_final) +
  labs(x = "Fine", y = NULL, fill = NULL, title = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),  # Keep grid minimal
    plot.title = element_blank()
  )

plot

# Export the plot with adjusted dimensions to reduce empty space
ggsave(
  filename = file.path(output_dir, "Lymphoid_stacked.png"),
  plot = plot,
  width = 13,        # Reduce width
  height = 53,       # Reduce height
  units = "mm",
  dpi = 1200         # Keep high resolution
)

# Panel d - Immune heatmap ----

# Create the extended immune feature list
cell_type_immune_extended <- paste0(cell_type_immune, "_over_all_immune_count_prop")

# Remove "Endothelial_cells_over_all_immune_count_prop", "Immune_unassigned_over_all_immune_count_prop", and "Unassigned_over_all_immune_count_prop" from the list
cell_type_immune_extended <- cell_type_immune_extended[!cell_type_immune_extended %in% c("Endothelial_cells_over_all_immune_count_prop", "Immune_unassigned_over_all_immune_count_prop", "Unassigned_over_all_immune_count_prop")]

# Filter, modify, reshape, and sort the data
heatmap_data <- master_feature_table_filtered %>%
  filter(!is.na(who_diagnosis)) %>%
  filter(feature_variable %in% cell_type_immune_extended) %>%
  filter(immunotherapy == "no") %>%
  pivot_wider(names_from = feature_variable, values_from = feature_value) %>%
  group_by(who_diagnosis) %>%
  summarise(across(all_of(cell_type_immune_extended), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(who_diagnosis)

# Prepare data for ggplot
heatmap_long <- heatmap_data %>%
  pivot_longer(cols = -who_diagnosis, names_to = "feature_variable", values_to = "average_abundance")

# Apply power transformation to average_abundance
heatmap_long <- heatmap_long %>%
  mutate(scaled_abundance = average_abundance^0.5)

# Create data_matrix_z
data_matrix <- as.matrix(heatmap_data[, -1])
data_matrix_z <- (scale((data_matrix)))
colnames(data_matrix_z) <- colnames(data_matrix)
rownames(data_matrix_z) <- gsub("_over_all_immune_count_prop", "", rownames(data_matrix_z))
rownames(data_matrix_z) <- gsub("_", " ", (heatmap_data$who_diagnosis))

# Convert the z-scored matrix to a long format for ggplot
heatmap_z_long <- as.data.frame(as.table(data_matrix_z))
colnames(heatmap_z_long) <- c("who_diagnosis", "feature_variable", "z_score_average_abundance")

# Clean column names to ensure they match
colnames(heatmap_long) <- gsub(" ", "_", colnames(heatmap_long))
colnames(heatmap_z_long) <- gsub(" ", "_", colnames(heatmap_z_long))

# Merge the abundance and z-score data
tumor_ab_in_merged_heatmap <- merge(heatmap_long, heatmap_z_long, by = c("feature_variable", "who_diagnosis"))

tumor_ab_in_merged_heatmap$feature_variable <- gsub("_over_all_immune_count_prop", "", tumor_ab_in_merged_heatmap$feature_variable)
#tumor_ab_in_merged_heatmap$feature_variable <- gsub("_", " ", tumor_ab_in_merged_heatmap$feature_variable)

# Define the desired order
desired_order <- c(
  "4 - GBM",
  "4 - Astrocytoma",
  "3 - Astrocytoma",
  "2 - Astrocytoma",
  "3 - Oligodendroglioma",
  "2 - Oligodendroglioma",
  "3 - PXA",
  "2 - PXA",
  "Pediatric DIPG",
  "Pediatric HGG (other)"
)

# Set who_diagnosis as a factor with the specified levels
tumor_ab_in_merged_heatmap$who_diagnosis <- factor(tumor_ab_in_merged_heatmap$who_diagnosis, levels = desired_order)


# Update the feature_variable with the correct feature names
tumor_ab_in_merged_heatmap <- tumor_ab_in_merged_heatmap %>%
  mutate(feature_variable = case_when(
    feature_variable == "DC_Mac_CD209" ~ "DC/Mac CD209+",
    feature_variable == "Macrophage_CD206" ~ "Mac CD206+",
    feature_variable == "Macrophage_CD68" ~ "Mac CD68+",
    feature_variable == "Macrophage_CD68_CD163" ~ "Mac CD68+CD163+",
    feature_variable == "Microglia_CD163" ~ "Microglia CD163+",
    feature_variable == "Myeloid_CD11b_HLADRplus" ~ "Myeloid CD11b+HLADR+",
    feature_variable == "Myeloid_CD11b_HLADRminus" ~ "Myeloid CD11b+HLADR-",
    feature_variable == "Myeloid_CD14" ~ "Myeloid CD14+",
    feature_variable == "Myeloid_CD141" ~ "Myeloid CD141+",
    feature_variable == "Myeloid_CD14_CD163" ~ "Myeloid CD14+CD163+",
    feature_variable == "Tcell_CD4" ~ "T CD4+",
    feature_variable == "Tcell_CD8" ~ "T CD8+",
    feature_variable == "Tcell_FoxP3" ~ "Treg",
    feature_variable == "Mast_cells" ~ "Mast",
    feature_variable == "Bcells" ~ "B",
    feature_variable == "Endothelial_cells" ~ "Endothelial",
    TRUE ~ feature_variable
  ))


# Define the mapping for feature_variable to expressions using bquote
feature_variable_expressions <- c(
  "DC/Mac CD209+" = bquote("DC/Mac CD209"^"+"),
  "Mac CD206+" = bquote("Mac CD206"^"+"),
  "Mac CD68+" = bquote("Mac CD68"^"+"),
  "Mac CD68+CD163+" = bquote("Mac CD68"^"+" ~ "CD163"^"+"),
  "Microglia CD163+" = bquote("Microglia CD163"^"+"),
  "Myeloid CD11b+HLADR+" = bquote("Myeloid CD11b"^"+" ~ "HLADR"^"+"),
  "Myeloid CD11b+HLADR-" = bquote("Myeloid CD11b"^"+" ~ "HLADR"^"-"),
  "Myeloid CD14+" = bquote("Myeloid CD14"^"+"),
  "Myeloid CD141+" = bquote("Myeloid CD141"^"+"),
  "Myeloid CD14+CD163+" = bquote("Myeloid CD14"^"+" ~ "CD163"^"+"),
  "T CD4+" = bquote("T CD4"^"+"),
  "T CD8+" = bquote("T CD8"^"+"),
  "Treg" = bquote("Treg"),
  "Mast" = bquote("Mast"),
  "B" = bquote("B"),
  "Endothelial" = bquote("Endothelial")
)


# Set the levels in the order you want
tumor_ab_in_merged_heatmap$feature_variable <- factor(
  tumor_ab_in_merged_heatmap$feature_variable,
  levels = unique(tumor_ab_in_merged_heatmap$feature_variable)
)

# Prepare the data matrix for clustering
heatmap_matrix <- tumor_ab_in_merged_heatmap %>%
  dplyr::select(feature_variable, who_diagnosis, z_score_average_abundance) %>%
  pivot_wider(names_from = who_diagnosis, values_from = z_score_average_abundance) %>%
  column_to_rownames("feature_variable")

# Convert to a matrix
heatmap_matrix <- as.matrix(heatmap_matrix)

# Perform hierarchical clustering
hc <- hclust(dist(heatmap_matrix, method = "euclidean"), method = "complete")

# Extract the order of rows from the clustering
ordered_features <- rownames(heatmap_matrix)[hc$order]

# Reorder the feature_variable column in tumor_ab_in_merged_heatmap
tumor_ab_in_merged_heatmap$feature_variable <- factor(tumor_ab_in_merged_heatmap$feature_variable, levels = ordered_features)

# Convert the dendrogram to a ggplot-friendly format and rotate it
dend_data <- as.dendrogram(hc) %>% dendro_data()

# Adjust the dendrogram segments to flip it correctly
dend_data$segments <- dend_data$segments %>% 
  mutate(y = -y, yend = -yend)  # Flip the y coordinates

# Create the ggplot dendrogram with adjusted height
dend_plot <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.15) +  # Adjust line size for thinner lines
  coord_flip() +  # Rotate the dendrogram
  scale_x_reverse() +  # Reverse the x-axis to face the heatmap
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.45, 0, 0.2, 0.1), "cm")  # Adjust margins to align with heatmap
  )

# Rename the who_diagnosis column values to the new desired format
tumor_ab_in_merged_heatmap$who_diagnosis <- recode(tumor_ab_in_merged_heatmap$who_diagnosis,
                                                   "4 - GBM" = "GBM - 4",
                                                   "4 - Astrocytoma" = "Astrocytoma - 4",
                                                   "3 - Astrocytoma" = "Astrocytoma - 3",
                                                   "2 - Astrocytoma" = "Astrocytoma - 2",
                                                   "3 - Oligodendroglioma" = "Oligodendroglioma - 3",
                                                   "2 - Oligodendroglioma" = "Oligodendroglioma - 2",
                                                   "3 - PXA" = "PXA - 3",
                                                   "2 - PXA" = "PXA - 2",
                                                   "Pediatric DIPG" = "Pediatric DIPG",
                                                   "Pediatric HGG (other)" = "Pediatric HGG (other)")

heatmap_plot <- ggplot(tumor_ab_in_merged_heatmap, aes(x = who_diagnosis, y = feature_variable)) +
  geom_tile(aes(width = 1, height = 1), color = "black", linewidth = 0.1, fill = NA) +  # Add boxes around each point
  geom_point(
    aes(size = scaled_abundance, fill = z_score_average_abundance),
    shape = 21,  # Shape 21 allows for both fill and border
    color = "black",  # Outline color
    stroke = 0.1      # Outline thickness
  ) +
  scale_size_area(max_size = 4) +
  scale_fill_gradient2(
    low = "#00429D", mid = "white", high = "#B81B4A", midpoint = 0, limits = c(-1.5, 1.5),
    breaks = c(-1.5, 0, 1.5),  # Set the breaks at the limits and midpoint
    labels = c("-1.5", "0", "1.5"),  # Label the breaks
    oob = scales::squish  # Cap the values beyond the limits
  ) +
  labs(
    x = " ",
    y = " ",
    title = " "
  ) +
  scale_y_discrete(labels = feature_variable_expressions, position = "right") +  # Apply expressions to y-axis labels
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 6),  # Apply Helvetica and size 6 to all text
    axis.text.x = element_blank(),  # Rotate and resize x-axis text
    axis.text.y = element_blank(),  
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = "none",  # Remove the legend by setting position to "none"
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Adjust margins to align with dendrogram
  )

# Combine the dendrogram and heatmap plots with adjusted height
combined_plot <- grid.arrange(dend_plot, heatmap_plot, ncol = 2, widths = c(0.5, 4))

# Display the combined plot
print(combined_plot)


# Save the heatmap plot
ggsave(
  filename = file.path(output_dir, "heatmap_immune_cells.png"),
  plot = combined_plot,
  width = 55,
  height = 90,
  units = "mm",
  dpi = 1200
)

# Save the data as CSV
write_csv(tumor_ab_in_merged_heatmap, file.path(output_dir, "tumor_ab_in_merged_heatmap.csv"))

# Panel f - Tumor heatmap ----


# Filter the dataframe and summarize by patient_id and feature_variable, selecting only the necessary columns
abundance_df <- master_feature_table_filtered %>%
  filter(feature_variable %in% tumor_func_columns_extended) %>%
  select(patient_id, feature_variable, feature_value) %>%
  group_by(patient_id, feature_variable) %>%
  summarise(feature_value = mean(feature_value, na.rm = TRUE)) %>%
  rename(abundance = feature_value) %>%
  mutate(feature_variable = str_replace(feature_variable, "_func_over_all_tumor_count_prop", ""))



# Filter the dataframe and summarize by patient_id and feature_variable, selecting only the necessary columns
intensity_df <- master_feature_table_filtered %>%
  filter(cell_meta_cluster_final == "Tumor_cells") %>%
  filter(feature_variable %in% tumor_columns_core) %>%
  select(patient_id, feature_variable, feature_value) %>%
  group_by(patient_id, feature_variable) %>%
  summarise(feature_value = mean(feature_value, na.rm = TRUE))%>%
  rename(intensity = feature_value)

# Perform a full join on 'patient_id' and 'source'
tumor_ab_in_merged <- full_join(abundance_df, intensity_df, by = c("patient_id", "feature_variable"))

tumor_ab_in_merged <- tumor_ab_in_merged %>%
  left_join(metadata, by = "patient_id")

# Combine final_diagnosis and who_grade into a new factor variable
tumor_ab_in_merged_heatmap <- tumor_ab_in_merged %>%
  group_by(feature_variable, who_diagnosis) %>%
  summarize(
    average_intensity = mean(intensity, na.rm = TRUE),
    average_abundance = mean(abundance, na.rm = TRUE)
  ) %>%
  ungroup()  %>%
  mutate(
    who_diagnosis = factor(who_diagnosis, levels = desired_order_who_diagnosis)
  ) %>%
  group_by(feature_variable) %>%
  mutate(
    z_score_intensity = (average_intensity - mean(average_intensity, na.rm = TRUE)) / sd(average_intensity, na.rm = TRUE)
  ) %>%
  ungroup()


tumor_ab_in_merged_heatmap <- tumor_ab_in_merged_heatmap[!is.na(tumor_ab_in_merged_heatmap$who_diagnosis), ]

tumor_ab_in_merged_heatmap$feature_variable <- gsub('GM2_GD2', 'GM2/GD2', tumor_ab_in_merged_heatmap$feature_variable)

# Prepare the data matrix for clustering
heatmap_matrix <- tumor_ab_in_merged_heatmap %>%
  dplyr::select(feature_variable, who_diagnosis, average_abundance) %>%
  pivot_wider(names_from = who_diagnosis, values_from = average_abundance) %>%
  column_to_rownames("feature_variable")

# Convert to a matrix
heatmap_matrix <- as.matrix(heatmap_matrix)

# Perform hierarchical clustering
hc <- hclust(dist(heatmap_matrix, method = "euclidean"), method = "complete")

# Extract the order of rows from the clustering
ordered_features <- rownames(heatmap_matrix)[hc$order]

# Reorder the feature_variable column in tumor_ab_in_merged_heatmap
tumor_ab_in_merged_heatmap$feature_variable <- factor(tumor_ab_in_merged_heatmap$feature_variable, levels = ordered_features)

# Convert the dendrogram to a ggplot-friendly format and rotate it
dend_data <- as.dendrogram(hc) %>% dendro_data()

# Adjust the dendrogram segments to flip it correctly
dend_data$segments <- dend_data$segments %>% 
  mutate(y = -y, yend = -yend)  # Flip the y coordinates

# Create the ggplot dendrogram with adjusted height
dend_plot <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.15) +  # Adjust line size for thinner lines
  coord_flip() +  # Rotate the dendrogram
  scale_x_reverse() +  # Reverse the x-axis to face the heatmap
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.55, 0, 0.45, 0.1), "cm")  # Adjust margins to align with heatmap
  )

# Rename the who_diagnosis column values to the new desired format
tumor_ab_in_merged_heatmap$who_diagnosis <- recode(tumor_ab_in_merged_heatmap$who_diagnosis,
                                                   "4 - GBM" = "GBM",
                                                   "4 - Astrocytoma" = "Astrocytoma - 4",
                                                   "3 - Astrocytoma" = "Astrocytoma - 3",
                                                   "2 - Astrocytoma" = "Astrocytoma - 2",
                                                   "3 - Oligodendroglioma" = "Oligodendroglioma - 3",
                                                   "2 - Oligodendroglioma" = "Oligodendroglioma - 2",
                                                   "3 - PXA" = "PXA - 3",
                                                   "2 - PXA" = "PXA - 2",
                                                   "Pediatric DIPG" = "Pediatric DIPG",
                                                   "Pediatric HGG (other)" = "Pediatric HGG (other)")


heatmap_plot <- ggplot(tumor_ab_in_merged_heatmap, aes(x = who_diagnosis, y = feature_variable)) +
  geom_tile(aes(width = 1, height = 1), color = "black", linewidth = 0.1, fill = NA) +  # Add boxes around each point
  geom_point(
    aes(size = average_abundance, fill = z_score_intensity),
    shape = 21,  # Shape 21 allows for both fill and border
    color = "black",  # Outline color
    stroke = 0.1      # Outline thickness
  ) +
  scale_size_area(max_size = 4) +
  scale_fill_gradient2(
    low = "#00429D", mid = "white", high = "#B81B4A", midpoint = 0, limits = c(-1.5, 1.5),
    breaks = c(-1.5, 0, 1.5),  # Set the breaks at the limits and midpoint
    labels = c("-1.5", "0", "1.5"),  # Label the breaks
    oob = scales::squish  # Cap the values beyond the limits
  ) +
  labs(
    x = " ",
    y = " ",
    title = " "
  ) +
  scale_y_discrete(labels = feature_variable_expressions, position = "right") +  # Apply expressions to y-axis labels
  #scale_size_area(max_size = 10, breaks = c(0.1, 0.3, 0.5, 0.7)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 6),  # Apply Helvetica and size 6 to all text
    axis.text.x = element_blank(),  # Rotate and resize x-axis text
    axis.text.y = element_blank(),  # Resize y-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = "none",  # Remove legend
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Adjust margins to align with dendrogram
  )


# Combine the dendrogram and heatmap plots with adjusted height
combined_plot <- grid.arrange(dend_plot, heatmap_plot, ncol = 2, widths = c(0.5, 4))

# Display the combined plot
print(combined_plot)

# Save the heatmap plot
ggsave(
  filename = file.path(output_dir, "heatmap_tumor_cells.png"),
  plot = combined_plot,
  width = 55,
  height = 40,
  units = "mm",
  dpi = 1200
)

annotation_colors <- c(
  "Astrocytoma - 2" = "#B1DFDB",  # Light Sky Blue
  "Oligodendroglioma - 2" = "#B1DFDB",  # Light Sky Blue
  "PXA - 2" = "#B1DFDB",  # Light Sky Blue
  "Astrocytoma - 3" = "#618FBF",  # Dodger Blue
  "Oligodendroglioma - 3" = "#618FBF",  # Dodger Blue
  "PXA - 3" = "#618FBF",  # Dodger Blue
  "Astrocytoma - 4" = "#00429D",  # Dark Blue
  "GBM" = "#00429D",  # Dark Blue
  "Pediatric DIPG" = "white",
  "Pediatric HGG (other)" = "white"
)

# Create the annotation plot
annotation_plot <- ggplot(tumor_ab_in_merged_heatmap, aes(x = who_diagnosis, y = 1)) +
  geom_tile(aes(fill = who_diagnosis), color = "black", linewidth =0.2 ) +
  scale_fill_manual(values = annotation_colors) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

annotation_plot

ggsave(
  filename = "annotation_plot.png",
  plot = annotation_plot,
  device = "png",
  width = 48,  # width in mm
  height = 3,  # height in mm
  units = "mm"
)


###### Figure 3 ################################################################
# Panel a - Shannon diversity index ----

# Define the extended tumor function columns
tumor_func_columns_extended <- c(
  "B7H3_func_over_all_tumor_count_prop", "EGFR_func_over_all_tumor_count_prop",
  "GM2_GD2_func_over_all_tumor_count_prop", "GPC2_func_over_all_tumor_count_prop",
  "HER2_func_over_all_tumor_count_prop", "NG2_func_over_all_tumor_count_prop",
  "VISTA_func_over_all_tumor_count_prop"
)

# Function to calculate Shannon diversity index
calculate_shannon_index <- function(proportions) {
  proportions <- proportions[proportions > 0] # Remove zero proportions
  -sum(proportions * log(proportions))
}

# Filter, pivot, and calculate Shannon diversity index
diversity_data <- master_feature_table_filtered %>%
  filter(feature_variable %in% tumor_func_columns_extended) %>%
  select(sample_id, feature_variable, feature_value) %>%
  spread(key = feature_variable, value = feature_value) %>%
  rowwise() %>%
  mutate(shannon_index = calculate_shannon_index(c_across(all_of(tumor_func_columns_extended))))


# Display the resulting data
diversity_data %>%
  dplyr::select(sample_id, shannon_index) %>%
  ungroup() %>%
  as.data.frame() %>%
  head()

tumor_antigen_segment_df_metadata <- combine_data_metadata(diversity_data,metadata,"sample_id")
# Merge the diversity index data with the original dataframe to get final_diagnosis_simple
merged_data <- tumor_antigen_segment_df_metadata %>%
  dplyr::filter(immunotherapy == "no") %>%
  dplyr::select(sample_id, final_diagnosis_simple) %>%
  dplyr::distinct() %>%
  dplyr::right_join(diversity_data, by = "sample_id") %>%
  dplyr::filter(!is.na(final_diagnosis_simple))

# Define the order for final_diagnosis_simple
diagnosis_order <- c("GBM", "Astrocytoma", "Oligodendroglioma", "PXA", "Pediatric DIPG", "Pediatric HGG (other)")

# Convert final_diagnosis_simple to a factor with the specified order
merged_data <- merged_data %>%
  dplyr::mutate(final_diagnosis_simple = factor(final_diagnosis_simple, levels = diagnosis_order))

merged_data <- merged_data %>%
  dplyr::filter(!is.na(final_diagnosis_simple))

# Define colors for each diagnosis
diagnosis_colors <- c("GBM" = "#c52a52", "Astrocytoma" = "#325da9", "Oligodendroglioma" = "#228b22", 
                      "PXA" = "#faa51e", "Pediatric DIPG" = "#b2dfdb", "Pediatric HGG (other)" = "#505150")

p <- ggplot(merged_data, aes(x = final_diagnosis_simple, y = shannon_index, fill = final_diagnosis_simple)) +
  geom_boxplot(color = "black", outlier.shape = NA, size = 0.25) +
  geom_jitter(width = 0.2, size = 0.05, alpha = 1, color = "black") +
  scale_fill_manual(values = diagnosis_colors) +
  labs(title = "",
       x = NULL,
       y = NULL) +
  theme_minimal() +
  theme(text = element_text(family = "Helvetica", size = 6, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, family = "Helvetica", color = "black"),
        axis.title.y = element_text(size = 6, family = "Helvetica", color = "black"),
        plot.title = element_text(size = 6, family = "Helvetica", color = "black"),
        panel.grid = element_blank(),     # Remove grid lines
        axis.line = element_line(color = "black", linewidth = 0.25),  # Add axis lines
        axis.ticks = element_line(color = "black", linewidth = 0.15),
        axis.ticks.length = unit(0.5, "mm"),
        legend.position = "none")  # This line removes the legend

p

ggsave(
  filename = "shannon_diversity_index_plot.pdf",
  plot = p,
  height = 33,
  width = 52,
  units = "mm"
)

# Panel b - tumor ag coverage  -----

column_of_interest <- "final_diagnosis_simple"
col_order <- c("PXA", "Pediatric DIPG", "Pediatric HGG (other)")
col_colors <- c("PXA" = "#faa51e", "Pediatric DIPG" = "#b2dfdb", "Pediatric HGG (other)" = "#505150")  # Example colors

generate_cumulative_line_plot(
  data = data_file,
  column_of_interest = column_of_interest,
  col_order = col_order,
  col_colors = col_colors,
  output_filename = file.path(output_dir, "no_change_tumor_coverage.pdf"),  # Use output_dir for file path
  width = 45,   # Width in millimeters
  height = 35   # Height in millimeters
)

# Panel c - tumor ag coverage  -----

col_order <- c("GBM", "Astrocytoma", "Oligodendroglioma")
col_colors <- c("GBM" = "#c52a52", "Astrocytoma" = "#325da9", "Oligodendroglioma" = "#228b22")  # Example colors

generate_cumulative_line_plot(
  data = data_file,
  column_of_interest = column_of_interest,
  col_order = col_order,
  col_colors = col_colors,
  output_filename = file.path(output_dir, "change_tumor_coverage.pdf"),  # Use output_dir for file path
  width = 45,   # Width in millimeters
  height = 35   # Height in millimeters
)



# Panel g - antigen negative tumor cells violin plot  -----
col_order <- c("GBM", "Astrocytoma", "Oligodendroglioma","PXA", "Pediatric DIPG", "Pediatric HGG (other)")
col_colors <- c(
  "GBM" = "#c52a52", "Astrocytoma" = "#325da9", "Oligodendroglioma" = "#228b22","PXA" = "#faa51e", "Pediatric DIPG" = "#b2dfdb", "Pediatric HGG (other)" = "#505150"
)  # Customize colors as needed


plot <- generate_violin_plot(
  data = data_file,
  column_of_interest = column_of_interest,
  col_order = col_order,
  col_colors = col_colors,
  output_filename = file.path(output_dir, "violin_plot.pdf"),  # Use output_dir for file path
  width = 55,   # Width in millimeters
  height = 32   # Height in millimeters
)

# Display the plot
print(plot)




# Panel d - top tumor ag for coverage, % patient heatmap ----

# Example usage
generate_marker_pair_summary_heatmap(
  data_file = data_file,
  metadata_file = metadata,
  column_names = c("B7H3", "EGFR", "HER2", "NG2", "GM2_GD2", "GPC2", "VISTA"),
  column_of_interest = "final_diagnosis_simple",
  col_order = c("GBM", "Astrocytoma", "Oligodendroglioma")
)

# Panel e - python notebook - Figure_Generation/Figure_3_panel_e.ipynb ----
###### Figure 4 ################################################################
# Panel a - b dyanmic range /mean intensity/coverage bubble plot ----

cell_table_dt <- as.data.table(cell_table_all_merged_thresholded)

# calculated for panel a
fold_change_dt <- combine_data_metadata(fold_change_dt, metadata, merge_by = "sample_id", remove_duplicates = FALSE)
cell_table_dt <- combine_data_metadata(cell_table_dt, metadata, merge_by = "sample_id", remove_duplicates = FALSE)

# Define tumor columns
tumor_columns_core <- c("B7H3", "EGFR", "GM2_GD2", "GPC2", "HER2", "NG2", "VISTA")

# Convert to data.table if not already
setDT(fold_change_dt)

dynamic_ranges <- cell_table_dt[
  , lapply(.SD, function(x) {
    valid_values <- x[x != 100.199]  # Exclude 100.199
    quantile(valid_values, 0.95, na.rm = TRUE) - quantile(valid_values, 0.05, na.rm = TRUE)
  }), 
  by = .(patient_id, final_diagnosis_simple), 
  .SDcols = tumor_columns_core
]

dynamic_ranges <- dynamic_ranges[!is.na(final_diagnosis_simple) & final_diagnosis_simple != "GBM_other"]


# Step 3: Y-axis - Calculate the average dynamic range within each patient for each tumor marker
y_axis_avg_within_patient <- dynamic_ranges[, lapply(.SD, mean, na.rm = TRUE), 
                                            by = .(patient_id, final_diagnosis_simple), 
                                            .SDcols = tumor_columns_core]

# Step 4: Summarize the average within-patient dynamic ranges per final_diagnosis_simple for each tumor marker
y_axis_dynamic_range <- y_axis_avg_within_patient[, lapply(.SD, mean, na.rm = TRUE), 
                                                  by = final_diagnosis_simple, 
                                                  .SDcols = tumor_columns_core]

coverage <- master_feature_table_filtered %>% filter(feature_type == "Cell_Abundance", bio_feature_type == "Relative_to_all_tumor_cells", feature_variable == "B7H3_func_over_all_tumor_count_prop")


coverage <- coverage %>%
  filter(tumor_region %in% c("tumor_core", "other"))


# Summarizing coverage by patient_id
coverage_summary <- coverage %>%
  group_by(patient_id) %>%
  summarise(coverage_value = mean(feature_value, na.rm = TRUE))

# Summarizing antigen_intensity by patient_id for B7H3 ---


antigen_intensity_summary <- antigen_intensity %>%
  filter(feature_variable == "B7H3", feature_value != 100.199) %>%
  group_by(patient_id) %>%
  summarise(intensity_value = mean(feature_value, na.rm = TRUE))

# Merging the datasets
merged_df <- y_axis_avg_within_patient %>%
  inner_join(coverage_summary, by = "patient_id") %>%
  inner_join(antigen_intensity_summary, by = "patient_id")

# Filter out 'Pediatric DIPG' and 'Pediatric HGG (other)' from final_diagnosis_simple
filtered_df <- merged_df %>%
  filter(!final_diagnosis_simple %in% c("Pediatric DIPG", "Pediatric HGG (other)", "PXA"))


# Custom function to generate y-axis labels as 1X, 2X, etc.
custom_log_labels <- function(breaks) {
  min_value <- min(filtered_df$intensity_value, na.rm = TRUE)  # Get the minimum break value
  sapply(breaks, function(x) paste0(round(x / min_value, 1), "X"))  # Create labels
}

# Determine the maximum y-value and expand by a factor
max_y <- max(filtered_df$intensity_value, na.rm = TRUE) * 11

# Calculate the corresponding y-values for 1X and 50X
min_value <- min(filtered_df$intensity_value, na.rm = TRUE)
y_1x <- min_value  # 1X
y_50x <- min_value * 50  # 50X

# Generate log-scaled breaks, including 1X and 50X positions
y_breaks <- c(min_value, min_value * 10, min_value * 50, max_y)

# remove and patients where <10% of cells were identified as tumor


# Assuming metadata is your dataframe containing patient_id and sample_id
tumor_count_result_summarized <- tumor_count_result %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(patient_id) %>%
  summarise(all_tumor_count = sum(all_tumor_count, na.rm = TRUE))

# View the summarized result
print(tumor_count_result_summarized)


# Assuming `metadata` contains `patient_id` and `sample_id`
all_cell_count_tumor_FOV_result_summarized <- all_cell_count_tumor_FOV_result %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(patient_id) %>%
  summarise(all_cell_count_tumor_FOV = sum(all_cell_count_tumor_FOV, na.rm = TRUE))

# View the summarized result
print(all_cell_count_tumor_FOV_result_summarized)


# Assuming the two summarized dataframes and filtered_df are available
filtered_df <- filtered_df %>%
  left_join(all_cell_count_tumor_FOV_result_summarized, by = "patient_id") %>%
  left_join(tumor_count_result_summarized, by = "patient_id")

# View the updated filtered_df
print(filtered_df)

# Add the perc_tumor column and filter
filtered_df <- filtered_df %>%
  mutate(perc_tumor = (all_tumor_count / all_cell_count_tumor_FOV) * 100) %>%
  rowwise() %>% 
  mutate(removed = ifelse(perc_tumor < 10,TRUE,FALSE))

removed_df <- filtered_df %>% 
  filter(removed==TRUE) 

retained_filtered <- filtered_df %>% 
  filter(!removed)

# View the unique patient_id and perc_tumor values for removed rows
print(removed_df %>% select(patient_id, perc_tumor))


ggplot(filtered_df, aes(x = coverage_value, y = intensity_value, color = final_diagnosis_simple)) +
  geom_point(aes(size = B7H3), alpha = 0.7) +
  geom_hline(yintercept = y_1x, linetype = "dashed", color = "gray") +  # 1X line
  geom_hline(yintercept = y_50x, linetype = "dashed", color = "gray") +  # 50X line
  geom_vline(xintercept = 0.72, linetype = "dotted", color = "black") +  # Vertical line at 0.72
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +     # Vertical line at 1
  scale_size_continuous(range = c(0.05, 4)) +
  xlim(0, 1) +
  labs(
    x = '', 
    y = '', 
    title = '', 
    color = 'Final Diagnosis', size = 'Antigen Intensity'
  ) +
  scale_color_manual(values = c("GBM" = "#c52a52", "Astrocytoma" = "#325da9", "Oligodendroglioma" = "#228b22")) +
  scale_y_log10(
    limits = c(min_value / 2, max_y),  # Properly set y-axis limits
    breaks = y_breaks,  # Set breaks to align with lines
    labels = custom_log_labels(y_breaks)  # Apply correct labels
  ) +
  theme_minimal(base_size = 6, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_text(color = "black", size = 6),
    axis.text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 6, hjust = 0.5)
  )


# Export as SVG with the specified size
ggsave(
  filename = file.path(output_dir, "scatter_plot_B7H3_linear.svg"),
  width = 100,
  height = 50,
  units = "mm"
)

# Summarizing antigen_intensity by patient_id for EGFR ---

coverage <- master_feature_table_filtered %>% filter(feature_type == "Cell_Abundance", bio_feature_type == "Relative_to_all_tumor_cells", feature_variable == "EGFR_func_over_all_tumor_count_prop")

coverage <- coverage %>%
  filter(tumor_region %in% c("tumor_core", "other"))

# Summarizing coverage by patient_id
coverage_summary <- coverage %>%
  group_by(patient_id) %>%
  summarise(coverage_value = mean(feature_value, na.rm = TRUE))


antigen_intensity_summary <- antigen_intensity %>%
  filter(feature_variable == "EGFR", feature_value != 100.199) %>%
  group_by(patient_id) %>%
  summarise(intensity_value = mean(feature_value, na.rm = TRUE))

# Merging the datasets
merged_df <- y_axis_avg_within_patient %>%
  inner_join(coverage_summary, by = "patient_id") %>%
  inner_join(antigen_intensity_summary, by = "patient_id")

# Filter out 'Pediatric DIPG' and 'Pediatric HGG (other)' from final_diagnosis_simple
filtered_df <- merged_df %>%
  filter(!final_diagnosis_simple %in% c("Pediatric DIPG", "Pediatric HGG (other)", "PXA"))


# Custom function to generate y-axis labels as 1X, 2X, etc.
custom_log_labels <- function(breaks) {
  min_value <- min(filtered_df$intensity_value, na.rm = TRUE)  # Get the minimum break value
  sapply(breaks, function(x) paste0(round(x / min_value, 1), "X"))  # Create labels
}

# Determine the maximum y-value and expand by a factor
max_y <- max(filtered_df$intensity_value, na.rm = TRUE) * 2.7

# Calculate the corresponding y-values for 1X and 50X
min_value <- min(filtered_df$intensity_value, na.rm = TRUE)
y_1x <- min_value  # 1X
y_50x <- min_value * 50  # 50X

# Generate log-scaled breaks, including 1X and 50X positions
y_breaks <- c(min_value, min_value * 10, min_value * 50, max_y)


# remove and patients where <10% of cells were identified as tumor

# Assuming metadata is your dataframe containing patient_id and sample_id
tumor_count_result_summarized <- tumor_count_result %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(patient_id) %>%
  summarise(all_tumor_count = sum(all_tumor_count, na.rm = TRUE))

# View the summarized result
print(tumor_count_result_summarized)


# Assuming `metadata` contains `patient_id` and `sample_id`
all_cell_count_tumor_FOV_result_summarized <- all_cell_count_tumor_FOV_result %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(patient_id) %>%
  summarise(all_cell_count_tumor_FOV = sum(all_cell_count_tumor_FOV, na.rm = TRUE))

# View the summarized result
print(all_cell_count_tumor_FOV_result_summarized)


# Assuming the two summarized dataframes and filtered_df are available
filtered_df <- filtered_df %>%
  left_join(all_cell_count_tumor_FOV_result_summarized, by = "patient_id") %>%
  left_join(tumor_count_result_summarized, by = "patient_id")

# View the updated filtered_df
print(filtered_df)

# Add the perc_tumor column and filter
filtered_df <- filtered_df %>%
  mutate(perc_tumor = (all_tumor_count / all_cell_count_tumor_FOV) * 100) %>%
  rowwise() %>% 
  mutate(removed = ifelse(perc_tumor < 10,TRUE,FALSE))

removed_df <- filtered_df %>% 
  filter(removed==TRUE) 

retained_filtered <- filtered_df %>% 
  filter(!removed)

# View the unique patient_id and perc_tumor values for removed rows
print(removed_df %>% select(patient_id, perc_tumor))

ggplot(filtered_df, aes(x = coverage_value, y = intensity_value, color = final_diagnosis_simple)) +
  geom_point(aes(size = EGFR), alpha = 0.7) +
  geom_hline(yintercept = y_1x, linetype = "dashed", color = "gray") +  # 1X line
  geom_hline(yintercept = y_50x, linetype = "dashed", color = "gray") +  # 50X line
  geom_vline(xintercept = 0.72, linetype = "dotted", color = "black") +  # Vertical line at 0.72
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +     # Vertical line at 1
  scale_size_continuous(range = c(0.05, 4)) +
  xlim(0, 1) +
  labs(
    x = '', 
    y = '', 
    title = '', 
    color = 'Final Diagnosis', size = 'Antigen Intensity'
  ) +
  scale_color_manual(values = c("GBM" = "#c52a52", "Astrocytoma" = "#325da9", "Oligodendroglioma" = "#228b22")) +
  scale_y_log10(
    limits = c(min_value / 2, max_y),  # Properly set y-axis limits
    breaks = y_breaks,  # Set breaks to align with lines
    labels = custom_log_labels(y_breaks)  # Apply correct labels
  ) +
  theme_minimal(base_size = 6, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_text(color = "black", size = 6),
    axis.text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 6, hjust = 0.5)
  )

# Export as SVG with the specified size
ggsave(
  filename = file.path(output_dir, "scatter_plot_EGFR_linear.svg"),
  width = 100,
  height = 50,
  units = "mm"
)



# Panel c - tumor ag dynamic range single cell plot -----

# Define the tumor columns
tumor_columns_core <- c("B7H3", "EGFR", "GM2_GD2", "GPC2", "HER2", "NG2", "VISTA")

# Convert the data frame to a data.table
cell_table_dt <- as.data.table(cell_table_all_merged_thresholded)

# Step 1: Calculate the median for each marker after excluding 0 and 100.199
medians <- sapply(tumor_columns_core, function(col) {
  valid_values <- cell_table_dt[[col]][cell_table_dt[[col]] != 0 & cell_table_dt[[col]] != 100.199]
  median(valid_values, na.rm = TRUE)
})
print(medians)
# Name the medians correctly
names(medians) <- tumor_columns_core

# Step 2: Create a new fold change dataframe including sample_id, label, and segment columns
fold_change_dt <- cell_table_dt[, c("sample_id", "label", tumor_columns_core), with = FALSE]

# Exclude rows that are all NA for the tumor columns
fold_change_dt <- fold_change_dt[!apply(fold_change_dt[, ..tumor_columns_core], 1, function(row) all(is.na(row)))]

# Set values of 100.199 to 0
fold_change_dt[, (tumor_columns_core) := lapply(.SD, function(x) ifelse(x == 100.199, 0, x)), .SDcols = tumor_columns_core]

# Exclude rows where the sum of all tumor columns is 0
fold_change_dt <- fold_change_dt[rowSums(fold_change_dt[, ..tumor_columns_core], na.rm = TRUE) != 0]


# Save the parquet file
arrow::write_parquet(fold_change_dt, file.path(output_dir, "fold_change_dt.parquet"))

# Calculate log2 fold change for each column individually
for (col in tumor_columns_core) {
  median_value <- medians[col]
  fold_change_dt[, (col) := ifelse(get(col) == 0, NA, log2(get(col) / median_value))]
}

# Step 3: Randomly sample 10,000 cells for each marker, excluding 0 values
set.seed(123) # For reproducibility

sampled_data_list <- lapply(tumor_columns_core, function(col) {
  # Filter the data table for non-NA values in the fold_change column
  # Keep 0s but remove NAs
  sampled_dt <- fold_change_dt[!is.na(get(col)), 
                               .(sample_id, label, fold_change = get(col))]
  
  # Randomly sample 10,000 rows or fewer if there are not enough rows
  sampled_dt <- sampled_dt[sample(.N, min(10000, .N), replace = TRUE)]
  
  # Add the marker column
  sampled_dt[, marker := col]
  
  return(sampled_dt)
})

# Combine the sampled data
sampled_data_dt <- rbindlist(sampled_data_list)


color_final <- c(
  "#607d8b", "#eb6574",  # static requirements
  "#00429d", "#85b7ce", "#b1dfdb", "#ffcab9", "#fd9291", "#e75d6f", 
  "#c52a52", "#93003a", "#c0eade", "#6694c1", "#fb8a8c", "#d5405e",
  "#4e78b5", "#9dced6", "#ffb3a7", "#00429d", "#b81b4a", "#ffdac4",  "#80b1cc", "#325da9"
)

# Define the plot
plot <- ggplot(sampled_data_dt, aes(x = marker, y = fold_change)) +
  #geom_violin(alpha = 0.7, fill = "darkgreen", width = 0.7) +
  geom_jitter(aes(color = fold_change), alpha = 1, width = 0.2, size = 0.1) +
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  scale_color_gradientn(colors = c("black",  "orange", "#c52a52"), 
                        values = scales::rescale(c(-1.5, 0, 1.5)),
                        limits = c(-1.5, 1.5), oob = scales::squish, 
                        breaks = c(-1.5, 0, 1.5), labels = c("-1.5", "0", "1.5")) +
  scale_x_discrete(labels = function(x) gsub("_", "/", x)) +
  scale_y_continuous(limits = c(-1.5, 6.5), breaks = seq(-1.5, 6.5, by = 1)) +  # Set y-axis labels every 0.5 units
  labs(title = "",
       x = "",
       y = "Fold change (log2) from median") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "black", size = 12),
    axis.title = element_text(family = "Helvetica", color = "black", size = 12),
    axis.text = element_text(family = "Helvetica", color = "black", size = 12),
    legend.position = "none",  # Remove the legend
    legend.key.size = unit(1.5, 'lines'),  # Increase legend key size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(size = 0.5, color = "black"),  # Add axis lines
    panel.border = element_blank()  # Remove panel border
  )
plot

# Save the plot as a PNG file
ggsave(
  filename = file.path(output_dir, "sampled_data_plot.png"),
  plot = plot,
  dpi = 1200,
  width = 320,
  height = 100,
  units = "mm"
)

d
###### Figure 5 ################################################################
# Panel c - d - QUICHE network plots generation python notebook -----
# Panel e - DE cells for longitudinal LGG ----


longitudinal_feature_table <- master_feature_table_filtered %>%
  filter(
    site == "UCSF" &
      final_diagnosis_simple %in% c("Astrocytoma", "Oligodendroglioma", "PXA") &
      immunotherapy == "no" &
      paired == "yes" &
      recurrence_status_paired %in% c("primary", "recurrence_1")
  )

longitudinal_feature_table <- longitudinal_feature_table %>%
  mutate(new_column = paste(cell_meta_cluster_final, feature_type, bio_feature_type, feature_variable, sep = "-")) %>%
  select(new_column, everything())

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
immune_cells_with_suffix <- paste0(immune_cells, "_over_all_immune_count_prop")

# Combine the modified immune_cells list with the other_cells list
combined_list <- c(immune_cells_with_suffix, other_cells)

# Print the combined list
combined_list

longitudinal_feature_table <- longitudinal_feature_table[longitudinal_feature_table$feature_variable %in% combined_list, ]

# Step 1: Aggregate data and filter for over_all_immune or over_all_cell
aggregated_data <- longitudinal_feature_table %>%
  group_by(new_column, feature_type, feature_variable) %>%
  filter(recurrence_status_paired %in% c("primary", "recurrence_1")) %>%
  group_by(new_column, feature_type, feature_variable, patient_id, recurrence_status_paired) %>%
  summarize(feature_value = mean(feature_value, na.rm = TRUE), .groups = "drop") %>%
  spread(recurrence_status_paired, feature_value) #%>%

# Reshape the data to long format for plotting
long_data <- aggregated_data %>%
  gather(key = "recurrence_status_paired", value = "feature_value", primary, recurrence_1)


unique(aggregated_data$feature_variable)

# Calculate t-test p-values for each group
p_values <- aggregated_data %>%
  group_by(new_column, feature_type, feature_variable) %>%
  summarize(p_value = t.test(primary, recurrence_1, paired = FALSE)$p.value)

# Merge p-values with long_data for annotation
long_data_with_pvals <- long_data %>%
  left_join(p_values, by = c("new_column", "feature_type", "feature_variable"))


# Create a new column to modify the feature_variable for clean titles
long_data_with_pvals <- long_data_with_pvals %>%
  mutate(
    feature_variable_clean = gsub("_over_all_immune_count_prop", "", feature_variable)
  )

# Step 1: Calculate mean, standard deviation, and p-value in two columns for primary and recurrence
mean_data <- long_data_with_pvals %>%
  group_by(feature_variable_clean) %>%
  summarize(
    mean_primary = mean(feature_value[recurrence_status_paired == "primary"], na.rm = TRUE),
    mean_recurrence = mean(feature_value[recurrence_status_paired == "recurrence_1"], na.rm = TRUE),
    sd_primary = sd(feature_value[recurrence_status_paired == "primary"], na.rm = TRUE),
    sd_recurrence = sd(feature_value[recurrence_status_paired == "recurrence_1"], na.rm = TRUE),
    mean_p_value = mean(p_value, na.rm = TRUE),  # Mean p-value
    .groups = "drop"
  )

# Step 2: Calculate log2FC and propagate errors
log2fc_data <- mean_data %>%
  mutate(
    log2FC = log2(mean_recurrence / mean_primary),  # Calculate log2FC
    se_log2FC = sqrt(
      (sd_primary / (mean_primary * log(2)))^2 +
        (sd_recurrence / (mean_recurrence * log(2)))^2
    )  # Error propagation for log2FC
  ) %>%
  filter(!is.na(log2FC) & !is.infinite(log2FC)) %>%  # Remove invalid values
  filter(!feature_variable_clean %in% c("Unassigned", "Immune","Undetermined", "Immune_unassigned"))  # Exclude unwanted variables

log2fc_plot <- ggplot(log2fc_data, aes(x = log2FC, y = reorder(feature_variable_clean, log2FC))) +
  # Draw the vertical line first to place it at the back
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.25) +
  # Plot points with conditional color for significance
  geom_point(aes(
    size = -log10(mean_p_value),
    color = ifelse(mean_p_value < 0.05, "significant", "non-significant")
  )) +
  # Horizontal error bars with log2FC error (thickness 0.1)
  geom_errorbarh(aes(xmin = log2FC - se_log2FC, xmax = log2FC + se_log2FC), 
                 height = 0.2, color = "black", size = 0.1) +
  # Customize theme: no grid, axis lines with shorter tick marks
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(size = 0.25),  # Thin axis lines
    axis.ticks = element_line(size = 0.25),  # Thin tick marks
    axis.ticks.length = unit(0.5, "mm"),  # Short tick marks (0.5 mm)
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    legend.position = "none"  # Remove legend
  ) +
  # Adjust size scale with smaller max size
  scale_size_continuous(range = c(1, 3)) +  # Smaller max size for points
  # Customize color scale without legend
  scale_color_manual(
    values = c("significant" = "#B81B4A", "non-significant" = "grey")
  )

log2fc_plot

# Save the ggplot as a PDF
ggsave(
  filename = file.path(output_dir, "log2fc_plot_fig_4.pdf"),
  plot = log2fc_plot,  # Use the ggplot object directly
  device = "pdf",
  width = 25 / 25.4,  # Convert mm to inches
  height = 75 / 25.4,  # Convert mm to inches
  units = "in"
)


# Panel f - volcano plot functional marker + cell type DE ----

longitudinal_feature_table <- master_feature_table_filtered %>%
  filter(
    site %in% c("UCSF") &
      final_diagnosis_simple %in% c( "Astrocytoma","Oligodendroglioma") &
      recurrence_status_paired == "recurrence_1" &
      immunotherapy %in% c("no", "yes") 
  )

longitudinal_feature_table <- longitudinal_feature_table %>%
  mutate(new_column = paste(cell_meta_cluster_final, feature_type, bio_feature_type, feature_variable, sep = "-")) %>%
  select(new_column, everything())


# Filter the dataframe for specified conditions
longitudinal_feature_table <- longitudinal_feature_table %>%
  filter(bio_feature_type == "Functional_marker_positivity")

# Display the filtered dataframe
print(longitudinal_feature_table)

# Filter out rows where cell_meta_cluster_final is "Unassigned"
longitudinal_feature_table <- longitudinal_feature_table %>%
  filter(cell_meta_cluster_final != "Unassigned")

# Display the filtered data
print(longitudinal_feature_table)


# Step 1: Aggregate data and filter for over_all_immune or over_all_cell
longitudinal_feature_table <- longitudinal_feature_table %>%
  group_by(new_column, feature_type, feature_variable) %>%
  filter(immunotherapy %in% c("no", "yes")) %>%
  group_by(new_column, feature_type, feature_variable, patient_id, immunotherapy, cell_meta_cluster_final) %>%
  summarize(feature_value = mean(feature_value, na.rm = TRUE), .groups = "drop") 


log2_fc_pval_table <- longitudinal_feature_table %>%
  group_by(cell_meta_cluster_final, feature_variable) %>%
  filter(
    sum(immunotherapy == "yes", na.rm = TRUE) >= 2 &
      sum(immunotherapy == "no", na.rm = TRUE) >= 2
  ) %>%
  summarise(
    mean_yes = mean(feature_value[immunotherapy == "yes"], na.rm = TRUE),
    mean_no = mean(feature_value[immunotherapy == "no"], na.rm = TRUE),
    sdev_yes = sd(feature_value[immunotherapy == "yes"], na.rm = TRUE),
    sdev_no = sd(feature_value[immunotherapy == "no"], na.rm = TRUE),
    log2_fold_change = ifelse(
      mean_yes >= 0.05 & mean_no >= 0.05,
      log2(mean_yes / mean_no),
      NA
    ),
    p_value = ifelse(
      mean_yes >= 0.05 & mean_no >= 0.05,
      t.test(
        feature_value[immunotherapy == "yes"],
        feature_value[immunotherapy == "no"],
        paired = FALSE
      )$p.value,
      NA
    ),
    .groups = "drop"
  ) %>%
  select(
    cell_meta_cluster_final,
    feature_variable,
    log2_fold_change,
    sdev_yes,
    sdev_no,
    p_value
  )

# Display the result
print(log2_fc_pval_table)


# Create a new dataframe excluding rows with more than one "func" in feature_variable from log2_fc_pval_table
filtered_log2_fc_pval_table <- log2_fc_pval_table %>%
  filter(
    str_count(feature_variable, "func") <= 1
  )

# Display the new dataframe
print(filtered_log2_fc_pval_table)

library(ggplot2)

# Adjust the dataframe for significance mapping
significant_log2_fc_pval_table <- filtered_log2_fc_pval_table %>%
  mutate(
    significance = ifelse(
      abs(log2_fold_change) > 1 & p_value < 0.0599,
      "Significant",
      "Not Significant"
    )
  )

# Aggregate data to remove duplicates
aggregated_table <- significant_log2_fc_pval_table %>%
  group_by(log2_fold_change, p_value, significance) %>%
  summarise(count = n(), .groups = "drop")

volcano_plot <- ggplot(aggregated_table, aes(x = log2_fold_change, y = -log10(p_value))) +
  geom_point(aes(
    color = significance,
    stroke = 0,
    alpha = ifelse(significance == "Not Significant", 0.7, 1)
  ), size = 2) + # Set point size to 1
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey", size = 0.2) + # Add vertical lines at ±1
  geom_hline(yintercept = -log10(0.0599), linetype = "dashed", color = "grey", size = 0.2) + # Add horizontal line at p-value = 0.1
  scale_color_manual(
    values = c("Significant" = "#B81B4A", "Not Significant" = "grey"),
    guide = "none"
  ) +
  scale_alpha_identity() + # Ensures the alpha values are taken directly from the data
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.title = element_blank(), # Remove axis titles
    axis.text = element_text(size = 7, color = "black"), # Set axis labels to Arial, size 7
    axis.ticks = element_line(size = 0.25), # Set tick mark thickness to 0.25
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

# Display the plot
print(volcano_plot)

# Export the plot as a PDF
ggsave(
  filename = file.path(output_dir, "log2fc_plot_fig_4_GBM_immuno.pdf"),
  plot = volcano_plot,
  device = "pdf",
  dpi = 600,
  width = 75 / 25.4, # Convert mm to inches
  height = 60 / 25.4 # Convert mm to inches
)



###### Figure 6 ################################################################
# Panel a - glycan class bar plot ----


MALDI_data <- MALDI_data[MALDI_data$site == "Stanford", ]

# Identify the exact positions of columns starting with "H2N2F1" and ending with "H6N7F1"
start_col <- grep("^H2N2F1$", colnames(MALDI_data))
end_col <- grep("^H6N7F1$", colnames(MALDI_data))

# Check if the start and end columns exist and the range is valid
if (length(start_col) == 1 && length(end_col) == 1 && start_col < end_col) {
  # Subset column names before, in the range, and after
  cols_before <- colnames(MALDI_data)[1:(start_col - 1)]
  cols_in_range <- colnames(MALDI_data)[start_col:end_col]
  cols_after <- colnames(MALDI_data)[(end_col + 1):ncol(MALDI_data)]
  
  # Remove columns with any NA values within the range
  valid_cols_in_range <- cols_in_range[colSums(is.na(MALDI_data[, cols_in_range])) == 0]
  
  # Combine all columns: before, valid in range, and after
  final_cols <- c(cols_before, valid_cols_in_range, cols_after)
  
  # Create the updated dataframe with the filtered columns
  MALDI_data <- MALDI_data[, final_cols]
} else {
  stop("Could not find the specified range of columns or the range is invalid.")
}

# Display the first few rows of the filtered dataframe
head(MALDI_data)

classifications <- colnames(classes)[3:12]

# Assuming 'classes' is your dataframe
glycan_counts <- data.frame(Column = character(), Num_True = integer(), Max_Possible_True = integer(), stringsAsFactors = FALSE)

# Iterate over columns 3 to 12
for (i in 3:12) {
  column_name <- colnames(classes)[i]
  
  # Filter rows where the current column equals 1
  filtered_df <- classes %>% filter(.[[i]] == 1)
  
  # Get unique values in the 'composition' column
  unique_compositions <- unique(filtered_df$composition)
  
  # Check if 'Stanford_MALDI' column names match the unique composition values
  matched_columns <- colnames(MALDI_data) %in% unique_compositions
  
  # Count the number of TRUE values
  num_true <- sum(matched_columns)
  
  # Get the maximum possible number of TRUE values
  max_possible_true <- length(unique_compositions)
  
  # Add the results to the dataframe
  glycan_counts <- rbind(glycan_counts, data.frame(Column = column_name, Num_True = num_true, Max_Possible_True = max_possible_true, stringsAsFactors = FALSE))
}

# Print the results dataframe
print(glycan_counts)

# Create a mapping of old names to new names with capitalized first letters and spaces
name_mapping <- c(
  "fucosylated" = "Fucosylated",
  "sialylated" = "Sialylated",
  "highMannose" = "High Mannose",
  "hybrid" = "Hybrid",
  "paucimannose" = "Paucimannose",
  "agalactosylated" = "Agalactosylated",
  "biantennary" = "Biantennary",
  "triantennary" = "Triantennary",
  "tetraantennary" = "Tetraantennary",
  "polylacnac" = "Polylacnac"
)

# Apply the name mapping to the classification column
glycan_counts$classification <- factor(glycan_counts$Column, levels = names(name_mapping))
glycan_counts$classification <- plyr::revalue(glycan_counts$classification, name_mapping)

glycan_counts$percent <- (glycan_counts$Num_True / 69)*100


# Plot the bar graph
p <- ggplot(glycan_counts, aes(x = reorder(classification, -percent), y = percent)) +
  geom_bar(stat = "identity", fill = "#85B7CE") +
  geom_text(aes(label = sprintf("%.0f%%", percent)), vjust = -0.5, size = 2, angle = 0)+
  theme_minimal() +
  labs(
    x = NULL,
    y = "Glycan (% of all)"
  ) +
  theme(
    axis.title.x = element_text(size = 6, family = "Helvetica", color = "black"),
    axis.title.y = element_text(size = 6, family = "Helvetica", color = "black"),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, family = "Helvetica", color = "black"),
    axis.text.y = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines inside
    panel.border = element_blank(),  # Remove outer border
    axis.line = element_line(color = "black", linewidth = 0.25)  # Add axis lines with specified linewidth
  ) +
  scale_y_continuous(limits = c(0, 85))  # Adjust y-axis to avoid cutting off the bars

print(p)

# Export as PNG file with specified dimensions and resolution
ggsave(
  filename = file.path(output_dir, "glycan_bar_graph.pdf"),
  plot = p,
  dpi = 1200,
  width = 60,
  height = 40,
  units = "mm"
)


# Panel b - glycan class heatmap ----

# Define the name mapping
name_mapping <- c(
  "fucosylated" = "Fucosylated",
  "sialylated" = "Sialylated",
  "highMannose" = "High Mannose",
  "hybrid" = "Hybrid",
  "paucimannose" = "Paucimannose",
  "agalactosylated" = "Agalactosylated",
  "biantennary" = "Biantennary",
  "triantennary" = "Triantennary",
  "tetraantennary" = "Tetraantennary",
  "polylacnac" = "Polylacnac"
)

# Filter the dataframe based on the given conditions
filtered_df <- master_feature_table_filtered %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"),
         grepl("glycans", bio_feature_type),
         feature_variable %in% names(name_mapping))

# Summarize the data by feature_variable and who_grade
summary_df <- filtered_df %>%
  group_by(feature_variable, who_grade) %>%
  summarize(mean_value = mean(feature_value, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = who_grade, values_from = mean_value, names_prefix = "grade_")


# Convert summary_df to a matrix and z-score row-wise
heatmap_matrix <- as.matrix(summary_df %>% column_to_rownames(var = "feature_variable"))
z_scored_matrix <- t(apply(heatmap_matrix, 1, scale))

# Rename columns to 1, 2, 3 and reorder to 2, 3, 4
colnames(z_scored_matrix) <- c("2", "3", "4")
z_scored_matrix <- z_scored_matrix[, c("2", "3", "4")]

# Rename rows based on the name_mapping
rownames(z_scored_matrix) <- name_mapping[rownames(z_scored_matrix)]

# Create the transposed heatmap with the specified settings
ht <- Heatmap(t(z_scored_matrix),  # Transpose the matrix
              name = "Z-Score",
              row_names_gp = gpar(fontsize = 10, fontfamily = "Helvetica", col = "black"),
              column_names_gp = gpar(fontsize = 10, fontfamily = "Helvetica", col = "black"),
              cluster_rows = FALSE,  # Disable row clustering since columns are now transposed to rows
              cluster_columns = TRUE,  # Enable column clustering (original rows)
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_names_rot = 90,  # Rotate column names by 90 degrees
              rect_gp = gpar(col = "black", lwd = 0.25),  # Add black outline around each box
              show_heatmap_legend = FALSE,  # Remove legend
              row_names_side = "left",  # Place row names on the left
              col = colorRamp2(c(-1.2, 0, 1.2), c("blue", "white", "#b81b4a")))  # Change color scheme

# Draw the transposed heatmap
draw(ht)

# Export as PNG
png(
  filename = file.path(output_dir, "heatmap_glycan_classes.png"),
  width = 80,
  height = 70,
  units = "mm",
  res = 600
)
draw(ht)
dev.off()




# Panel c - individual glycan box plots ----

filtered_df <- master_feature_table_filtered %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"),
         feature_variable %in% c("H3N5F1", "H5N4S1","H5N2"))

wide_df <- filtered_df %>%
  pivot_wider(names_from = feature_variable,
              values_from = feature_value,
              id_cols = c(sample_id, who_grade))

# Assuming `wide_df` is your dataframe
unique_features <- c("H5N4S1", "H3N5F1","H5N2")

# Loop through each feature, create the plot, and save it with the specified settings
for (feature in unique_features) {
  wide_df <- wide_df %>%
    mutate(
      point_size = 0.05,  # Fixed point size for all points
      point_color = "black"  # Fixed point color for all points
    )
  
  p <- ggplot(wide_df, aes(x = who_grade, y = .data[[feature]])) +
    geom_boxplot(fill = "#618fbf", color = "black", outlier.shape = NA, size = 0.25) +
    geom_jitter(aes(color = point_color, size = point_size), width = 0.2, alpha = 1) +
    scale_color_identity() +
    scale_size_identity() +
    labs(title = feature,  # Set the figure title to the feature name
         x = NULL,
         y = "Relative Intensity") +  # Set y-axis title to "Relative Intensity"
    theme_minimal() +
    theme(text = element_text(family = "Helvetica", size = 6, color = "black"),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, family = "Helvetica", color = "black"),
          axis.text.y = element_text(size = 6, family = "Helvetica", color = "black"),
          axis.title.y = element_text(size = 6, family = "Helvetica", color = "black"),
          plot.title = element_blank(),
          panel.grid = element_blank(),     # Remove grid lines
          axis.line = element_line(color = "black", linewidth = 0.25))  # Add axis lines
  
  # Save the plot
  ggsave(
    filename = file.path(output_dir, paste0("boxplot_", feature, ".png")),
    plot = p,
    dpi = 600,
    height = 35,
    width = 35,
    units = "mm"
  )
  print(p)
}


# Panel d - glycan + gene DE heatmap ----

genes_to_filter <- c("FUT1", "FUT2", "FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT8", "FUT9",
                     "ST6GAL1", "ST6GAL2", "ST3GAL4", "ST3GAL6", "MGAT1", "MGAT2", "MGAT4a", "MGAT4b", "MGAT5",
                     "B3GNT2", "MAN1B1", "MAN1A1", "MAN1A2", "MAN1C1", "B3GNT3", "B3GNT4", "B3GNT7", "B3GNT8", "B3GNT9",
                     "B4GALT1", "B4GALT2", "B4GALT3", "B4GALT4", "B3GALT1", "B3GALT2", "B3GALT5", "GCNT2", "GCNT7",
                     "MGAT3", "MAN2A1", "MAN2A2")

# Define the name mapping
name_mapping <- c(
  "fucosylated" = "Fucosylated",
  "sialylated" = "Sialylated",
  "highMannose" = "High Mannose",
  "hybrid" = "Hybrid",
  "paucimannose" = "Paucimannose",
  "agalactosylated" = "Agalactosylated",
  "biantennary" = "Biantennary",
  "triantennary" = "Triantennary",
  "tetraantennary" = "Tetraantennary",
  "polylacnac" = "Polylacnac"
)

# Get unique feature_variable values based on the given conditions
glycan_cols <- master_feature_table_filtered %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"),
         grepl("glycans", bio_feature_type),
         !feature_variable %in% names(name_mapping),
         feature_variable != "sum_intensity") %>%
  pull(feature_variable) %>%
  unique()
glycan_cols


filtered_master_feature_table <- master_feature_table_filtered %>%
  filter(
    grepl("Stanford", sample_id) &
      (feature_variable %in% glycan_cols |feature_variable %in% genes_to_filter ))


filtered_master_feature_table <- filtered_master_feature_table %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"))

master_feature_pivot <- filtered_master_feature_table %>%
  dplyr::select(patient_id, feature_variable, feature_value, who_grade, tumor_region) %>%
  pivot_wider(names_from = feature_variable, values_from = feature_value, values_fn = mean) %>%
  dplyr::select(patient_id, everything())



# Merge master_feature_pivot with metadata using sample_id
merged_data <- master_feature_pivot %>%
  left_join(metadata, by = c("patient_id", "who_grade"))

data_glycan <- merged_data 



# Perform differential expression analysis for glycan data
results_glycan <- list()
for (glycan in glycan_cols) {
  formula <- as.formula(paste(glycan, "~ who_grade"))
  model <- aov(formula, data = data_glycan)
  p_value <- summary(model)[[1]]["Pr(>F)"][1]
  results_glycan[[glycan]] <- p_value
}


# Initialize a list to store the p-values
# Initialize a vector to store significant glycans
significant_glycans <- c()

# Loop through results_glycan to find significant glycans
for (glycan in names(results_glycan)) {
  p_value <- results_glycan[[glycan]]["Pr(>F)"][[1]]
  
  # Check if p-value is less than 0.05
  if (!is.na(p_value) && p_value < 0.1) {
    significant_glycans <- c(significant_glycans, glycan)
  }
}

# Print the significant glycans
print(significant_glycans)


significant_glycans <- significant_glycans[!is.na(significant_glycans)]

# Subset data for significant glycans
data_significant_glycan <- data_glycan %>% dplyr::select(all_of(significant_glycans), who_grade)

# Prepare data for heatmap and summarize duplicates by taking the mean
data_long_glycan <- data_significant_glycan %>%
  pivot_longer(-who_grade, names_to = "feature", values_to = "expression") %>%
  group_by(feature, who_grade) %>%
  summarise(expression = mean(expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = who_grade, values_from = expression)


# Merge master_feature_pivot with metadata using sample_id
merged_data <- master_feature_pivot %>%
  left_join(metadata, by = c("patient_id", "who_grade"))



# Assuming master_feature_pivot is available
selected_columns <- c(genes_to_filter)

# Find columns that are present in merged_data and also in columns_to_exclude
selected_columns <- intersect(selected_columns, colnames(merged_data))

data_features <- merged_data 


# Perform differential expression analysis for master_feature_pivot data
results_features <- list()
for (feature in selected_columns) {
  formula <- as.formula(paste(feature, "~ who_grade"))
  model <- aov(formula, data = data_features)
  p_value <- summary(model)[[1]]["Pr(>F)"][1]
  results_features[[feature]] <- p_value
}

# Extract p-values from results_features
p_values <- sapply(results_features, function(x) x["who_grade", "Pr(>F)"])

# Filter significant features based on p-values
significant_features <- names(p_values)[p_values < 0.15]

# Display the significant features
significant_features


significant_features <- significant_features[!is.na(significant_features)]

# Subset data for significant features
data_significant_features <- data_features %>% dplyr::select(all_of(significant_features), who_grade)

# Prepare data for heatmap and summarize duplicates by taking the mean
data_long_features <- data_significant_features %>%
  pivot_longer(-who_grade, names_to = "feature", values_to = "expression") %>%
  group_by(feature, who_grade) %>%
  summarise(expression = mean(expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = who_grade, values_from = expression)

# Combine both datasets
data_long_glycan$feature_info <- paste0(data_long_glycan$feature, "_glycan")

# Modify data_long_features to handle features ending in _prop
data_long_features$feature_info <- ifelse(grepl("_prop$", data_long_features$feature), 
                                          paste0(data_long_features$feature, "_cells"), 
                                          paste0(data_long_features$feature, "_gene"))

# Combine the datasets
combined_data_long <- bind_rows(data_long_glycan, data_long_features)



# Create a new dataframe with row-wise z-scores
combined_data_long_z <- combined_data_long %>%
  mutate(across(`2`:`4`, as.numeric)) %>%  # Ensure the columns are numeric
  rowwise() %>%
  mutate(across(`2`:`4`, ~ (. - mean(c_across(`2`:`4`))) / sd(c_across(`2`:`4`)))) %>%
  ungroup()
combined_data_long
combined_data_long_z


# Save combined data as CSV
write.csv(combined_data_long, file.path(output_dir, "combined_data_long_heatmap_DE_glycans.csv"), row.names = FALSE)

# Read input data
combined_data_long_z <- read.csv(file.path(base_dir, "gbm_mean_by_stage_table_norm_gly_enz_new.csv")) # downloadable at https://bruce.parkerici.org/pages/raw-data-access.html


colnames(combined_data_long_z) <- c("feature", "2", "3", "4", "k_cluster")




combined_data_long_z$feature_info <- ifelse(grepl("^H", combined_data_long_z$feature), 
                                            paste0(combined_data_long_z$feature, "_glycan"), 
                                            paste0(combined_data_long_z$feature, "_gene"))



combined_data_long_z

# Create annotation data frame
annotation <- data.frame(
  Type = ifelse(grepl("_glycan$", combined_data_long_z$feature_info), 
                "Glycan", 
                ifelse(grepl("_cells$", combined_data_long_z$feature_info), 
                       "Cells", 
                       "Gene")),
  row.names = combined_data_long_z$feature
)

# Set rownames for heatmap
heatmap_data <- as.data.frame(combined_data_long_z %>% dplyr::select(-feature)%>% dplyr::select(-feature_info))
rownames(heatmap_data) <- combined_data_long_z$feature

# Ensure the data is numeric
heatmap_data[] <- lapply(heatmap_data, as.numeric)

# Convert to numeric matrix for heatmap
heatmap_matrix <- as.matrix(heatmap_data)

library(ComplexHeatmap)
library(dplyr)

# # Set the number of clusters
# K <- 6

# Perform row Z-score normalization
heatmap_matrix_scaled <- t((t(heatmap_matrix)))

# Perform hierarchical clustering
row_hclust <- hclust(dist(heatmap_matrix_scaled, method = "euclidean"), method = "average")

row_hclust


# Ensure the column order
column_order <- c("2", "3", "4")
heatmap_matrix_scaled <- heatmap_matrix_scaled[, column_order]


# Assuming heatmap_matrix_scaled is a matrix with row names
row_names <- rownames(heatmap_matrix_scaled)

# Remove "_func_over_all_tumor_prop" and "_over_all_immune_count_prop"
row_names <- str_replace_all(row_names, "_func_over_all_tumor_prop", "")
row_names <- str_replace_all(row_names, "_over_all_immune_count_prop", "")

# Replace "_" before "Na" with " + "
row_names <- str_replace_all(row_names, "_(?=Na)", " + ")

# Change "plus" to "+"
row_names <- str_replace_all(row_names, "plus", "+")

# Replace any remaining "_" with a space
row_names <- str_replace_all(row_names, "_", " ")

# Set the new row names to the matrix
rownames(heatmap_matrix_scaled) <- row_names

# Print the new row names to verify
print(rownames(heatmap_matrix_scaled))


# Set annotation colors
annotation_colors <- list(Type = c("Glycan" = "#b1dfdb", "Cells" = "black", "Gene" = "#00429d"))

# Update annotation based on transposed structure (Type will now be a row annotation)
annotation <- data.frame(Type = annotation$Type)

# Define custom color function
col_fun <- colorRamp2(c(-1.2, 0, 1.2), c("#00429D", "white", "#B81B4A"))
set.seed(2)

# Transpose row split values to apply along columns instead
column_split_values <- heatmap_data$k_cluster  # This is now used to split columns

bottom_annotation <- columnAnnotation(
  Type = anno_simple(annotation$Type, 
                     col = annotation_colors$Type, 
                     border = TRUE,  # Add border to each annotation box
                     gp = gpar(col = "black", lwd = 0.25)),
  show_annotation_name = FALSE  # Hide the annotation title
)

# Create the heatmap with column splitting based on k_cluster
heatmap_result <- Heatmap(t(heatmap_matrix_scaled),  # Transpose the matrix
                          name = "Z scored",
                          clustering_distance_columns = "euclidean",  # Cluster columns (original rows)
                          cluster_columns = TRUE,  # Enable column clustering
                          cluster_rows = FALSE,  # Disable row clustering
                          show_column_names = TRUE,
                          show_row_names = TRUE,
                          border = TRUE,
                          column_split = column_split_values,  # Use k_cluster for column splitting
                          column_title = NULL,  # Remove k_cluster column labels
                          row_title = NULL,  # Remove k_cluster row labels
                          show_column_dend = TRUE,  # Show the column dendrogram
                          bottom_annotation = bottom_annotation,  # Apply the bottom annotation with outlines
                          col = col_fun,
                          rect_gp = gpar(col = "black", lwd = 0.35),  # Add black outline around each heatmap box
                          column_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica", color = "black"),
                          row_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica", color = "black"),
                          row_names_side = "left")  # Place row names on the left

# Draw the transposed heatmap without any legends
draw(heatmap_result, 
     show_heatmap_legend = FALSE,  # Disable the heatmap legend
     show_annotation_legend = FALSE)  # Disable the annotation legend



# Export as PNG
png(
  filename = file.path(output_dir, "heatmap_glycan_DE.png"),
  width = 180,
  height = 90,
  units = "mm",
  res = 600
)
draw(
  heatmap_result,
  show_heatmap_legend = FALSE,  # Disable the heatmap legend in the exported image
  show_annotation_legend = FALSE  # Disable the annotation legend in the exported image
)
dev.off()



# Panel e - Cell - glycan - RNA correlation illustration ----


# Load necessary libraries
library(dplyr)
library(tidygraph)
library(igraph)

# Filter the significant data for positive correlations and specific glycans
filtered_data <- significant_data %>%
  filter(Correlation > 0, 
         Glycan %in% c("sialylated", "fucosylated", "triantennary", "agalactosylated"))

normalize <- function(x) {
  x / max(x, na.rm = TRUE)
}

merged_all_top_go_df <- merged_all_top_go_df %>%
  mutate(
    Description = case_when(
      # Neuronal development and neurotransmitter signaling
      Original_Descriptions %in% c(
        "modulation of chemical synaptic transmission_fucosylated / regulation of trans-synaptic signaling_fucosylated",
        "synapse organization_fucosylated",
        "regulation of neurotransmitter levels_fucosylated / neurotransmitter secretion_fucosylated / signal release from synapse_fucosylated / neurotransmitter transport_fucosylated",
        "signal release_fucosylated",
        "learning or memory_fucosylated / cognition_fucosylated",
        "regulation of synaptic plasticity_agalactosylated",
        "regulation of neurotransmitter levels_agalactosylated / neurotransmitter secretion_agalactosylated / signal release from synapse_agalactosylated / neurotransmitter transport_agalactosylated",
        "glutamate receptor signaling pathway_agalactosylated",
        "social behavior_agalactosylated / biological process involved in intraspecies interaction between organisms_agalactosylated",
        "signal release_agalactosylated",
        "modulation of chemical synaptic transmission_agalactosylated / regulation of trans-synaptic signaling_agalactosylated",
        "neuron cell-cell adhesion_agalactosylated"
      ) ~ "Neuronal development and neurotransmitter signaling",
      
      # Vesicle mediated transport
      Original_Descriptions %in% c(
        "synaptic vesicle cycle_fucosylated / vesicle-mediated transport in synapse_fucosylated",
        "regulated exocytosis_fucosylated / exocytosis_fucosylated",
        "synaptic vesicle exocytosis_fucosylated",
        "synaptic vesicle cycle_agalactosylated / vesicle-mediated transport in synapse_agalactosylated",
        "regulated exocytosis_agalactosylated"
      ) ~ "Vesicle mediated transport",
      
      # Antigen presentation (MHC class II)
      Original_Descriptions %in% c(
        "antigen processing and presentation of peptide antigen_sialylated / antigen processing and presentation_sialylated / antigen processing and presentation of exogenous peptide antigen_sialylated / antigen processing and presentation of peptide antigen via MHC class II_sialylated / antigen processing and presentation of peptide or polysaccharide antigen via MHC class II_sialylated / antigen processing and presentation of exogenous antigen_sialylated / antigen processing and presentation of exogenous peptide antigen via MHC class II_sialylated",
        "antigen processing and presentation of peptide antigen via MHC class II_triantennary / antigen processing and presentation of peptide or polysaccharide antigen via MHC class II_triantennary / antigen processing and presentation of peptide antigen_triantennary / MHC class II protein complex assembly_triantennary / peptide antigen assembly with MHC class II protein complex_triantennary / antigen processing and presentation of exogenous peptide antigen via MHC class II_triantennary / antigen processing and presentation_triantennary / MHC protein complex assembly_triantennary / peptide antigen assembly with MHC protein complex_triantennary / antigen processing and presentation of exogenous peptide antigen_triantennary"
      ) ~ "Antigen presentation (MHC class II)",
      
      # Leukocyte mediated immune activation
      Original_Descriptions %in% c(
        "leukocyte mediated immunity_sialylated / lymphocyte mediated immunity_sialylated",
        "regulation of immune effector process_sialylated",
        "leukocyte mediated immunity_triantennary / lymphocyte mediated immunity_triantennary / adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains_triantennary / immunoglobulin mediated immune response_triantennary / B cell mediated immunity_triantennary"
      ) ~ "Leukocyte mediated immune activation",
      
      # Negative regulation of immune processes
      Original_Descriptions == "negative regulation of immune system process_sialylated" ~ "Negative regulation of immune processes",
      
      # Leukocyte mediated cytotoxicity
      Original_Descriptions %in% c(
        "cell killing_sialylated / leukocyte mediated cytotoxicity_sialylated"
      ) ~ "Leukocyte mediated cytotoxicity",
      
      # Cell adhesion and wound healing
      Original_Descriptions %in% c(
        "positive regulation of cell adhesion_sialylated",
        "wound healing_sialylated"
      ) ~ "Cell adhesion and wound healing",
      
      # Default case: keep the original description
      TRUE ~ Original_Descriptions
    )
  )


# Normalize gene_count weights in merged_all_top_go_df
merged_all_top_go_df <- merged_all_top_go_df %>%
  mutate(Gratio = normalize(Gratio))


# Merge the filtered data with descriptions and normalized gene counts
plot_data <- filtered_data %>%
  inner_join(merged_all_top_go_df %>% select(glycan = glycan, Description, Gratio), 
             by = c("Glycan" = "glycan"))

plot_data <- plot_data %>%
  mutate(
    Feature = sub("_.*", "", Feature),                # Extract text before the first underscore
    Glycan = tools::toTitleCase(as.character(Glycan)) # Convert to character and capitalize
  )

# View the updated `plot_data`
print(plot_data)

# Prepare nodes with appropriate layers
nodes <- tibble(
  name = unique(c(plot_data$Feature, plot_data$Glycan, plot_data$Description)),
  layer = case_when(
    name %in% plot_data$Feature ~ 1,
    name %in% plot_data$Glycan ~ 2,
    TRUE ~ 3
  ),
  type = case_when(
    name %in% plot_data$Glycan ~ "Glycan",
    TRUE ~ "Other"
  )
)

# Create edges with normalized weights only for Glycan -> Description
edges <- bind_rows(
  # Feature -> Glycan edges (Correlation weight unchanged)
  plot_data %>% 
    select(from = Feature, to = Glycan, weight = Correlation) %>%
    mutate(type = "Feature-to-Glycan"),
  
  # Glycan -> Description edges (Normalized gene_count)
  plot_data %>% 
    select(from = Glycan, to = Description, weight = Gratio) %>%
    mutate(type = "Glycan-to-Description")
) %>%
  filter(from %in% nodes$name, to %in% nodes$name)

# Debug: Print unique weights to verify scaling
print(unique(edges$weight))

glycan_colors <- c(
  "Sialylated" = "#7e57c2",    # Purple shade
  "Fucosylated" = "#ffd54f",   # Yellow shade
  "Triantennary" = "#4fc3f7",  # Light blue shade
  "Agalactosylated" = "#ff69b4" # Brown shade
)

# Map edge colors based on Glycan nodes they connect to
edges <- edges %>%
  mutate(
    edge_color = case_when(
      to %in% names(glycan_colors) ~ glycan_colors[to],
      from %in% names(glycan_colors) ~ glycan_colors[from],
      TRUE ~ "black"
    )
  )

# Create the graph object
graph <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)




# Adjust node positions to spread out vertically within each layer
layout <- nodes %>%
  group_by(layer) %>%
  mutate(
    y = seq(-1, 1, length.out = n()),  # Spread nodes evenly within each layer
    x = layer
  )

layout <- layout %>%
  mutate(
    x = case_when(
      name == "Endothelial" ~ 1,  # Adjust to be in the Glycan layer
      name == "Fucosylated" ~ 2,  # Place in the pathway layer
      TRUE ~ x
    ),
    y = case_when(
      # Tumor, Neurons, Immune at the top
      name == "Tumor" ~ 0.7,
      name == "Neurons" ~ 0.6,
      name == "Immune" ~ 0.5,
      
      # Agalactosylated, Fucosylated, Sialylated, and Triantennary in order
      name == "Agalactosylated" ~ 0.85,
      name == "Fucosylated" ~ 0.65,
      name == "Sialylated" ~ 0.45,
      name == "Triantennary" ~ 0.25,
      
      # Remaining nodes sorted and spread equally from 1 to 0.1
      name == "Neuronal development and neurotransmitter signaling" ~ 1.0,
      name == "Vesicle mediated transport" ~ 0.85,
      name == "Negative regulation of immune processes" ~ 0.7,
      name == "Leukocyte mediated cytotoxicity" ~ 0.55,
      name == "Cell adhesion and wound healing" ~ 0.4,
      name == "Antigen presentation (MHC class II)" ~ 0.25,
      name == "Leukocyte mediated immune activation" ~ 0.1,
      
      # Default case
      TRUE ~ y
    )
  )



# Export as PNG with transparent background
png(
  filename = file.path(output_dir, "Feature_Glycan_Description_Connections.png"),
  bg = "transparent",
  width = 150,
  height = 125,
  units = "mm",
  res = 600
)

ggraph(graph, layout = "manual", x = layout$x, y = layout$y) +
  geom_edge_link(aes(width = weight, color = edge_color, 
                     linetype = ifelse(type == "Glycan-to-Description", "dashed", "solid")),
                 alpha = 1, lineend = "round") +  scale_edge_width(range = c(0.2, 3)) +
  scale_edge_color_identity() +
  geom_node_point(aes(color = ifelse(type == "Glycan", name, "black")), size = 10) +
  geom_node_point(color = "white", size = 8) +  # Inner white center for donut effect
  scale_color_manual(values = c(glycan_colors, black = "black")) +
  # Set node text to blank
  #geom_node_text(aes(label = ""), hjust = 0, vjust = -2, size = 3, nudge_x = 0.1) +
  #geom_node_text(aes(label = name), hjust = 0, vjust = -2, size = 3, nudge_x = 0.1) +
  theme_graph() +
  #labs(title = "Feature-Glycan-Description Connections") +
  theme(
    plot.margin = margin(10, 10, 10, 10),  # Increase right margin to avoid cut-off
    legend.position = "none"
  ) +
  # Expand plot limits to ensure no overlap or cut-off
  expand_limits(x = 5)

dev.off()

# Panel f - glycan to gene direct correlation heatmap ----


# Define the list of glycan names to loop through
glycan_names_to_extract <- c("fucosylated", "sialylated", "triantennary", "agalactosylated")#,"hybrid","biantennary","highMannose")

# Initialize an empty dataframe to store results
all_glycan_genes_correlation_df <- data.frame()

# Loop through each glycan name
for (glycan_name in glycan_names_to_extract) {
  
  # Find the index of the current glycan name
  glycan_index <- which(sapply(go_results, function(x) as.character(attr(x, "glycan_name"))) == glycan_name)
  
  # Check if the index was found and extract the result
  if (length(glycan_index) > 0) {
    glycan_result <- go_results[[glycan_index]]
    
    # Convert the result into a dataframe if it has data
    if (nrow(glycan_result@result) > 0) {
      # Create a temporary dataframe to process the data
      temp_df <- data.frame(
        Pathway_Name = glycan_result@result$Description,
        Genes_Found = sapply(glycan_result@result$geneID, function(x) paste(unlist(strsplit(x, "/")), collapse = ", ")),
        p_value = glycan_result@result$pvalue,
        Count = glycan_result@result$Count
      )
      
      # Select the top 5 pathways based on Count
      top_5_paths <- temp_df %>%
        arrange(desc(Count)) %>%
        head(5)
      
      # Extract individual genes from Genes_Found column
      unique_genes <- unique(unlist(strsplit(paste(top_5_paths$Genes_Found, collapse = ", "), ", ")))
      
      # Filter correlation matrix for the current glycan and genes
      if (glycan_name %in% colnames(correlation_matrix)) {
        glycan_correlation <- correlation_matrix[unique_genes, glycan_name, drop = FALSE]
        
        # Create a dataframe with Gene, Correlation, and Glycan_Name
        glycan_gene_correlation_df <- data.frame(
          Gene = rownames(glycan_correlation),
          Correlation = as.numeric(glycan_correlation[, 1]),
          Glycan_Name = glycan_name
        )
        
        # Append the results to the main dataframe
        all_glycan_genes_correlation_df <- bind_rows(all_glycan_genes_correlation_df, glycan_gene_correlation_df)
        
      } else {
        cat(glycan_name, "not found in correlation matrix.\n")
      }
      
    } else {
      cat("No enriched pathways found for", glycan_name, "\n")
    }
  } else {
    cat(glycan_name, "not found among glycan names.\n")
  }
}

# Display the final dataframe with genes, correlations, and glycan names
print(all_glycan_genes_correlation_df)


# Initialize an empty list to store gene-glycan correlation data
unique_genes <- unique(all_glycan_genes_correlation_df$Gene)
glycan_columns <- glycan_names_to_extract

# Create a new dataframe to hold each gene and its correlation to all glycans
gene_glycan_correlation_df <- data.frame(matrix(ncol = length(glycan_columns), nrow = length(unique_genes)))
colnames(gene_glycan_correlation_df) <- glycan_columns
rownames(gene_glycan_correlation_df) <- unique_genes

# Fill the dataframe by looping through each unique gene
for (gene in unique_genes) {
  if (gene %in% rownames(correlation_matrix)) {
    # Use only specified glycan columns to extract the correlation for the gene
    gene_glycan_correlation_df[gene, glycan_columns] <- correlation_matrix[gene, glycan_columns]
  } else {
    # Fill with NA if the gene is not found in correlation_matrix
    gene_glycan_correlation_df[gene, glycan_columns] <- NA
  }
}

# Display the final dataframe
print(gene_glycan_correlation_df)


# Sort all_glycan_genes_correlation_df by Correlation in descending order and get the top 25 genes
top_25_genes <- all_glycan_genes_correlation_df %>%
  arrange(desc(Correlation)) %>%
  head(25) %>%
  pull(Gene)

# Filter gene_glycan_correlation_filtered_df to include only the top 25 genes
top_25_gene_glycan_correlation_df <- gene_glycan_correlation_df[rownames(gene_glycan_correlation_df) %in% top_25_genes, ]

# Display the filtered dataframe
print(top_25_gene_glycan_correlation_df)


# Load ComplexHeatmap package
library(ComplexHeatmap)
library(circlize)  # For colorRamp2

# Define the matrix for the heatmap (correlation values)
heatmap_matrix <- as.matrix(top_25_gene_glycan_correlation_df)


# Define color scale with specific colors for desired ranges
col_fun <- colorRamp2(
  c(min(heatmap_matrix, na.rm = TRUE), 0, 0, 0.4, 0.6), 
  c("white", "white", "white", "white", "#B81B4A")
)

# Plot the heatmap with ComplexHeatmap
ht <- Heatmap(
  heatmap_matrix,
  name = "Correlation",                # Name for the heatmap legend
  col = col_fun,                       # Apply color scale
  cluster_rows = TRUE,                 # Enable row clustering
  cluster_columns = TRUE,              # Enable column clustering
  show_row_names = TRUE,               # Show row names (gene names)
  show_column_names = TRUE,            # Show column names (glycan names)
  border = TRUE,                       # Enable borders
  border_gp = gpar(col = "black", lwd = 0.35),  # Set border color and thickness
  rect_gp = gpar(col = "black", lwd = 0.35),    # Black outline around each heatmap box
  row_dend_width = unit(5, "mm"),      # Set row dendrogram width to narrow it
  column_dend_height = unit(5, "mm"),  # Set column dendrogram height to narrow it
  #row_names_gp = gpar(fontfamily = "Arial", size = 10),    # Set row names font to Arial
  show_heatmap_legend = FALSE         # Show the heatmap legend
)

# Draw the heatmap
draw(ht)

# Save the heatmap with specified dimensions and resolution
png(
  filename = file.path(output_dir, "heatmap_output_panel_e_mod.png"),
  width = 60,
  height = 155,
  units = "mm",
  res = 600
)
draw(ht)
dev.off()




# Panel g - h - gene to glycan correlation scatter plot ----

# Define the features to use in the plot
feature_x <- "sialylated"  # Replace with desired feature for x-axis
feature_y <- "GRN"          # Replace with desired feature for y-axis
feature_size <- "DC_Mac_CD209_over_all_cell_count_FOV_immune_prop"  # Replace with desired feature for size
#feature_size <- "HLADR"  # Replace with desired feature for size


# # Define the features to use in the plot
feature_x <- "fucosylated"  # Replace with desired feature for x-axis
feature_y <- "NRXN1"          # Replace with desired feature for y-axis
feature_size <- "Neurons_over_all_tumor_count_prop"  # Replace with desired feature for size

# Filter and aggregate duplicate values by patient_id
filtered_data <- master_feature_table_filtered %>%
  filter(
    immunotherapy == "no",
    recurrence == "no",
    site == "Stanford",
    feature_variable %in% c(feature_x, feature_y, feature_size),
    !is.na(who_grade) & who_grade != "unknown"
  )


filtered_data <- master_feature_table_filtered %>%
  filter(site == "Stanford",
         who_grade %in% c("2", "3", "4"),
         tumor_region %in% c("tumor_core", "other"),
         feature_variable %in% c(feature_x, feature_y, feature_size))



filtered_data_wide <- filtered_data %>%
  pivot_wider(
    id_cols = c(patient_id, who_grade),  # Include patient_id and who_grade as identifiers
    names_from = feature_variable,
    values_from = feature_value,
    values_fn = mean  # Aggregate duplicates by mean if needed
  )

# Remove rows where feature_x is 0
filtered_data_wide <- filtered_data_wide %>%
  filter(!is.na(.data[[feature_x]]) & .data[[feature_x]] != 0)
# Remove rows with NA in any column
filtered_data_wide <- filtered_data_wide %>%
  filter(complete.cases(.))

# Calculate unweighted correlation coefficient
correlation_unweighted <- round(cor(
  filtered_data_wide[[feature_x]],
  filtered_data_wide[[feature_y]],
  use = "complete.obs"
), 2)

# Calculate weighted correlation coefficient
# Function to calculate weighted correlation
weighted_correlation <- function(x, y, w) {
  # Center x and y by their weighted means
  mx <- weighted.mean(x, w)
  my <- weighted.mean(y, w)
  x <- x - mx
  y <- y - my
  
  # Compute weighted covariance and variance
  cov_xy <- sum(w * x * y)
  var_x <- sum(w * x^2)
  var_y <- sum(w * y^2)
  
  # Calculate weighted correlation
  cov_xy / sqrt(var_x * var_y)
}

# Apply the function to get the weighted correlation
correlation_weighted <- round(weighted_correlation(
  filtered_data_wide[[feature_x]], 
  filtered_data_wide[[feature_y]], 
  filtered_data_wide[[feature_size]]
), 2)

correlation_unweighted
correlation_weighted

ggplot(filtered_data_wide, aes_string(x = feature_x, y = feature_y, size = feature_size, color = "who_grade")) +
  geom_point(alpha = 1) +
  geom_smooth(
    method = "lm", 
    aes(weight = .data[[feature_size]]),  # Add weights for size
    color = "black", 
    se = FALSE, 
    show.legend = FALSE, 
    size = 0.5
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_color_manual(
    values = c("2" = "#4D9A60", "3" = "#D4854C", "4" = "#9E3850")  # Assign specific colors
  )+
  labs(
    x = NULL,
    y = NULL,
    size = feature_size
  ) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.25, color = "black"),
    axis.ticks = element_line(size = 0.25),
    axis.text = element_text(color = "black"),   # Make axis text black
    axis.title = element_text(color = "black"), # Make axis titles black
    plot.title = element_blank(),
    legend.position = "none"
  )

file_path <- file.path(
  output_dir,
  sprintf("scatter_plot_%s_vs_%s_size_%s.png", feature_x, feature_y, feature_size)
)

# Export the plot as a PNG with the specified path, dimensions, and resolution
ggsave(file_path, width = 45, height = 30, units = "mm", dpi = 600)




###### Figure 7 ################################################################
# Panel b - AUC plot (WHO grade) ----


# Rename column 'Label' to 'Type'
names(auc_df_WHO)[names(auc_df_WHO) == "Label"] <- "Type"

# Rename values in the 'Type' column from 'Randomized' to 'Permuted'
auc_df_WHO$Type[auc_df_WHO$Type == "Randomized"] <- "Permuted"
auc_df_WHO
# Perform a paired t-test using auc_df
t_test_result <- t.test(AUC ~ Type, data = auc_df_WHO, paired = TRUE)

# Extract the p-value
p_value <- t_test_result$p.value

# Determine if the p-value is significant
significance_label <- ifelse(p_value < 0.05, "*", "")
auc_plot <- ggplot(auc_df_WHO, aes(x = Type, y = AUC, fill = Type)) +
  geom_boxplot(outlier.shape = NA, aes(alpha = 1), size = 0.25) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.7, size = 0.2) +
  labs(x = "", y = "AUC") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.25),  # Show x and y axis lines with thickness 0.25
    axis.text.x = element_text(size = 7, family = "Helvetica", color = "black"),  # Set font for x-axis text
    axis.text.y = element_text(size = 7, family = "Helvetica", color = "black"),  # Set font for y-axis text
    axis.title.y = element_text(size = 7, family = "Helvetica"),  # Set font for y-axis label
    axis.title.x = element_blank(),
    legend.position = "none",
    #plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.border = element_blank(),  # Remove outer border
    plot.title = element_text(size = 7, family = "Helvetica"),  # Set font for plot title (if used)
    plot.subtitle = element_text(size = 7, family = "Helvetica"),  # Set font for plot subtitle (if used)
    plot.caption = element_text(size = 7, family = "Helvetica"),  # Set font for plot caption (if used)
    axis.ticks = element_line(size = 0.25),  # Set thickness of tick marks
    axis.ticks.length = unit(0.5, "mm")  # Set length of tick marks
  ) +
  scale_fill_manual(values = c("Original" = "#85B7CE", "Permuted" = "white")) +
  scale_alpha_identity() +
  geom_segment(aes(x = 1, y = 0.99, xend = 2, yend = 0.99), color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 0.99, label = significance_label, size = 3, family = "Helvetica", hjust = 0.5) +
  ylim(NA, 1.0)

# Display the plot
print(auc_plot)

ggsave(
  filename = file.path(output_dir, "AUC_who_grade.png"),
  plot = auc_plot,
  dpi = 1200,
  width = 35,
  height = 35,
  units = "mm"
)

# Panel c - median feature importance plot (WHO grade) ----


# Filter out rows where Feature starts with "Unassigned_"
importance_df_filtered <- importance_df_filtered %>%
  filter(!grepl("^Unassigned_", Feature))

importance_df_filtered <- importance_df_filtered %>%
  mutate(combined_feature_type = Feature_Type,  # Copy everything from Feature_Type
         combined_feature_type = ifelse(Feature_Type == "RNA" & grepl("Low", Feature), 
                                        "CD45- cell RNA", 
                                        ifelse(Feature_Type == "RNA" & grepl("Immune", Feature), 
                                               "CD45+ cell RNA", 
                                               combined_feature_type)))
# First, calculate the median importance for each 'Features_broad' in 'importance_df_filtered'
importance_df_filtered <- importance_df_filtered %>%
  group_by(Features_broad) %>%
  mutate(median_importance = median(Importance, na.rm = TRUE)) %>%
  ungroup()


# Now proceed with 'importance_df_filtered_detailed' creation and remove duplicates
importance_df_filtered_detailed <- importance_df_filtered %>%
  mutate(Features_broad = gsub("_", " ", Features_broad),
         Features_broad = ifelse(Feature_Type == "Glycan" & is.na(Features_broad), 
                                 "Complex", 
                                 ifelse(is.na(Features_broad), Feature_type_two, Features_broad))) %>%
  group_by(Features_broad) %>%
  mutate(median_importance = median(Importance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Features_broad = reorder(Features_broad, median_importance)) %>%
  distinct(Feature, .keep_all = TRUE)  # Remove duplicates based on 'Feature' column

# View the updated dataframe
head(importance_df_filtered_detailed)

importance_df_filtered_detailed <- importance_df_filtered_detailed %>%
  group_by(combined_feature_type) %>%
  mutate(sum_importance = sum(Importance, na.rm = TRUE)) %>%
  ungroup()


top_75 <- importance_df_filtered_detailed %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 75) %>%
  select(Importance, Feature_Type, Feature_type_two, Features_broad, combined_feature_type)

write.csv(
  top_75,
  file = file.path(output_dir, "top_75.csv"),
  row.names = FALSE
)


# Define the custom colors for each Feature_Type
annotation_colors <- c(
  "MIBI" = "#d5405e",
  "Glycan" = "#b1dfdb",
  "CD45- cell RNA" = "#00429d",
  "CD45+ cell RNA" = "#a3c1d9"  # Similar shade of blue to "Selected RNA"
)

# View the updated colors
annotation_colors

# Identify the top 25 Importance scores and create a new column
importance_df_filtered_detailed <- importance_df_filtered_detailed %>%
  arrange(desc(Importance)) %>%
  mutate(is_top_25 = ifelse(row_number() <= 25, TRUE, FALSE))  # Create a flag for the top 25 rows

importance_df_filtered_median <- importance_df_filtered_detailed %>%
  arrange(desc(Importance))%>%
  slice(1:75)


# Check if "Glycan" exists in combined_feature_type and add it if not
if (!"Glycan" %in% importance_df_filtered_median$combined_feature_type) {
  importance_df_filtered_median <- importance_df_filtered_median %>%
    add_row(combined_feature_type = "Glycan", Importance = 0, is_top_25 = FALSE)
}


# Update the plot
plot <- ggplot(importance_df_filtered_median, aes(x = reorder(combined_feature_type, Importance, FUN = median), y = Importance, fill = combined_feature_type)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, size = 0.2, color = "black") +  # Add box plot
  geom_point(aes(color = is_top_25), position = position_jitter(width = 0.2), size = 0.1, show.legend = FALSE) +  # Add jittered points, color top 25 red
  coord_flip() +  # Flip the coordinates
  labs(x = "", y = "Log10 Importance Score", title = "", fill = "Feature Type") +  # Rename legend title
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +  # Custom color palette for points
  scale_fill_manual(values = annotation_colors) +  # Custom color palette for box plots
  scale_y_log10() +  # Apply log scale to the y-axis (Importance)
  theme(
    text = element_text(size = 7),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "none",  # Remove the legend
    panel.grid.major.y = element_line(color = "black", size = 0.05),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

plot


ggsave(
  filename = file.path(output_dir, "detailed_features_who_grade.png"),
  plot = plot,
  width = 80,
  height = 42,
  units = "mm",
  dpi = 1200
)

# Panel d - top 25 feature importance plot (WHO grade) ----



importance_df_filtered_detailed_plot <- importance_df_filtered_detailed %>%
  arrange(desc(Importance))%>%
  slice(1:25)

importance_df_filtered_detailed_plot <- importance_df_filtered_detailed_plot[order(importance_df_filtered_detailed_plot$Importance),]

importance_df_filtered_detailed_plot$Feature_type_two <- factor(importance_df_filtered_detailed_plot$Feature_type_two, levels = unique(importance_df_filtered_detailed_plot$Feature_type_two))

# Update the plot to use a linear scale and sort by Importance
plot <- ggplot(importance_df_filtered_detailed_plot, aes(x = Feature_type_two, y = Importance, fill = combined_feature_type)) +
  geom_bar(stat = "identity", alpha = 1, size = 0.2) +  # Use geom_bar for a bar plot
  coord_flip() +  # Flip the coordinates to make it easier to read the feature names
  labs(x = "", y = "Importance score", title = "", fill = "Feature Type") +  # Rename legend title
  theme_minimal() +
  scale_fill_manual(values = annotation_colors) +  # Use custom color palette for bars
  theme(
    text = element_text(size = 7),  # Adjust font size for all text
    axis.title = element_text(size = 7, color = "black"),  # Adjust font size for axis titles
    axis.text = element_text(size = 7, color = "black"),  # Adjust font size for axis text
    axis.text.y = element_blank(),
    legend.position = "none",  # Remove the legend
    panel.grid.major.y = element_line(color = "black", size = 0.05),  # Add gray y-axis grid lines
    panel.grid.minor.y = element_blank(),  # Remove minor y-axis grid lines
    panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor x-axis grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Add outline around the plot
  )

plot

ggsave(
  filename = file.path(output_dir, "type_features_who_grade.png"),
  plot = plot,
  width = 40,
  height = 90,
  units = "mm",
  dpi = 1200
)

# Order the dataframe by descending Importance
importance_df_filtered_detailed_plot <- importance_df_filtered_detailed_plot %>%
  arrange(desc(Importance))

# Update the levels of Combined_Feature to reflect the sorted order
importance_df_filtered_detailed_plot$Combined_Feature <- factor(
  importance_df_filtered_detailed_plot$Combined_Feature,
  levels = rev(importance_df_filtered_detailed_plot$Combined_Feature)
)

# Define the custom colors for each combined_feature_type
annotation_colors <- c(
  "MIBI" = "#d5405e",
  "Glycan" = "#b1dfdb",
  "CD45- cell RNA" = "#00429d",
  "CD45+ cell RNA" = "#a3c1d9"
)

# Create the annotation plot
annotation_plot <- ggplot(importance_df_filtered_detailed_plot, aes(y = Combined_Feature, x = 1)) +
  geom_tile(aes(fill = combined_feature_type), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = annotation_colors) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

# Display the plot
print(annotation_plot)


ggsave(
  filename = file.path(output_dir, "annotation_plot_who_grade.png"),
  plot = annotation_plot,
  width = 5,       # Width in mm
  height = 76,     # Height in mm
  units = "mm",    # Units of measurement
  dpi = 1200       # Resolution in DPI
)


# Panel e - CD31 and VEGFA boxplot ----



filtered_df <- master_feature_table_filtered %>%
  filter(who_grade %in% c("2","3","4"),
         site == "Stanford",
         tumor_region %in% c("tumor_core", "other"),
         bio_feature_type %in% c("protein_intensity"),
         cell_meta_cluster_final == "Endothelial_cells",
         feature_variable == "CD31")


# Plot the Shannon diversity index by final_diagnosis_simple
p <- ggplot(filtered_df, aes(x = who_grade, y = feature_value)) +
  geom_boxplot(fill = "#d5405e", color = "black", outlier.shape = NA, size = 0.25) +
  geom_jitter(width = 0.2, size = 0.05, alpha = 1, color = "black") +
  labs(title = "",
       x = "WHO Grade",
       y = "Intensity") +
  theme_minimal() +
  theme(text = element_text(family = "Helvetica", size = 7, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5,family = "Helvetica",  color = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 7,family = "Helvetica",  color = "black"),
        panel.grid = element_blank(),     # Remove grid lines
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(size = 0.25),  # Set thickness of tick marks
        axis.ticks.length = unit(0.5, "mm"))  # Set length of tick marks  # Add axis lines

p

ggsave(
  filename = file.path(output_dir, "CD31_intensity_who_grade_plot.png"),
  plot = p,
  dpi = 600,
  height = 43,
  width = 40,
  units = "mm"
)





filtered_df <- master_feature_table_filtered %>%
  filter(who_grade %in% c("2","3","4"),
         site == "Stanford",
         tumor_region %in% c("tumor_core", "other"),
         bio_feature_type %in% c("spatial_RNA"),
         feature_type == "Immune_Low",
         feature_variable == "VEGFA")




# Box plot with specific fill colors for Non and Immune feature_source, jitter points, and no legend
p <- ggplot(filtered_df, aes(x = who_grade, y = feature_value, fill = feature_type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), size = 0.25) +
  geom_jitter(aes(color = feature_type), size = 0.1, alpha = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)) +  # Jitter with dodge
  scale_fill_manual(values = c("Immune_Low" = "#00429d", "Immune_High" = "#a3c1d9")) +  # Custom fill colors for Non and Immune
  scale_color_manual(values = c("Immune_Low" = "black", "Immune_High" = "black")) +  # Match jitter point colors with fill
  labs(title = "", x = "WHO Grade", y = "Counts") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 7, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, family = "Helvetica", color = "black"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 6, family = "Helvetica", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", linewidth = 0.25),  # Add axis lines
    axis.ticks = element_line(size = 0.25),  # Set thickness of tick marks
    axis.ticks.length = unit(0.5, "mm")  # Set length of tick marks  # Add axis lines
  ) +
  guides(fill = "none", color = "none")  # Remove legends for both fill and color

p

ggsave(
  filename = file.path(output_dir, "VEGFA_who_grade_plot.png"),
  plot = p,
  dpi = 600,
  height = 43,
  width = 40,
  units = "mm"
)





# Panel h - AUC plot (survival) ----


# Rename column 'Label' to 'Type'
names(auc_df_GBM)[names(auc_df_GBM) == "Label"] <- "Type"

# Rename values in the 'Type' column from 'Randomized' to 'Permuted'
auc_df_GBM$Type[auc_df_GBM$Type == "Randomized"] <- "Permuted"
auc_df_GBM
# Perform a paired t-test using auc_df
t_test_result <- t.test(AUC ~ Type, data = auc_df_GBM, paired = TRUE)

# Extract the p-value
p_value <- t_test_result$p.value

# Determine if the p-value is significant
significance_label <- ifelse(p_value < 0.05, "*", "")
auc_plot <- ggplot(auc_df_GBM, aes(x = Type, y = AUC, fill = Type)) +
  geom_boxplot(outlier.shape = NA, aes(alpha = 1), size = 0.25) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.7, size = 0.2) +
  labs(x = "", y = "AUC") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.25),  # Show x and y axis lines with thickness 0.25
    axis.text.x = element_text(size = 7, family = "Helvetica", color = "black"),  # Set font for x-axis text
    axis.text.y = element_text(size = 7, family = "Helvetica", color = "black"),  # Set font for y-axis text
    axis.title.y = element_text(size = 7, family = "Helvetica"),  # Set font for y-axis label
    axis.title.x = element_blank(),
    legend.position = "none",
    #plot.margin = unit(c(1, 1, 2, 1), "lines"),
    panel.border = element_blank(),  # Remove outer border
    plot.title = element_text(size = 7, family = "Helvetica"),  # Set font for plot title (if used)
    plot.subtitle = element_text(size = 7, family = "Helvetica"),  # Set font for plot subtitle (if used)
    plot.caption = element_text(size = 7, family = "Helvetica"),  # Set font for plot caption (if used)
    axis.ticks = element_line(size = 0.25),  # Set thickness of tick marks
    axis.ticks.length = unit(0.5, "mm")  # Set length of tick marks
  ) +
  scale_fill_manual(values = c("Original" = "#85B7CE", "Permuted" = "white")) +
  scale_alpha_identity() +
  geom_segment(aes(x = 1, y = 0.99, xend = 2, yend = 0.99), color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = 0.99, label = significance_label, size = 3, family = "Helvetica", hjust = 0.5) +
  ylim(NA, 1.0)

# Display the plot
print(auc_plot)

ggsave(
  filename = file.path(output_dir, "AUC_survival.png"),
  plot = auc_plot,
  dpi = 1200,
  width = 35,
  height = 35,
  units = "mm"
)

# Panel i - median feature importance plot (survival) ----


importance_df_filtered_survival <- importance_df_filtered_survival %>%
  mutate(combined_feature_type = Feature_Type,  # Copy everything from Feature_Type
         combined_feature_type = ifelse(Feature_Type == "RNA" & grepl("Low", Feature), 
                                        "CD45- cell RNA", 
                                        ifelse(Feature_Type == "RNA" & grepl("Immune", Feature), 
                                               "CD45+ cell RNA", 
                                               combined_feature_type)))
# First, calculate the median importance for each 'Features_broad' in 'importance_df_filtered_survival'
importance_df_filtered_survival <- importance_df_filtered_survival %>%
  group_by(Features_broad) %>%
  mutate(median_importance = median(Importance, na.rm = TRUE)) %>%
  ungroup()


# Now proceed with 'importance_df_filtered_detailed' creation and remove duplicates
importance_df_filtered_detailed <- importance_df_filtered_survival %>%
  mutate(Features_broad = gsub("_", " ", Features_broad),
         Features_broad = ifelse(Feature_Type == "Glycan" & is.na(Features_broad), 
                                 "Complex", 
                                 ifelse(is.na(Features_broad), Feature_type_two, Features_broad))) %>%
  group_by(Features_broad) %>%
  mutate(median_importance = median(Importance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Features_broad = reorder(Features_broad, median_importance)) %>%
  distinct(Feature, .keep_all = TRUE)  # Remove duplicates based on 'Feature' column

# View the updated dataframe
head(importance_df_filtered_detailed)

importance_df_filtered_detailed <- importance_df_filtered_detailed %>%
  group_by(combined_feature_type) %>%
  mutate(sum_importance = sum(Importance, na.rm = TRUE)) %>%
  ungroup()


importance_df_filtered_detailed <- importance_df_filtered_detailed %>%
  group_by(combined_feature_type) %>%
  mutate(sum_importance = sum(Importance, na.rm = TRUE)) %>%
  ungroup()


top_75 <- importance_df_filtered_detailed %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 75) %>%
  select(Importance, Feature_Type, Feature_type_two, Features_broad, combined_feature_type)


write.csv(
  top_75,
  file = file.path(output_dir, "top_75_survival.csv"),
  row.names = FALSE
)



# Define the custom colors for each Feature_Type
annotation_colors <- c(
  "MIBI" = "#d5405e",
  "Glycan" = "#b1dfdb",
  "CD45- cell RNA" = "#00429d",
  "CD45+ cell RNA" = "#a3c1d9"  # Similar shade of blue to "Selected RNA"
)

# View the updated colors
annotation_colors

# Identify the top 25 Importance scores and create a new column
importance_df_filtered_detailed <- importance_df_filtered_detailed %>%
  arrange(desc(Importance)) %>%
  mutate(is_top_25 = ifelse(row_number() <= 25, TRUE, FALSE))  # Create a flag for the top 25 rows


importance_df_filtered_median <- importance_df_filtered_detailed %>%
  arrange(desc(Importance))%>%
  slice(1:75)

# Check if "Glycan" exists in combined_feature_type and add it if not
if (!"Glycan" %in% importance_df_filtered_median$combined_feature_type) {
  importance_df_filtered_median <- importance_df_filtered_median %>%
    add_row(combined_feature_type = "Glycan", Importance = 0, is_top_25 = FALSE)
}





# Update the plot
plot <- ggplot(importance_df_filtered_median, aes(x = reorder(combined_feature_type, Importance, FUN = median), y = Importance, fill = combined_feature_type)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, size = 0.2, color = "black") +  # Add box plot
  geom_point(aes(color = is_top_25), position = position_jitter(width = 0.2), size = 0.1, show.legend = FALSE) +  # Add jittered points, color top 25 red
  coord_flip() +  # Flip the coordinates
  labs(x = "", y = "Log10 Importance Score", title = "", fill = "Feature Type") +  # Rename legend title
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +  # Custom color palette for points
  scale_fill_manual(values = annotation_colors) +  # Custom color palette for box plots
  scale_y_log10() +  # Apply log scale to the y-axis (Importance)
  theme(
    text = element_text(size = 7),
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "none",  # Remove the legend
    panel.grid.major.y = element_line(color = "black", size = 0.05),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

plot


ggsave(
  filename = file.path(output_dir, "detailed_features_survival_status.png"),
  plot = plot,
  width = 80,
  height = 42,
  units = "mm",
  dpi = 1200
)

# Panel j - PCA plot of survival ----


# Plotting the scatter plot with formatting
p <- ggplot(PCA, aes(x = PC1, y = PC2, color = survival_status)) +
  geom_point(size = 1) +
  theme_minimal() +  # Set font to Helvetica and size to 6
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line.x = element_line(size = 0.25),  # Add X axis line
    axis.line.y = element_line(size = 0.25),  # Add Y axis line
    axis.title.x = element_text(family = "Helvetica", size = 7),  # X axis title
    axis.title.y = element_text(family = "Helvetica", size = 7),  # Y axis title
    axis.text = element_text(family = "Helvetica", size = 7),  # Axis text
    plot.title = element_blank(),  # Remove title
    axis.ticks = element_line(size = 0.25),  # Set thickness of tick marks
    axis.ticks.length = unit(0.5, "mm"),  # Set length of tick marks
    legend.position = "none"  # Remove legend
  ) +
  scale_color_manual(values = c("long" = "#66c2a5", "short" = "#8e3b88")) +
  labs(title = "",
       x = "PC1 (27%)",
       y = "PC2 (9%)")

p

ggsave(
  filename = file.path(output_dir, "PCA_scatter_plot.png"),
  plot = p,
  width = 40,
  height = 38,
  units = "mm",
  dpi = 1200
)

