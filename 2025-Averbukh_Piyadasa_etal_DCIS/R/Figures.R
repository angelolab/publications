#### Misc. ----
# Figures scripts by Hadeesha Piyadasa & Inna Averbukh 
# December 9, 2024
rm(list = ls(all.names = TRUE))
gc()

output_dir <- "your_path/Figures"
#### Load required packages ----
library(data.table)
library(dplyr)
library(scales)
library(tibble)
library(circlize)
library(ggbreak)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(readr)
library(plotly)
library(htmlwidgets)
library(webshot)
library(writexl)
library(ggridges)
library(purrr)

#### Load relavent files ----

cell_table <- fread("your_path/cell_table_cleaned.csv")
metadata <- fread("your_path/DCIS_Metadata_clean.csv")
test_tab <- read.csv("your_path/test_tab_w_all_metadata.csv")

#### Preprocessing ----

# Perform a left join to merge metadata into cell_table
cell_table <- merge(cell_table, metadata, by.x = "fov", by.y = "MIBI_ID", all.x = TRUE)

cell_table_filtered <- copy(cell_table)

cell_table_filtered <- cell_table_filtered %>%
  filter(sample_type == "primary_event")

cell_table_filtered <- cell_table_filtered %>%
  rename(cell_meta_cluster = concise_meta_cluster)


#### Figure 1: Panel C ----

cell_table_filtered <- copy(cell_table)

cell_table_filtered <- cell_table_filtered %>%
  filter(sample_type == "primary_event")

cell_table_filtered <- cell_table_filtered %>%
  rename(cell_meta_cluster = concise_meta_cluster)

filtered_data <- cell_table_filtered %>%
  mutate(cell_meta_cluster_final = case_when(
    grepl("KRT", cell_meta_cluster, ignore.case = TRUE) ~ "Epithelial",
    cell_meta_cluster == "Epithelial" ~ "Epithelial",
    cell_meta_cluster == "Epithelial_Epcam" ~ "Epithelial",
    cell_meta_cluster == "Stromal" ~ "Stromal",
    cell_meta_cluster == "Myoep" ~ "Myoep",
    cell_meta_cluster == "Endothelial" ~ "Endothelial",
    TRUE ~ cell_meta_cluster
  )) %>%
  filter(!(cell_meta_cluster_final %in% c("Undetermined", "Immune_other"))) %>%
  filter(sample_type == "primary_event") %>%
  filter(event_recur_type %in% c("Non-progressor", "Invasive_Ipsilateral"))

filtered_data <- filtered_data %>%
  mutate(cell_meta_cluster_final_broad = case_when(
    cell_meta_cluster_final %in% c(
      "Mono_Mac_CD68", "Mono_Mac_CD163_206", "Mono_Mac_CD14",
      "Mono_Mac_DC_CD209", "Mono_Mac_CD11c", "Neutrophils", "APC", "Mast_cells"
    ) ~ "myeloid",
    
    cell_meta_cluster_final %in% c(
      "Tcell_CD4", "Tcell_CD8", "Tcell_other", "Tregs"
    ) ~ "lymphoid",
    
    cell_meta_cluster_final == "Epithelial" ~ "epithelial",
    cell_meta_cluster_final == "Endothelial" ~ "endothelial",
    cell_meta_cluster_final == "Myoep" ~ "myoepithelial",
    cell_meta_cluster_final == "Stromal" ~ "stromal",
    
    TRUE ~ "other"
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
  "endothelial",
  "epithelial", "stromal","myoepithelial","myeloid","lymphoid"
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


color_broad <- c("#66CDAA", "#C83296", "#0FCB75",  "#FF9600", "#4270E4","#CA0000")


pull_values_broad <- ifelse(summary_data_broad$cell_meta_cluster_final_broad %in% c("lymphoid"), 
                            0.0, 
                            ifelse(summary_data_broad$cell_meta_cluster_final_broad %in% c("myeloid"), 
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

# Save with transparent background
htmlwidgets::saveWidget(as_widget(pie), "fig2.html", selfcontained = FALSE)
webshot("fig2.html", file = file.path(output_dir, "broad_pie.png"), zoom = 10, vwidth = 600, vheight = 600, selector = "div.plot-container")

# Filter data for Immune broad category
immune_data <- filtered_data %>%
  filter(cell_meta_cluster_final_broad == "myeloid") %>%
  count(cell_meta_cluster_final) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(n)) %>%
  mutate(cell_meta_cluster_final = factor(cell_meta_cluster_final, levels = cell_meta_cluster_final))


# Define the color palette
color_final <- c(
  "#b1dfdb",  # Mono_Mac_CD163_206
  "#a3c1d9",  # APC
  "#1f4e79",  # Mono_Mac_CD11c
  "#4a90d9",  # Neutrophils
  "#3cb371",  # Mono_Mac_CD68
  "#6694c1",  # Mast_cells
  "#00429d",  # Mono_Mac_CD14
  "#228b22"   # Mono_Mac_DC_CD209
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
  filename = file.path(output_dir, "mymphoid_stacked.png"),
  plot = plot,
  width = 13,        # Reduce width
  height = 53,       # Reduce height
  units = "mm",
  dpi = 1200         # Keep high resolution
)


# Filter data for Immune broad category
immune_data <- filtered_data %>%
  filter(cell_meta_cluster_final_broad == "lymphoid") %>%
  count(cell_meta_cluster_final) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(n)) %>%
  mutate(cell_meta_cluster_final = factor(cell_meta_cluster_final, levels = cell_meta_cluster_final))


color_final <- c(
  "#ffb3a7",  # T CD4+
  "#b81b4a",  # T CD8+
  "#93003a",  # DN T cells
  "#d5405e"   # Tregs
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
  filename = file.path(output_dir, "lymphoid_stacked.png"),
  plot = plot,
  width = 13,        # Reduce width
  height = 53,       # Reduce height
  units = "mm",
  dpi = 1200         # Keep high resolution
)

#### Figure 1: Panel D ----

columns_of_interest_phenotype <- c(
  "ANXA1", 
  "CD10", "CD11c", "CD14", "CD163", "CD206", "CD209", "CD31", "CD3e", "CD4", "CD45", 
  "CD68", "CD8", "COL1A1", "Calponin1", "Calprotectin", "Chym_Tryp", "Ecadherin", "EpCAM",
   "FAP", "FOXP3",  "Fibronectin", "HLADR", 
  "SMA", "Vimentin")


columns_of_interest_phenotype <- paste0(columns_of_interest_phenotype, "_pred")

# Remove rows where cell_meta_cluster is "Undetermined", "Immune_other", or "Tcell_other"
cell_table_filtered <- cell_table_filtered[
  !cell_meta_cluster %in% c("Undetermined", "Immune_other", "Tcell_other")
]

# Rename specified clusters to "Epithelial"
clusters_to_rename <- c(
  "Epcam_Ecad_KRT18", "Epithelial", "Epcam_Ecad_KRT7_18", 
  "Epithelial_Epcam", "Epcam_Ecad_KRT7_15_18_81", 
  "Epcam_Ecad_KRT18_81", "Epcam_Ecad_KRT7_18_81"
)

cell_table_filtered[cell_meta_cluster %in% clusters_to_rename, cell_meta_cluster := "Epithelial"]


# Then proceed with subsetting and summarizing
subset_df <- cell_table_filtered %>%
  select(cell_meta_cluster, fov, all_of(columns_of_interest_phenotype))


# Aggregate the primary filtered data by `cell_meta_cluster_final` for columns of interest
heatmap_data <- subset_df %>%
  group_by(cell_meta_cluster) %>%
  select(all_of(columns_of_interest_phenotype)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))


# Convert heatmap_data_primary to matrix
heatmap_matrix <- heatmap_data %>%
  ungroup() %>%               # Ensure the data is ungrouped
  select(-cell_meta_cluster) %>%  # Remove the grouping column
  as.matrix()


rownames(heatmap_matrix) <- heatmap_data$cell_meta_cluster

# Rescale each column of heatmap_matrix
heatmap_matrix_scaled <- t(apply(heatmap_matrix, 1, function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}))

# Count the number of cells per cell_meta_cluster
counts_df <- cell_table_filtered %>%
  count(cell_meta_cluster)


png("your_path/Figure_1_panel_A_bar_plot_pred.png",
    width = 25, height = 60, units = "mm", res = 600)

ggplot(counts_df, aes(y = reorder(cell_meta_cluster, n), x = n)) +
  geom_bar(stat = "identity", fill = "#00429D", width = 0.5) +
  scale_y_discrete(position = "right") +
  scale_x_reverse() +
  scale_x_break(c(75000, 75001), scales = c(0.2, 0.8)) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.3, fill = NA),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.bottom = element_line(color = "black", size = 0.1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(0.05, "cm")
  )

dev.off()


# Make sure cell_meta_cluster is a character vector in counts_df
counts_df$cell_meta_cluster <- as.character(counts_df$cell_meta_cluster)

# Arrange counts_df in descending order of counts
counts_df_ordered <- counts_df %>% arrange(desc(n))

# Extract the new order of clusters based on counts and ensure it's a character vector
new_order <- as.character(counts_df_ordered$cell_meta_cluster)

# Check if all elements in new_order are present as rownames in heatmap_matrix_scaled
# If not, ensure they match or adjust accordingly
all(new_order %in% rownames(heatmap_matrix_scaled)) 

# Reorder the rows of the heatmap_matrix_scaled according to the new order
heatmap_matrix_scaled <- heatmap_matrix_scaled[new_order, , drop = FALSE]

ht <- Heatmap(
  heatmap_matrix_scaled, 
  name = " ",
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  col = colorRamp2(c(0, 0.7, 1), c("white", "#DAE9F5", "#00429D")),
  row_dend_width = unit(2, "mm"),
  column_dend_height = unit(2, "mm"),
  gap = unit(5, "mm"),
  rect_gp = gpar(col = "black", lwd = 0.25),
  row_dend_gp = gpar(lwd = 0.25),
  column_dend_gp = gpar(lwd = 0.25),
  show_heatmap_legend = FALSE
)

draw(ht)

pdf("your_path/Figure_1_panel_A_pred.pdf", 
    width = 75/25.4, height = 50/25.4)
draw(ht)

dev.off()


#### Figure 2: Panel C ----


# Define markers, column names, colors, and labels
markers <- c("SMA", "KRT18", "ANXA1")
col_names <- paste0("bin_", markers, "_pred")
colors <- c("#FCB421", "#C24E9D", "#55C9E9")
labels <- paste0(markers, "+")

# Initialize an empty data frame to store combined data
combined_data <- data.frame()

# Loop to prepare data for each marker
for (i in seq_along(markers)) {
  col <- col_names[i]
  marker_label <- labels[i]
  
  # Filter data for the current marker
  filtered_data <- cell_table_filtered %>%
    filter(fov == "TA536_R10C7",
           cell_meta_cluster %in% c("Epithelial", "Myoep"),
           DuctNumber == 1,
           !!sym(col) == 1) %>%
    mutate(Marker = factor(marker_label, levels = labels))  # Ensure consistent factor levels
  
  # Combine with previous data
  combined_data <- bind_rows(combined_data, filtered_data)
}

xlim <- c(0.2, 1.1)

# keep only values inside plotting range to avoid the warning
df <- combined_data %>% filter(r >= xlim[1], r <= xlim[2])

# per-marker densities, scaled to max = 1 within each marker
dens_df <- df %>%
  group_by(Marker) %>%
  reframe({
    d <- density(r, from = xlim[1], to = xlim[2], n = 512)
    tibble(x = d$x, y = d$y / max(d$y))
  })

p <- ggplot(dens_df, aes(x = x, y = Marker, height = y, fill = Marker)) +
  geom_ridgeline(scale = 0.9, alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = setNames(colors, labels)) +
  scale_x_continuous(
    limits = c(0.2, 1.1),
    breaks = seq(0.2, 1.1, by = 0.3)
  ) +
  coord_cartesian(ylim = c(1.5, 3.5)) +  # control vertical spacing only
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = NA),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 6),
    plot.title = element_blank(),
    legend.position = "none"
  )

print(p)
# Optional: Save the plot to file
ggsave(
  filename = "your_path/TA536_R10C7_density_plot.png", 
  plot = p, width = 42, height = 30, units = "mm",dpi = 600)


# Define markers, column names, colors, and labels
markers <- c("SMA", "Ecadherin", "KRT5")
col_names <- paste0("bin_", markers, "_pred")
colors <- c("#FCB421", "#C24E9D", "#55C9E9")
labels <- paste0(markers, "+")


# Initialize an empty data frame to store combined data
combined_data <- data.frame()

# Loop to prepare data for each marker
for (i in seq_along(markers)) {
  col <- col_names[i]
  marker_label <- labels[i]
  
  # Filter data for the current marker
  filtered_data <- cell_table_filtered %>%
    filter(fov == "TA536_R6C3",
           cell_meta_cluster %in% c("Epithelial", "Myoep"),
           DuctNumber == 6,
           !!sym(col) == 1) %>%
    mutate(Marker = factor(marker_label, levels = labels))  # Ensure consistent factor levels
  
  # Combine with previous data
  combined_data <- bind_rows(combined_data, filtered_data)
}



xlim <- c(0.2, 1.1)

# keep only values inside plotting range to avoid the warning
df <- combined_data %>% filter(r >= xlim[1], r <= xlim[2])

# per-marker densities, scaled to max = 1 within each marker
dens_df <- df %>%
  group_by(Marker) %>%
  reframe({
    d <- density(r, from = xlim[1], to = xlim[2], n = 512)
    tibble(x = d$x, y = d$y / max(d$y))
  })

p <- ggplot(dens_df, aes(x = x, y = Marker, height = y, fill = Marker)) +
  geom_ridgeline(scale = 0.9, alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = setNames(colors, labels)) +
  scale_x_continuous(
    limits = c(0.2, 1.1),
    breaks = seq(0.2, 1.1, by = 0.3)
  ) +
  coord_cartesian(ylim = c(1.5, 3.5)) +  # control vertical spacing only
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = NA),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 6),
    plot.title = element_blank(),
    legend.position = "none"
  )

print(p)

# Optional: Save the plot to file
ggsave(
  filename = "your_path/TA536_R6C3_density_plot.png", 
  plot = p, width = 42, height = 30, units = "mm",dpi = 600)


# Define markers, column names, colors, and labels
markers <- c("SMA", "KRT7", "KRT17")
col_names <- paste0("bin_", markers, "_pred")
colors <- c("#FCB421", "#C24E9D", "#55C9E9")
labels <- paste0(markers, "+")


# Initialize an empty data frame to store combined data
combined_data <- data.frame()

# Loop to prepare data for each marker
for (i in seq_along(markers)) {
  col <- col_names[i]
  marker_label <- labels[i]
  
  # Filter data for the current marker
  filtered_data <- cell_table_filtered %>%
    filter(fov == "TA583_R9C4",
           cell_meta_cluster %in% c("Epithelial", "Myoep"),
           DuctNumber == 3,
           !!sym(col) == 1) %>%
    mutate(Marker = factor(marker_label, levels = labels))  # Ensure consistent factor levels
  
  # Combine with previous data
  combined_data <- bind_rows(combined_data, filtered_data)
}



xlim <- c(0.2, 1.1)

# keep only values inside plotting range to avoid the warning
df <- combined_data %>% filter(r >= xlim[1], r <= xlim[2])

# per-marker densities, scaled to max = 1 within each marker
dens_df <- df %>%
  group_by(Marker) %>%
  reframe({
    d <- density(r, from = xlim[1], to = xlim[2], n = 512)
    tibble(x = d$x, y = d$y / max(d$y))
  })

p <- ggplot(dens_df, aes(x = x, y = Marker, height = y, fill = Marker)) +
  geom_ridgeline(scale = 0.9, alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = setNames(colors, labels)) +
  scale_x_continuous(
    limits = c(0.2, 1.1),
    breaks = seq(0.2, 1.1, by = 0.3)
  ) +
  coord_cartesian(ylim = c(1.5, 3.5)) +  # control vertical spacing only
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = NA),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 6),
    plot.title = element_blank(),
    legend.position = "none"
  )

print(p)

# Optional: Save the plot to file
ggsave(
  filename = "your_path/TA583_R9C4_density_plot.png", 
  plot = p, width = 42, height = 30, units = "mm",dpi = 600)



# Define markers, column names, colors, and labels
markers <- c("SMA", "KRT14", "KRT81")
col_names <- paste0("bin_", markers, "_pred")
colors <- c("#FCB421", "#C24E9D", "#55C9E9")
labels <- paste0(markers, "+")


# Initialize an empty data frame to store combined data
combined_data <- data.frame()

# Loop to prepare data for each marker
for (i in seq_along(markers)) {
  col <- col_names[i]
  marker_label <- labels[i]
  
  # Filter data for the current marker
  filtered_data <- cell_table_filtered %>%
    filter(fov == "TA536_R12C4",
           cell_meta_cluster %in% c("Epithelial", "Myoep"),
           DuctNumber == 14,
           !!sym(col) == 1) %>%
    mutate(Marker = factor(marker_label, levels = labels))  # Ensure consistent factor levels
  
  # Combine with previous data
  combined_data <- bind_rows(combined_data, filtered_data)
}



xlim <- c(0.2, 1.1)

# keep only values inside plotting range to avoid the warning
df <- combined_data %>% filter(r >= xlim[1], r <= xlim[2])

# per-marker densities, scaled to max = 1 within each marker
dens_df <- df %>%
  group_by(Marker) %>%
  reframe({
    d <- density(r, from = xlim[1], to = xlim[2], n = 512)
    tibble(x = d$x, y = d$y / max(d$y))
  })

p <- ggplot(dens_df, aes(x = x, y = Marker, height = y, fill = Marker)) +
  geom_ridgeline(scale = 0.9, alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = setNames(colors, labels)) +
  scale_x_continuous(
    limits = c(0.2, 1.1),
    breaks = seq(0.2, 1.1, by = 0.3)
  ) +
  coord_cartesian(ylim = c(1.5, 3.5)) +  # control vertical spacing only
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = NA),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 6),
    plot.title = element_blank(),
    legend.position = "none"
  )

print(p)

# Optional: Save the plot to file
ggsave(
  filename = "your_path/TA536_R12C4_density_plot.png", 
  plot = p, width = 42, height = 30, units = "mm",dpi = 600)



#### Figure 2: Panel D ----


# Read in the CSV file
test_tab <- read.csv("your_path/test_tab.csv")

# Print column names
colnames(test_tab)



# Vector of features (column names) you expect
features <- c("compartment_score_st")


# Check which columns are missing
missing_cols <- setdiff(features, colnames(test_tab))

if (length(missing_cols) > 0) {
  cat("The following columns do not exist in 'test_tab':\n")
  print(missing_cols)
} else {
  cat("All features exist in 'test_tab'.\n")
}


features <- c("compartment_score_st")

# Loop over features
for (feature in features) {
  if (feature %in% colnames(test_tab)) {
    
    # Create the plot
    p <- ggplot(test_tab, aes_string(x = feature)) +
      geom_histogram(bins = 30, fill = "black", color = "white", linewidth = 0.1) +
      xlim(-0.8, 0.8) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 22)) +   # y-axis starts at 0
      theme_minimal() +
      theme(
        text = element_text(size = 7),
        axis.title = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.25),
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length = unit(1, "mm")  # short ticks inside
      )
    
    # Define PNG export path
    out_file <- paste0(
      "your_path/",
      "Compartment_", feature, ".png"
    )
    
    # Save as high-resolution PNG (1200 dpi)
    ggsave(
      filename = out_file,
      plot = p,
      width = 45,
      height = 40,
      units = "mm",
      dpi = 1200
    )
  }
}


#### Figure 3: Panel A ----


# Define the path to the folder
folder_path <- "your_path/Lasso classifier tables/"

# List all relevant CSV files
file_names <- sprintf("Lasso_roc_table_run_%dwith_compartment_score_st.csv", 1:10)
file_paths <- file.path(folder_path, file_names)

# Read and combine all files
df_list <- lapply(1:10, function(i) {
  df <- read_csv(file_paths[i]) %>%
    mutate(run = paste0("Run ", i),  # Add run identifier
           minus_spec = 1 - specificity) %>%
    arrange(sensitivity)
  return(df)
})

df_combined <- bind_rows(df_list)

# Separate Run 7 from the rest
df_run7 <- df_combined %>% filter(run == "Run 7")
df_rest <- df_combined %>% filter(run != "Run 7")

# Define output PDF path
pdf("your_path/Figure_4_panel_A_ROC.pdf", 
    width = 35/25.4, height = 60/25.4)

# Plot all runs in light grey, then overlay Run 7 in black
ggplot() +
  # Individual runs in light grey
  geom_step(data = df_rest, aes(x = minus_spec, y = sensitivity, group = run), 
            direction = "vh", color = "lightgrey", size = 0.25) +
  # Run 7 ROC curve in black
  geom_step(data = df_run7, aes(x = minus_spec, y = sensitivity), 
            direction = "vh", color = "black", size = 0.4) +
  # Diagonal reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.25) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_line(size = 0.25),
    axis.text = element_text(size = 7),  # Set axis text size to 8
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    
    
  )

dev.off()

#### Figure 3: Panel B ----

lasso_features <- read.csv("your_path/classifier_coeffs_with_compartment_score.csv")


lasso_features <- lasso_features %>%
  arrange(desc(mean_coeff)) %>%
  #slice(1:90) %>% if you want to generate the supp table
  slice(1:25) %>% # for the figure
  mutate(Predictor = factor(Predictor, levels = (Predictor)))  # Ensures proper order in plot

# Data frame is 'lasso_features'
coeff_columns <- c("Coefficient_x", "Coefficient_y", "Coefficient_x_1", "Coefficient_y_1", 
                   "Coefficient_x_2", "Coefficient_y_2", "Coefficient_x_3", "Coefficient_y_3", 
                   "Coefficient_x_4", "Coefficient_y_4")

# Calculate the mean of non-zero coefficients for each row
lasso_features$mean <- apply(lasso_features[, coeff_columns], 1, function(row) {
  non_zero <- row[row != 0]
  if (length(non_zero) > 0) mean(non_zero) else 0
})

# Display the updated data frame with the new 'mean' column
print(lasso_features[, c("Predictor", "mean", coeff_columns)])


lasso_features <- lasso_features[order(lasso_features$mean),]


lasso_features <- lasso_features %>%
  mutate(type = case_when(
    grepl("_r_", Predictor) ~ "ductal_topology",
    grepl("_freq", Predictor) ~ "compartment_freq",
    grepl("^\\w+\\.\\w+$", Predictor) ~ "distance",  # e.g., APC.Stromal
    grepl("mean_r", Predictor) ~ "ductal_topology",
    TRUE ~ "Other"
  ))

lasso_features$type[lasso_features$Predictor == "median_myoep_EpCAM"] <- "myoep_features"
lasso_features$type[lasso_features$Predictor == "mean_major_axis_length_fap"] <- "morphology_features"
lasso_features$type[lasso_features$Predictor == "stroma_R_Mono_Mac_CD68_immune"] <- "compartment_freq"
lasso_features$type[lasso_features$Predictor == "mean_duct_alignment_fap"] <- "morphology_features"
lasso_features$type[lasso_features$Predictor == "R_Tcell_CD4_immune"] <- "compartment_freq"
lasso_features$type[lasso_features$Predictor == "R_Neutrophils_immune"] <- "compartment_freq"
lasso_features$type[lasso_features$Predictor == "compartment_score_st"] <- "ductal_topology"
lasso_features$type[lasso_features$Predictor == "mean_duct_alignment_col"] <- "morphology_features"


type_colors <- c(
  "ductal_topology" = "#00429d",     # Blue
  "compartment_freq" = "#a3c1d9",        # Green
  "distance" = "#d5405e",             # Red
  "myoep_features" = "#b1dfdb",          # Orange
  "morphology_features" = "#6a3d9a"     # Purple
)

# Sort by mean in descending order (positive first, then negative)
lasso_features <- lasso_features %>%
  arrange((mean))

# Create the bar plot
plot <- ggplot(lasso_features, aes(x = reorder(Predictor, -mean), y = mean, fill = type)) +
  geom_bar(stat = "identity", alpha = 1, size = 0.2, color = "black") +
  labs(x = "", y = "Importance score", title = "", fill = "Feature Type") +
  scale_fill_manual(values = type_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 7),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black", angle = 90),
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_line(color = "black", size = 0.05),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Print the plot
print(plot)

ggsave(
  filename = file.path(output_dir, "type_features_lasso.png"),
  plot = plot,
  width = 120,
  height = 40,
  units = "mm",
  dpi = 1200
)

# Sort by mean in ascending order (negative first, then positive)
lasso_features <- lasso_features %>%
  arrange(mean)

lasso_features <- lasso_features %>%
  mutate(compartment = case_when(
    grepl("duct_alignment", Predictor, ignore.case = TRUE) ~ "stroma",
    grepl("duct", Predictor, ignore.case = TRUE) ~ "duct",
    grepl("_d_", Predictor, ignore.case = TRUE) ~ "duct",
    grepl("compartment", Predictor, ignore.case = TRUE) ~ "duct",
    grepl("mean_r", Predictor, ignore.case = TRUE) ~ "duct",
    grepl("^\\w+\\.\\w+$", Predictor) ~ "whole_fov",
    grepl("stromal", Predictor, ignore.case = TRUE) ~ "stroma",
    grepl("stroma", Predictor, ignore.case = TRUE) ~ "stroma",
    grepl("myoep", Predictor, ignore.case = TRUE) ~ "myoep",
    grepl("mean_major", Predictor, ignore.case = TRUE) ~ "whole_fov",
    Predictor %in% c("R_Tcell_CD4_immune", "R_Neutrophils_immune") ~ "duct",
    TRUE ~ NA_character_
  ))

compartment_colors <- c(
  "duct" = "#f8766d",      # red-orange
  "stroma" = "#7cae00",    # lime green
  "myoep" = "#00bfc4",     # cyan
  "whole_fov" = "white"  # purple-pink
)

# Create the annotation plot
annotation_plot <- ggplot(lasso_features, aes(x = reorder(Predictor, -mean), y = 1)) +
  geom_tile(aes(fill = compartment), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = compartment_colors) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

# Print the plot
print(annotation_plot)


ggsave(
  filename = file.path(output_dir, "annotation_plot.png"),
  plot = annotation_plot,
  width = 110,       # Width in mm
  height = 2,     # Height in mm
  units = "mm",    # Units of measurement
  dpi = 1200       # Resolution in DPI
)



# Select specified columns
lasso_features_subset <- lasso_features %>%
  select(Predictor, mean, type, compartment, std, N_pos)

# Export to CSV
write.csv(lasso_features_subset, 
          file = "your_path/lasso_features_subset.csv", 
          row.names = FALSE)



#### Figure 3: Panel C,D,E ----

# Read in the CSV file
test_tab <- read.csv("your_path/test_tab_w_all_metadata.csv")

# Print column names
colnames(test_tab)



# Vector of features (column names) you expect
features <- c(#"Neutrophils.Mono_Mac_CD11c",
              "APC.Myoep",
              "pFreq_d_KRT14_KRT17_r_0",
              "pFreq_d_KRT18_KRT7_r_0.9")

# Check which columns are missing
missing_cols <- setdiff(features, colnames(test_tab))

if (length(missing_cols) > 0) {
  cat("The following columns do not exist in 'test_tab':\n")
  print(missing_cols)
} else {
  cat("All features exist in 'test_tab'.\n")
}


# Define custom colors
custom_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Pastel purple
                   "Non-progressor" = "#4C9A2A")       # Pastel green

# Define custom fill colors
custom_fill_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Reddish-purplish pastel
                        "Non-progressor" = "#4C9A2A")       # Darker pastel green


# Loop over each feature
for (feat in features) {
  
  # Run a Kruskal-Wallis test on the current feature grouped by 'event_recur_type'
  kw_result <- kruskal.test(
    x = test_tab[[feat]],
    g = test_tab[["event_recur_type"]]
  )
  
  # Print out the p-value to the console
  cat("Kruskal-Wallis p-value for", feat, ":", kw_result$p.value, "\n")
  
  # Construct output PDF file name based on the feature
  out_file <- paste0(
    "your_path/",
    "Figure_4_panel_B_", feat, ".pdf"
  )

  # Open PDF device
  pdf(out_file, width = 40/25.4, height = 40/25.4)
  
  # Build the plot
  p <- ggplot(test_tab, aes(
    x = .data[["event_recur_type"]],
    y = .data[[feat]],
    fill = .data[["event_recur_type"]],
    color = .data[["event_recur_type"]]  # Match outline to fill
  )) +
    geom_violin(trim = FALSE, size = 0.5, alpha = 0.6) +  # Fill & outline match, transparency applied
    geom_point(position = position_jitter(width = 0.2), size = 0.25, color = "black") +  # Keep dots black
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
    scale_color_manual(values = custom_colors) + # Apply same colors to outline
    labs(x = NULL, y = NULL, fill = NULL, color = NULL) +
    theme_minimal() +
    theme(
      text = element_text(family = "Helvetica", size = 8, face = "plain", colour = "black"),
      #axis.title.x = element_text(),
      #axis.title.y = element_text(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.ticks.y = element_line(size = 0.1, color = "black"),
      axis.ticks.length = unit(0.5, "mm"),
      legend.position = "none",
      plot.title = element_blank()
    )
  
  # Print the plot so it actually gets rendered
  print(p)
  
  # Close the PDF device
  dev.off()
}


#### Figure 3: Panel F ----


# Read in the CSV file
test_tab <- read.csv("your_path/test_tab_w_all_metadata.csv")

# Print column names
colnames(test_tab)



# Vector of features (column names) you expect
features <- c("compartment_score_st")


# Check which columns are missing
missing_cols <- setdiff(features, colnames(test_tab))

if (length(missing_cols) > 0) {
  cat("The following columns do not exist in 'test_tab':\n")
  print(missing_cols)
} else {
  cat("All features exist in 'test_tab'.\n")
}


# Define custom colors
custom_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Pastel purple
                   "Non-progressor" = "#4C9A2A")       # Pastel green

# Define custom fill colors
custom_fill_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Reddish-purplish pastel
                        "Non-progressor" = "#4C9A2A")       # Darker pastel green


# Loop over each feature
for (feat in features) {
  
  # Run a Kruskal-Wallis test on the current feature grouped by 'event_recur_type'
  kw_result <- kruskal.test(
    x = test_tab[[feat]],
    g = test_tab[["event_recur_type"]]
  )
  
  # Print out the p-value to the console
  cat("Kruskal-Wallis p-value for", feat, ":", kw_result$p.value, "\n")
  
  # Construct output PDF file name based on the feature
  out_file <- paste0(
    "your_path",
    "Figure_4_panel_B_", feat, ".pdf"
  )
  
  # Open PDF device
  pdf(out_file, width = 60/25.4, height = 50/25.4)
  
  # Build the plot
  p <- ggplot(test_tab, aes(
    x = .data[["event_recur_type"]],
    y = .data[[feat]],
    fill = .data[["event_recur_type"]],
    color = .data[["event_recur_type"]]  # Match outline to fill
  )) +
    geom_violin(trim = FALSE, size = 0.5, alpha = 0.6) +  # Fill & outline match, transparency applied
    geom_point(position = position_jitter(width = 0.2), size = 0.25, color = "black") +  # Keep dots black
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
    scale_color_manual(values = custom_colors) + # Apply same colors to outline
    labs(x = NULL, y = NULL, fill = NULL, color = NULL) +
    theme_minimal() +
    theme(
      text = element_text(family = "Helvetica", size = 10, face = "plain", colour = "black"),
      axis.title.x = element_text(),
      axis.title.y = element_text(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.ticks.y = element_line(size = 0.1, color = "black"),
      axis.ticks.length = unit(0.5, "mm"),
      legend.position = "none",
      plot.title = element_blank()
    )
  
  
  # Print the plot so it actually gets rendered
  print(p)
  
  # Close the PDF device
  dev.off()
}

p
#### Figure 3: Panel H ----


# Read in the CSV file
test_tab <- read.csv("your_path/test_tab.csv")

# Print column names
colnames(test_tab)



# Vector of features (column names) you expect
features <- c("stroma_Endothelial_freq")

# Check which columns are missing
missing_cols <- setdiff(features, colnames(test_tab))

if (length(missing_cols) > 0) {
  cat("The following columns do not exist in 'test_tab':\n")
  print(missing_cols)
} else {
  cat("All features exist in 'test_tab'.\n")
}


# Define custom colors
custom_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Pastel purple
                   "Non-progressor" = "#4C9A2A")       # Pastel green

# Define custom fill colors
custom_fill_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Reddish-purplish pastel
                        "Non-progressor" = "#4C9A2A")       # Darker pastel green


# Loop over each feature
for (feat in features) {
  
  # Run a Kruskal-Wallis test on the current feature grouped by 'event_recur_type'
  kw_result <- kruskal.test(
    x = test_tab[[feat]],
    g = test_tab[["event_recur_type"]]
  )
  
  # Print out the p-value to the console
  cat("Kruskal-Wallis p-value for", feat, ":", kw_result$p.value, "\n")
  
  # Construct output PDF file name based on the feature
  out_file <- paste0(
    "your_path",
    "Figure_4_panel_B_", feat, ".pdf"
  )
  
  # Open PDF device
  pdf(out_file, width = 60/25.4, height = 50/25.4)
  
  # Build the plot
  p <- ggplot(test_tab, aes(
    x = .data[["event_recur_type"]],
    y = .data[[feat]],
    fill = .data[["event_recur_type"]],
    color = .data[["event_recur_type"]]  # Match outline to fill
  )) +
    geom_violin(trim = FALSE, size = 0.5, alpha = 0.6) +  # Fill & outline match, transparency applied
    geom_point(position = position_jitter(width = 0.2), size = 0.25, color = "black") +  # Keep dots black
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
    scale_color_manual(values = custom_colors) + # Apply same colors to outline
    labs(x = NULL, y = NULL, fill = NULL, color = NULL) +
    theme_minimal() +
    theme(
      text = element_text(family = "Helvetica", size = 10, face = "plain", colour = "black"),
      axis.title.x = element_text(),
      axis.title.y = element_text(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.ticks.y = element_line(size = 0.1, color = "black"),
      axis.ticks.length = unit(0.5, "mm"),
      legend.position = "none",
      plot.title = element_blank()
    )
  
  
  # Print the plot so it actually gets rendered
  print(p)
  
  # Close the PDF device
  dev.off()
}


#### Figure 3: Panel J ----

tailed_quant <- read.csv("your_path/high_outliers_in_lasso_fts.csv")


# Assuming tailed_quant is already loaded
p <- ggplot(tailed_quant, aes(x = skew_cat, fill = skew_cat)) +
  geom_bar() +
  scale_x_discrete(limits = c("Progressors only", "Non-progressors only", "both", "none")) +
  scale_fill_manual(values = c("Non-progressors only" = "#4C9A2A", 
                               "both" = "black", 
                               "none" = "black", 
                               "Progressors only" = "#9B51E0")) +
  labs(x = NULL, y = "Count", title = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.25),  # Border around the plot
        axis.line = element_blank(),  # Remove default axis lines
        axis.ticks = element_line(size = 0.25),  # Make tick marks 0.25 thickness
        legend.position = "none")  # Remove legend since colors are specified

p

print(p)
# Save the plot as PDF
ggsave(filename = "your_path/skew_cat_barplot.pdf",
       plot = p,
       device = "pdf",
       width = 50/25.4,
       height = 30/25.4)

#### Figure 4: Panel K ----

# 1. Read in the existing 'test_tab' data (already done above, but shown here for completeness)
test_tab <- read.csv("your_path/test_tab_w_all_metadata.csv")

# 2. Read in the 'DCIS_meta' file
risk_score <- read.csv("your_path/confidence_score_table.csv")

# 3. Filter 'DCIS_meta' to only rows with 'MIBI_ID' matching the 'fov' column in 'test_tab'
risk_score_filt <- risk_score %>%
  filter(fov %in% test_tab$MIBI_ID_fov)

# Define custom colors
custom_colors <- c("Invasive_Ipsilateral" = "#9B51E0",  # Pastel purple
                   "Non-progressor" = "#4C9A2A")       # Pastel green


# combined violin plot 
p <- ggplot(risk_score_filt, aes(x = "all", y = score)) +
  geom_violin(
    fill = NA,
    color = "black",
    size = 0.3,
    trim = FALSE
  ) +
  geom_point(
  data = subset(risk_score_filt, event_recur_type == "Non-progressor"),
  aes(
    color = event_recur_type,
    size  = 1.2
  ),
  position = position_jitter(width = 0.5),
  alpha = 1,
  stroke = 0
) +
geom_point(
  data = subset(risk_score_filt, event_recur_type == "Invasive_Ipsilateral"),
  aes(
    color = event_recur_type,
    size  = 1.5
  ),
  position = position_jitter(width = 0.5),
  alpha = 1,
  stroke = 0
) +
  
  scale_color_manual(values = custom_colors) +
  scale_size_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 10, colour = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(size = 0.1, color = "black"),
    axis.ticks.length = unit(0.5, "mm"),
    legend.position = "none",
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background  = element_rect(fill = NA, colour = NA)
  )

p



ggsave(
  filename = "your_path/DCIS_meta_score_violin_combined.pdf",
  plot     = p,
  device   = "pdf",
  width    = 60/25.4,
  height   = 40/25.4,
  bg       = "transparent"
)


#### Figure 3: Panel L ----
quant_df <- risk_score_filt %>%
  mutate(
    score_q = ntile(score, 4)  # 4 quantiles
  ) %>%
  group_by(score_q) %>%
  summarise(
    n_total      = n(),
    n_invasive   = sum(event_recur_type == "Invasive_Ipsilateral"),
    pct_invasive = 100 * n_invasive / n_total,
    .groups = "drop"
  )

quant_ranges <- risk_score_filt %>%
  mutate(score_q = ntile(score, 4)) %>%
  group_by(score_q) %>%
  summarise(
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  )

quant_ranges

p_bar <- ggplot(quant_df, aes(x = factor(score_q, levels = 4:1), y = pct_invasive)) +
  geom_col(width = 0.8, color = "#9B51E0", fill = "#9B51E0" ) +
  scale_y_continuous(limits = c(0, 42), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 10, colour = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.25),
    axis.ticks.y = element_line(size = 0.1, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(0.5, "mm"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )

p_bar


ggsave(
  filename = "your_path/DCIS_risk_score_histo.pdf",
  plot     = p_bar,
  device   = "pdf",
  width    = 30/25.4,
  height   = 37/25.4
)


#### Figure 4: Panel B ----

gsea_res <- read.csv("your_path/gsea_res.csv")

gsea_res <- gsea_res %>%
  mutate(FDRQval_adj = ifelse(`FDRQval` == 0, 1e-3, `FDRQval`))

p <- ggplot(gsea_res, aes(x = NES, y = reorder(NAME, NES))) +
  geom_point(aes(size = SIZE, fill = FDRQval_adj), shape = 21, color = "black", stroke = 0.25, show.legend = TRUE) +
  scale_fill_gradientn(
    colours = c("black", "red", "orange", "yellow", "white"),
    trans = "log10"
  ) +
  scale_size_continuous(range = c(2, 6)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  labs(x = NULL, y = NULL) +
  geom_vline(xintercept = 0, color = "black", size = 0.25) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.title.x = element_blank(),
    axis.ticks = element_line(size = 0.25),
    axis.ticks.length = unit(1.5, "pt"),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

p <- p +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
p
ggsave(
  filename = "your_path/gsea_plot.png",
  plot = p,
  width = 35,
  height = 75,
  units = "mm",
  dpi = 1200,
  bg = "transparent"
)

