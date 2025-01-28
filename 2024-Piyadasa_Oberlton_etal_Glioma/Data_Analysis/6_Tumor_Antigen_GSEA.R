library(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(tidyr)
library(data.table)
library(grid)
library(caret)
library(glmnet)
library(mice)
library(broom)
library(randomForest)
library(cowplot)
library(pROC)
library(MLeval)
library(tools)
library(MLmetrics)
library(fgsea)

metadata <- read.csv('/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/metadata/metadata_complete.csv')
master_feature_table_filtered <- readRDS("/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/analysis/feature_tables/20240810_master_feature_table_na_removed.rds")


# merge metadata with feature table
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
master_feature_table_filtered <- combine_data_metadata(master_feature_table_filtered, metadata, merge_by = "sample_id", remove_duplicates = FALSE)




filtered_df <- master_feature_table_filtered %>%
  filter(final_diagnosis == "GBM",
         site == "Stanford",
         tumor_region %in% c("tumor_core", "tumor_infiltrating"),
         feature_variable %in% c("B7H3_func_over_all_tumor_count_prop", "GM2_GD2_func_over_all_tumor_count_prop"))


#### GSEA ------------------------------
# Step 1: Calculate the median for each feature_variable
median_values <- filtered_df %>%
  filter(feature_variable == "GM2_GD2_func_over_all_tumor_count_prop") %>% 
  group_by(feature_variable) %>%
  summarize(median_value = median(feature_value, na.rm = TRUE))

# Step 2: Merge median values back with the original dataframe and calculate log2 fold change
log2_fc_df <- filtered_df %>%
  inner_join(median_values, by = "feature_variable") %>%
  mutate(log2_fold_change = log2(feature_value / median_value))

# Step 3: Sort the dataframe by log2_fold_change in descending order
log2_fc_df_sorted <- log2_fc_df %>%
  arrange(desc(log2_fold_change))

# Step 4: Create a ranked list of samples by log2 fold change
ranked_list <- log2_fc_df_sorted %>%
  mutate(sample_id = gsub("^Stanford-", "", sample_id)) %>%
  arrange(desc(log2_fold_change)) %>%
  select(sample_id, log2_fold_change) %>%
  deframe()  # Converts to named vector for GSEA

# Step 4.1: Remove non-finite values from the ranked list
ranked_list <- ranked_list[is.finite(ranked_list)]

# Step 5: Define samples for tumor_core and tumor_infiltrating
tumor_core_samples <- log2_fc_df_sorted %>%
  filter(tumor_region == "tumor_core") %>%
  pull(sample_id)

tumor_infiltrating_samples <- log2_fc_df_sorted %>%
  filter(tumor_region == "tumor_infiltrating") %>%
  pull(sample_id)

tumor_core_samples <- gsub("^Stanford-", "", tumor_core_samples)
tumor_infiltrating_samples <- gsub("^Stanford-", "", tumor_infiltrating_samples)


# Create a list of samples
samples <- list(
  "Tumor_Core" = tumor_core_samples,
  "Tumor_Infiltrating" = tumor_infiltrating_samples
)

# Automatically calculate minSize and maxSize based on the gene sets
minSize <- min(sapply(samples, length)) 
maxSize <- max(sapply(samples, length))  

# Step 6: Perform GSEA using the fgseaMultilevel function
fgsea_results <- fgseaMultilevel(pathways = samples, 
                                 stats = ranked_list, 
                                 minSize = minSize,   #
                                 maxSize = maxSize)   

# Step 7: Display the GSEA results
print(fgsea_results)

# Step 8: Plot enrichment score for tumor_core
# Define the save location
save_location <- '/Users/boberlto/Library/CloudStorage/GoogleDrive-boberlto@stanford.edu/Shared drives/BRUCE_data/resource_paper_drafts/Figures/Supplement/Figure_3/GM2GD2_GSEA'

if (!endsWith(save_location, "/")) {
  save_location <- paste0(save_location, "/")
}

# Step 8: Plot enrichment score for tumor_core
tumor_core_plot <- plotEnrichment(samples[["Tumor_Core"]], ranked_list) +
  ggtitle("Tumor Core Enrichment") +
  theme_minimal()

# Save the plot
ggsave(paste0(save_location, "tumor_core_enrichment.png"), tumor_core_plot, width = 6, height = 4, dpi = 300)

# Step 9: Plot enrichment score for tumor_infiltrating
tumor_infiltrating_plot <- plotEnrichment(samples[["Tumor_Infiltrating"]], ranked_list) +
  ggtitle("Tumor Infiltrating Enrichment") +
  theme_minimal()

# Save the plot
ggsave(paste0(save_location, "tumor_infiltrating_enrichment.png"), tumor_infiltrating_plot, width = 6, height = 4, dpi = 300)

# Optional: Save the GSEA results as a CSV file
fgsea_results <- fgsea_results %>%
  mutate(across(where(is.list), ~ sapply(., toString)))
write.csv(fgsea_results, paste0(save_location, "GSEA_results.csv"), row.names = FALSE)
