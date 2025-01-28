# 20250121
# Authors - Benjamin Oberlton, Hadeesha Piyadasa

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


rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
set.seed(42)

group_to_classify = "survival_status" # toggle between survival_status, who_grade and idh_status

if (group_to_classify == "who_grade") {
  correlation_threshold <- 0.6
} else if (group_to_classify == "idh_status") {
  correlation_threshold <- 0.6
} else if (group_to_classify == "survival_status") {
  correlation_threshold <- 0.72
} else {
  stop("Invalid value for group_to_classify. Must be 'who_grade', 'idh_status', or 'survival_status'.")
}

# Import data -------------------------------------------------------------

modality = "total"
matrix = "nanostring"

remove_correlated_features = TRUE

hold_class_specific_features = FALSE
class_specific_cutoff = 0.7

top_feature_run = FALSE
top_feature_run_name = "total_features_no_other"

# Graph settings
# TRUE == Top number of features
# FALSE == Top percent of features
top_number_features_graphs = TRUE
top_n <- 50
top_percentage <- 0.01

run_pca_on_rf <- TRUE
top_n_pca <- 75
top_n_pca_bargraph <- 75

save_top_features = FALSE
# True is number False is percent
top_number_features_save= TRUE
number_top_features_to_save = 50
percent_top_features_to_save = 0.01
top_feature_save_name = "top_50_features_core_only"

save_plots = FALSE
# save_name = paste0("nanostring_", modality)
save_name = "total_features_total_RNA"

save_AUC_values = FALSE
AUC_append_path = "total_AUC_values_filtered_RNA"

use_selected_genes = FALSE




if(top_feature_run) {
  all_modalities_top_features <- fread(sprintf("your_path/Random_Forest_Model/%s.csv", top_feature_run_name))
}

if (group_to_classify %in% c("who_grade", "idh_status")) {
  file_path <- "20241124_nanostring_cor_matrix_who_grade.rds" # generated below
} else if (group_to_classify == "survival_status") {
  file_path <- "20241124_nanostring_cor_matrix_survival_status.rds" # generated below
} else {
  stop("Invalid value for group_to_classify. Must be 'who_grade', 'survival_status', or 'idh_status'.")
}

# Check if the file exists and load it if found
if (file.exists(file_path)) {
  cor_matrix <- readRDS(file_path)
  message(sprintf("Successfully loaded correlation matrix from %s.", file_path))
} else {
  stop(sprintf("File not found: %s. Please check the file path or classification type.", file_path))
}

# Import correlation matrix for data being analyzed
metadata <- fread("your_path/metadata_complete.csv") # download from www.bruce.parkerici.org

if (group_to_classify == "survival_status") {
  metadata_mod <- metadata %>%
    filter((!is.na(survival_duration_days_first_diag)) &
             (immunotherapy == "no") &
             (final_diagnosis == "GBM") &
             (recurrence == "no") &
             (alignment_status != 2) &
             (age >= 18) &
             (tumor_region != "other")) %>%
    mutate(death = 1,
           survival_duration_years_first_diag = survival_duration_days_first_diag / 365.25,
           sex = ifelse(sex == "M", 0, ifelse(sex == "F", 1, NA)),
           tumor_region = recode(tumor_region,
                                 "tumor_core" = 0,
                                 "tumor_core_to_infiltrating" = 1,
                                 "tumor_infiltrating" = 2),
           survival_status = ifelse(survival_duration_years_first_diag > 1.5, "long", "short"))
} else if (group_to_classify %in% c("who_grade", "idh_status")) {
  metadata_mod <- metadata %>%
    # Step 1: Initial Filtering
    filter(
      site == "Stanford",
      immunotherapy == "no",
      who_grade %in% c("2", "3", "4"),
      recurrence == "no",
      age >= 18,
      # Conditional Filtering Based on `who_grade` and `tumor_region`
      (
        (who_grade %in% c("2", "3") & tumor_region %in% c("other", "tumor_core")) |
          (who_grade == "4" & tumor_region == "tumor_core")
      )
    ) %>%
    # Step 2: Mutate and Recode
    mutate(
      death = 1,
      survival_duration_years_first_diag = survival_duration_days_first_diag / 365.25,
      sex = case_when(
        sex == "M" ~ 0,
        sex == "F" ~ 1,
        TRUE ~ NA_real_
      ),
      # Recode `tumor_region` for Selected Samples
      tumor_region = recode(tumor_region,
                            "tumor_core" = 0,
                            "other" = 1),
      survival_status = ifelse(survival_duration_years_first_diag > 1.5, "long", "short")
    )
} else {
  stop("Invalid value for group_to_classify. Must be 'survival_status', 'who_grade', or 'idh_status'.")
}

feature_table <- readRDS("your_path/20240810_master_feature_table_na_removed.rds") # download from www.bruce.parkerici.org

filtered_feature_table <- feature_table 

# Apply modality-specific filtering
filtered_feature_table <- if (modality == "nanostring") {
  filtered_feature_table %>% filter(bio_feature_type == "spatial_RNA")
} else if (modality == "maldi") {
  filtered_feature_table %>% filter(bio_feature_type == "glycans")
} else if (modality == "mibi") {
  filtered_feature_table %>% filter(!(bio_feature_type %in% c("spatial_RNA", "glycans"))) #%>% 
    #filter(feature_type != "Phenotype_marker_intensity")
  
} else {
  filtered_feature_table
}


sanitize_column_names <- function(names) {
  names %>%
    str_replace_all("\\+", "_plus_") %>%
    str_replace_all("-", "_minus_")
}

if (!top_feature_run) {
  # Continue with the rest of the transformations
  classifier_feature_table <- filtered_feature_table %>%
    mutate(feature_name = paste(cell_meta_cluster_final, feature_type, bio_feature_type, feature_variable, sep = "__")) %>% 
    mutate(feature_name = sanitize_column_names(feature_name)) %>% 
    semi_join(metadata_mod, by = "sample_id") %>%
    dplyr::select(-cell_meta_cluster_final, -feature_source, -feature_type, -feature_variable, -bio_feature_type, -fov, -source, -Broad_Feature_Type, -Feature_Class) %>% 
    # Adjust for sample_ids missing patient_ids
    left_join(metadata_mod %>% dplyr::select(sample_id, patient_id), by = "sample_id", suffix = c("", "_metadata")) %>%
    mutate(patient_id = coalesce(patient_id, patient_id_metadata)) %>%
    dplyr::select(-patient_id_metadata)
}

if (top_feature_run) {
  # Continue with the rest of the transformations
  classifier_feature_table <- filtered_feature_table %>%
    mutate(feature_name = paste(cell_meta_cluster_final, feature_type, bio_feature_type, feature_variable, sep = "__")) %>% 
    mutate(feature_name = sanitize_column_names(feature_name)) %>% 
    semi_join(metadata_mod, by = "sample_id") %>%
    dplyr::select(-cell_meta_cluster_final, -feature_source, -feature_type, -feature_variable, -bio_feature_type, -fov, -source,-Broad_Feature_Type, -Feature_Class) %>% 
    # Adjust for sample_ids missing patient_ids
    left_join(metadata_mod %>% dplyr::select(sample_id, patient_id), by = "sample_id", suffix = c("", "_metadata")) %>%
    mutate(patient_id = coalesce(patient_id, patient_id_metadata)) %>%
    dplyr::select(-patient_id_metadata)%>%
    filter(feature_name %in% all_modalities_top_features$Feature)
}

# Use only specified RNA genes ------------------------------------------------------------------------------

# Use specific gene sets
if (use_selected_genes) {
  
  classifier_feature_table <- classifier_feature_table %>% 
    pivot_wider(names_from = feature_name, values_from = feature_value, values_fn = list(feature_value = mean))
  
  # Read the CSV file without column names
  selected_genes <- read_csv("your_path/tme_gene_signatures.csv", col_names = FALSE) # download from www.bruce.parkerici.org
  
  # Convert dataframe to a single list and ignore NAs
  selected_genes_list <- unlist(lapply(selected_genes, function(x) x[!is.na(x)]))
  
  # Remove names
  names(selected_genes_list) <- NULL
  
  # Convert selected_genes_list to a character vector if it is not already
  selected_genes_list <- as.character(selected_genes_list)
  
  # Create a regex pattern for an exact match of the gene names
  gene_pattern <- paste0("^(", paste(selected_genes_list, collapse = "|"), ")$")
  
  # Create a logical vector for column selection
  selected_columns <- sapply(colnames(classifier_feature_table), function(col) {
    # Check if the column name contains "RNA"
    if (grepl("RNA", col)) {
      # Extract the part after the third '__'
      parts <- strsplit(col, "__")[[1]]
      if (length(parts) >= 4 && grepl(gene_pattern, parts[4])) {
        return(TRUE)  # The gene matches exactly
      } else {
        return(FALSE)  # No match found
      }
    } else {
      # If it doesn't contain "RNA", include the column
      return(TRUE)
    }
  })
  
  # Subset the dataframe based on the selected columns
  classifier_feature_table <- classifier_feature_table[, selected_columns]
  
  classifier_feature_table <- classifier_feature_table %>%
    pivot_longer(cols = -c(sample_id, patient_id),
                 names_to = "feature_name",
                 values_to = "feature_value")
}



length(unique(classifier_feature_table$feature_name))

# Merge in relevant metadata ----------------------------------------------
if (!hold_class_specific_features) {


  metadata_mod_longer <- metadata_mod %>%
    pivot_longer(cols = c(age, sex, tumor_region), names_to = "feature_name", values_to = "feature_value") %>%
    dplyr::select(sample_id, patient_id, feature_name, feature_value)

  classifier_feature_table_with_metadata <- bind_rows(classifier_feature_table, metadata_mod_longer)

  print(head(classifier_feature_table_with_metadata))
  # remove columns with only 0 and > 50% of values are NA, then impute NA with medians, before returning to longer format
  classifier_feature_table_with_metadata <- classifier_feature_table_with_metadata %>%
    dplyr::rename(feature_value_raw = feature_value) %>%
    pivot_wider(names_from = feature_name, values_from = feature_value_raw, values_fn = list(feature_value_raw = mean)) %>%
    dplyr::select(where(~ mean(is.na(.)) <= 0.5)) %>%
    dplyr::select(where(~ !(all(. == 0, na.rm = TRUE)))) %>%
    # mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    pivot_longer(
      cols = -c(sample_id, patient_id),  # Exclude identifier columns from pivoting
      names_to = "feature_name",
      values_to = "feature_value_raw"
    )



  length(unique(classifier_feature_table_with_metadata$feature_name))
}

# Merge in relevant metadata (hold class_specific_features) ----------------------------------------------
if (hold_class_specific_features) {
  # Transform survival_status to numeric
  metadata_mod <- metadata_mod %>%
    mutate(survival_status_numeric = ifelse(survival_status == "short", 0, 
                                            ifelse(survival_status == "long", 1, NA)))
  
  metadata_mod_longer <- metadata_mod %>% 
    pivot_longer(cols = c(age, sex, tumor_region, survival_status_numeric), names_to = "feature_name", values_to = "feature_value") %>% 
    dplyr::select(sample_id, patient_id, feature_name, feature_value) 
  
  classifier_feature_table_with_metadata <- bind_rows(classifier_feature_table, metadata_mod_longer)
  
  # Rename and pivot wider
  filter_70_data <- classifier_feature_table_with_metadata %>% 
    dplyr::rename(feature_value_raw = feature_value) %>%
    pivot_wider(names_from = feature_name, values_from = feature_value_raw, values_fn = list(feature_value_raw = mean))
  
  # Filter features that are >=50% NA
  features_50_percent_na <- filter_70_data %>% 
    dplyr::select(sample_id, patient_id, survival_status_numeric, where(~ mean(is.na(.)) >= class_specific_cutoff))
  
  # Check if >=70% of either category in survival_status is made up of non-NA samples
  features_to_keep_70_percent <- features_50_percent_na %>%
    pivot_longer(-c(sample_id, patient_id, survival_status_numeric), names_to = "feature_name", values_to = "feature_value") %>%
    group_by(feature_name, survival_status_numeric) %>%
    summarise(non_na_ratio = mean(!is.na(feature_value)), .groups = 'drop') %>%
    filter(non_na_ratio >= 0.5) %>%
    dplyr::select(feature_name) %>%
    distinct()
  
  features_to_keep_70_percent_list <- features_to_keep_70_percent$feature_name
  
  # Remove columns with only 0 and > 50% of values are NA, then impute NA with medians, before returning to longer format
  classifier_feature_table_with_metadata <- classifier_feature_table_with_metadata %>% 
    dplyr::rename(feature_value_raw = feature_value) %>%
    pivot_wider(names_from = feature_name, values_from = feature_value_raw, values_fn = list(feature_value_raw = mean)) %>% 
    dplyr::select(intersect(colnames(.), features_to_keep_70_percent_list) | where(~ mean(is.na(.)) <= 0.5)) %>%
    dplyr::select(where(~ !(all(. == 0, na.rm = TRUE)))) %>% 
    pivot_longer(
      cols = -c(sample_id, patient_id, survival_status_numeric),  # Exclude identifier columns from pivoting
      names_to = "feature_name",
      values_to = "feature_value_raw"
    ) %>% 
    dplyr::select(-survival_status_numeric)
  
  
  
  length(unique(classifier_feature_table_with_metadata$feature_name))
}
# Kaplan-Meier Curve ------------------------------------------------------
# survival data object
surv_obj <- Surv(time = metadata_mod$survival_duration_years_first_diag, event = metadata_mod$death)

# fit kaplan-meier curve
s1 <- survfit(surv_obj ~ 1, data = metadata_mod)

# plot kaplan-meier curve
km_curve <- ggsurvfit(s1, linewidth = 1) +
  labs(x = "Years", y = "overall survival") +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() 

km_curve


# Adjust tumor region ---------------------------------------------
classifier_feature_table_with_metadata_wide <- classifier_feature_table_with_metadata %>% 
  pivot_wider(names_from = feature_name, values_from = feature_value_raw)

if (group_to_classify == "survival_status") {
  tumor_region_percentages <- classifier_feature_table_with_metadata_wide %>%
    group_by(patient_id, tumor_region) %>%
    summarize(region_count = n(), .groups = "drop") %>%
    mutate(total_count = sum(region_count)) %>%
    ungroup() %>%
    mutate(percentage = region_count / total_count) %>%
    dplyr::select(patient_id, tumor_region, percentage) %>%
    pivot_wider(
      names_from = tumor_region,
      values_from = percentage,
      values_fill = 0,
      names_prefix = "tumor_"
    ) %>%
    dplyr::rename(
      tumor_core = tumor_0,
      tumor_core_to_infiltrating = tumor_1,
      tumor_infiltrating = tumor_2
    )
} else if (group_to_classify %in% c("who_grade", "idh_status")) {
  tumor_region_percentages <- classifier_feature_table_with_metadata_wide %>%
    group_by(patient_id, tumor_region) %>%
    summarize(region_count = n(), .groups = "drop") %>%
    mutate(total_count = sum(region_count)) %>%
    ungroup() %>%
    mutate(percentage = region_count / total_count) %>%
    dplyr::select(patient_id, tumor_region, percentage) %>%
    pivot_wider(
      names_from = tumor_region,
      values_from = percentage,
      values_fill = 0,
      names_prefix = "tumor_"
    ) %>%
    dplyr::rename(
      tumor_core = tumor_0,
      other = tumor_1
    )
} else {
  stop("Invalid value for group_to_classify. Must be 'survival_status', 'who_grade', or 'idh_status'.")
}

classifier_feature_table_with_metadata_wide <- classifier_feature_table_with_metadata_wide %>%
  left_join(tumor_region_percentages, by = "patient_id") %>% 
  dplyr::select(-tumor_region)


classifier_feature_table_with_metadata <- classifier_feature_table_with_metadata_wide %>%
  pivot_longer(cols = -c(sample_id, patient_id),
               names_to = "feature_name",
               values_to = "feature_value_raw")

# Normalize values ---------------------------------------------
classifier_feature_table_with_metadata_wide <- classifier_feature_table_with_metadata %>% 
  pivot_wider(names_from = feature_name, values_from = feature_value_raw)

# Define the columns you want to exclude
exclude_cols <- c("sample_id", "patient_id", "sex", 
                  "tumor_core", "tumor_core_to_infiltrating", 
                  "tumor_infiltrating","other")

# Select continuous columns by excluding the specified columns if they exist
continuous_columns <- classifier_feature_table_with_metadata_wide %>%
  dplyr::select(-any_of(exclude_cols)) %>%  # Exclude only existing columns
  names()                            # Retrieve the column names

classifier_feature_table_with_metadata_wide_normalized <- classifier_feature_table_with_metadata_wide %>%
  mutate(across(all_of(continuous_columns), ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

classifier_feature_table_with_metadata_long_normalized <- classifier_feature_table_with_metadata_wide_normalized %>%
  pivot_longer(
    cols = -c(sample_id, patient_id),  # Exclude identifier columns from pivoting
    names_to = "feature_name",
    values_to = "feature_value_normalized"
  )

# get average patient data
normalized_classifier_feature_table <- classifier_feature_table_with_metadata %>%
  left_join(classifier_feature_table_with_metadata_long_normalized, 
            by = c("sample_id", "patient_id", "feature_name"))  %>% 
  group_by(patient_id, feature_name) %>%
  summarize(mean_feature_value = mean(feature_value_normalized, na.rm = TRUE)) %>%
  ungroup()

# Calculate Correlations and Remove correlated features ------------------------------------------


# Read the CSV file without column names
selected_genes <- read_csv("your_path/tme_gene_signatures.csv", col_names = TRUE) # download from www.bruce.parkerici.org

# Assuming `selected_genes` is your dataframe
selected_genes <- selected_genes %>%
  rename_with(~ gsub("-", "_", .x)) %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  rename_with(~ gsub("[^A-Za-z0-9_]", "", .x))

# Display the modified column names
colnames(selected_genes)

# Convert all values in the dataframe to a single character vector and remove NAs
selected_genes_list <- unlist(selected_genes, use.names = FALSE)
selected_genes_list <- selected_genes_list[!is.na(selected_genes_list)] %>% as.character()

# Display the list
print(selected_genes_list)
selected_genes_list <- sort(selected_genes_list)

# Combine selected_genes_list into a single regex pattern for an exact match
gene_pattern <- paste0("^(", paste(selected_genes_list, collapse = "|"), ")$")

# Split the feature_name column once and create a new column for the fourth part
normalized_classifier_feature_table$part4 <- sapply(strsplit(normalized_classifier_feature_table$feature_name, "__"), function(parts) {
  if (length(parts) >= 4) parts[4] else NA
})

selected_genes_features <- unique(
  normalized_classifier_feature_table$feature_name[
    normalized_classifier_feature_table$part4 %in% selected_genes_list & 
      grepl("spatial_RNA", normalized_classifier_feature_table$feature_name)
  ]
)

# Optionally remove the new column
normalized_classifier_feature_table$part4 <- NULL

# Print the result
print(selected_genes_features)

if (remove_correlated_features) {
  # Step 1: Calculate the Correlation Matrix
  # Pivot the data to wide format and set to patient_data 
  normalized_wide <- normalized_classifier_feature_table %>%
    dplyr::select(patient_id, feature_name, mean_feature_value) %>%
    pivot_wider(names_from = feature_name, values_from = mean_feature_value) %>%
    dplyr::select(-patient_id) 
  
  # Identify features with zero standard deviation
  constant_features <- sapply(normalized_wide, function(x) sd(x, na.rm = TRUE) == 0)
  
  # Handle NA standard deviations separately
  na_features <- is.na(constant_features)
  
  # Remove constant features and NA standard deviation features from the dataset
  normalized_wide_filtered <- normalized_wide[, !(constant_features | na_features)]
  
  # #Calculate the correlation matrix
  # #Takes a long time to run, only needs to be run once to get the correlation matrix
  # # Replace pbapply with apply
  # cor_matrix <- apply(normalized_wide_filtered, 2, function(x) {
  #   sapply(normalized_wide_filtered, function(y) cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  # })
  # 
  # # # Calculate the correlation matrix with a progress bar
  # # cor_matrix <- pbapply::pbapply(normalized_wide_filtered, 2, function(x) {
  # #   sapply(normalized_wide_filtered, function(y) cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  # # })
  # 
  # # Save the correlation matrix
  # saveRDS(cor_matrix, sprintf("20241124_nanostring_cor_matrix_%s.rds", group_to_classify))

  # Step 2: Identify Highly Correlated Features

  
  # List of features to keep
  features_to_keep <- c("sex", "age", "tumor_core", "tumor_core_to_infiltrating", "tumor_infiltrating","other")
  if (modality == "total") {
    additional_features_to_keep <- unique(normalized_classifier_feature_table$feature_name[!grepl("spatial_RNA", normalized_classifier_feature_table$feature_name)])
    features_to_keep <- union(features_to_keep, additional_features_to_keep)
  }

  # Add features that contain any value from selected_genes_list

  features_to_keep <- union(features_to_keep, selected_genes_features)

  # Define a custom function to handle NAs in correlation matrix
  findCorrelationIgnoreNA <- function(cor_matrix, cutoff) {
    # Replace NA values with a value outside the correlation range (0 here)
    cor_matrix_no_na <- ifelse(is.na(cor_matrix), 0, cor_matrix)
    
    # Find highly correlated pairs using findCorrelation function from caret package
    highly_correlated <- findCorrelation(cor_matrix_no_na, cutoff = cutoff, names = TRUE)
    
    return(highly_correlated)
  }
  
  # Find pairs of features with correlation above the threshold
  highly_correlated <- findCorrelationIgnoreNA(cor_matrix, cutoff = correlation_threshold)
  
  # Exclude the features to keep from the list of highly correlated features
  highly_correlated <- setdiff(highly_correlated, features_to_keep)
  
  # Step 3: Remove Redundant Features
  # Remove one feature from each pair of highly correlated features
  selected_features <- setdiff(colnames(cor_matrix), highly_correlated)
  
  # Include the features to keep in the final selected features
  selected_features <- union(selected_features, features_to_keep)
  
  # Filter the original data to retain only the selected features
  final_normalized_classifier_feature_table <- normalized_classifier_feature_table %>%
    filter(feature_name %in% selected_features)

  
  
  # Check the resulting table
  print(final_normalized_classifier_feature_table)
  print(length(unique(final_normalized_classifier_feature_table$feature_name)))
}  else {
  final_normalized_classifier_feature_table <- normalized_classifier_feature_table 
  # remove 0 std dev columns
}


# Random Forest Classifier data setup---------------------------------------
# Ensure the feature table has one row per patient and the features as columns
feature_data <- final_normalized_classifier_feature_table %>%
  pivot_wider(names_from = feature_name, values_from = mean_feature_value) %>% 
  rename_with(~ make.names(.)) %>% 
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) 

feature_metadata <- metadata_mod %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  dplyr::select(patient_id, !!sym(group_to_classify), survival_duration_days_first_diag) # Replaces 'who_grade' with 'group_to_classify' value dynamically


# Merge feature data with survival status
data_for_model <- feature_data %>%
  left_join(feature_metadata, by = "patient_id")

# Remove patient_id column for modeling
data_for_model <- data_for_model %>%
  dplyr::select(-patient_id)

if (group_to_classify %in% c("who_grade", "idh_status")) {
  data_for_model <- data_for_model %>%
    dplyr::select(-any_of(c("tumor_core", "tumor_infiltrating", "tumor_core_to_infiltrating", "other", "survival_duration_days_first_diag")))
  
  data_for_model[[group_to_classify]] <- ifelse(
    grepl("^[0-9]+$", data_for_model[[group_to_classify]]),
    paste0("X", data_for_model[[group_to_classify]]),
    data_for_model[[group_to_classify]]
  )
} else if (group_to_classify == "survival_status") {
  data_for_model <- data_for_model %>%
    dplyr::select(-any_of("survival_duration_days_first_diag"))
}

# Integrated workflow ----

# Function to run the model and return AUC values and importance values

run_model <- function(data, label_column, num_runs = 10, randomize_labels = FALSE) {
  set.seed(123)
  auc_values <- numeric(num_runs)
  rf_gini_importance_list <- list()  # Gini importance (MeanDecreaseGini) from rf_model$finalModel
  oob_importance_list <- list()      # OOB importance (MeanDecreaseAccuracy)
  varimp_gini_importance_list <- list()  # Gini importance from varImp function
  
  for (i in 1:num_runs) {
    print(paste("Starting run", i, "of", num_runs))
    
    # Optional randomization of labels
    if (randomize_labels) {
      data[[label_column]] <- sample(data[[label_column]])
    }
    
    stratified_folds <- createFolds(data[[label_column]], k = 3, list = TRUE, returnTrain = TRUE)
    
    # Define trainControl with stratified folds
    train_control <- trainControl(method = "cv", 
                                  number = 3, 
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = TRUE,
                                  summaryFunction = multiClassSummary,
                                  index = stratified_folds) # Use stratified folds
    
    # Train the random forest model with ROC as the metric
    print("Training Random Forest model...")
    rf_model <- train(as.formula(paste(label_column, "~ .")), data = data,
                      method = "rf",
                      metric = "ROC",
                      trControl = train_control,
                      importance = TRUE)
    
    # Evaluate the model and extract AUC for the run
    x <- evalm(rf_model, plots = "r")
    auc_values[i] <- x$stdres$`Group 1`["AUC-ROC", "Score"]
    
    # Get Gini importance (MeanDecreaseGini) from rf_model
    print("Extracting Gini importance from rf_model (Mean Decrease in Gini)...")
    rf_gini_importance <- rf_model$finalModel$importance[, "MeanDecreaseGini", drop = FALSE]
    rf_gini_importance_list[[i]] <- rf_gini_importance
    
    # Get OOB importance (Mean Decrease in Accuracy) from rf_model
    print("Extracting OOB importance from rf_model (Mean Decrease in Accuracy)...")
    oob_importance <- rf_model$finalModel$importance[, "MeanDecreaseAccuracy", drop = FALSE]
    oob_importance_list[[i]] <- oob_importance
    
    # Get Gini importance from varImp function
    print("Extracting Gini importance from varImp function...")
    varimp_gini_importance <- varImp(rf_model, scale = FALSE)$importance[, 1, drop = FALSE]
    varimp_gini_importance_list[[i]] <- varimp_gini_importance
  }
  
  # Combine Gini importance from rf_model into a single data frame
  print("Combining Gini importance from rf_model (MeanDecreaseGini) from all runs...")
  rf_gini_importance_df <- do.call(cbind, rf_gini_importance_list)
  colnames(rf_gini_importance_df) <- paste0("Run_", 1:num_runs)
  rf_gini_importance_df <- as.data.frame(rf_gini_importance_df)
  rf_gini_importance_df$Median <- apply(rf_gini_importance_df, 1, median)
  rf_gini_importance_df$Feature <- rownames(rf_gini_importance_df)
  rownames(rf_gini_importance_df) <- NULL
  rf_gini_importance_df <- rf_gini_importance_df[, c("Feature", "Median", paste0("Run_", 1:num_runs))]
  
  # Combine OOB importance into a single data frame
  print("Combining OOB importance (MeanDecreaseAccuracy) from all runs...")
  oob_importance_df <- do.call(cbind, oob_importance_list)
  colnames(oob_importance_df) <- paste0("Run_", 1:num_runs)
  oob_importance_df <- as.data.frame(oob_importance_df)
  oob_importance_df$Median <- apply(oob_importance_df, 1, median)
  oob_importance_df$Feature <- rownames(oob_importance_list[[1]])
  rownames(oob_importance_df) <- NULL
  oob_importance_df <- oob_importance_df[, c("Feature", "Median", paste0("Run_", 1:num_runs))]
  
  # Combine Gini importance from varImp into a single data frame
  print("Combining Gini importance from varImp function from all runs...")
  varimp_gini_importance_df <- do.call(cbind, varimp_gini_importance_list)
  colnames(varimp_gini_importance_df) <- paste0("Run_", 1:num_runs)
  varimp_gini_importance_df <- as.data.frame(varimp_gini_importance_df)
  varimp_gini_importance_df$Median <- apply(varimp_gini_importance_df, 1, median)
  varimp_gini_importance_df$Feature <- rownames(varimp_gini_importance_df)
  rownames(varimp_gini_importance_df) <- NULL
  varimp_gini_importance_df <- varimp_gini_importance_df[, c("Feature", "Median", paste0("Run_", 1:num_runs))]
  
  print("Returning AUC values, Gini importance (rf_model and varImp), and OOB importance dataframes...")
  
  # Return the results
  list(auc_values = auc_values, 
       rf_gini_importance_df = rf_gini_importance_df,  # Gini from rf_model
       oob_importance_df = oob_importance_df,          # OOB importance
       varimp_gini_importance_df = varimp_gini_importance_df)  # Gini from varImp
}

original_results <- run_model(data_for_model, group_to_classify)

# Randomize the labels and run the model again
randomized_results <- run_model(data_for_model, group_to_classify, randomize_labels = TRUE)

# Combine AUC values into a single data frame
auc_data <- data.frame(
  AUC = c(original_results$auc_values, randomized_results$auc_values),
  Label = rep(c("Original", "Randomized"), each = 10)
)

file_path <- sprintf("../../tables/auc_data_%s.csv", group_to_classify)

if (file.exists(file_path)) {
  response <- readline(prompt = sprintf("The file '%s' already exists. Do you want to overwrite it? (yes/no): ", file_path))
  if (tolower(response) == "no") {
    message("File not overwritten.")
  } else {
    write_csv(auc_data, file_path)
    message(sprintf("File '%s' has been overwritten.", file_path))
  }
} else {
  write_csv(auc_data, file_path)
  message(sprintf("File '%s' has been saved.", file_path))
}

# Perform a t-test to compare the two groups
t_test_result <- t.test(AUC ~ Label, data = auc_data)
p_value <- t_test_result$p.value

# Create a combined box plot for comparison with custom aesthetics
auc_boxplot <- ggplot(auc_data, aes(x = Label, y = AUC, fill = Label)) +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.7) + # Hide outliers and set black outline
  scale_fill_manual(values = c("Original" = "#2E86C1", "Randomized" = "#E74C3C")) + # Custom fill colors
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "black") + # Add individual points with transparency
  xlab("Label Type") +
  ylab("AUC") +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "none" # Hide legend if not needed
  ) +
  annotate("text", x = 1.5, y = max(auc_data$AUC), label = paste("p-value:", signif(p_value, digits = 3)), size = 8, hjust = 0.5, vjust = -0.5)

print(auc_boxplot)

rf_gini_importance_df <- original_results$rf_gini_importance_df   # Gini importance from rf_model


# Create a new column based on the presence of specific terms in the Feature column
rf_gini_importance_df$Feature_Type <- ifelse(grepl("glycans", rf_gini_importance_df$Feature, ignore.case = TRUE), "Glycan",
                                     ifelse(grepl("RNA", rf_gini_importance_df$Feature, ignore.case = TRUE), "RNA", "MIBI"))


rf_gini_importance_df <- rf_gini_importance_df %>%
  dplyr::rename(Importance = Median) %>%
  arrange(desc(Importance)) 

extract_after_third_double_underscore <- function(feature) {
  parts <- unlist(strsplit(feature, "__"))
  if (grepl("intensity", feature)) {
    if (length(parts) >= 3) {
      return(paste(c(parts[1], parts[3:length(parts)]), collapse = "__"))
    } else {
      return("")
    }
  } else {
    if (length(parts) > 3) {
      return(paste(parts[4:length(parts)], collapse = "__"))
    } else {
      return("")
    }
  }
}

# Create the new column Feature_type_two
rf_gini_importance_df$Feature_type_two <- sapply(rf_gini_importance_df$Feature, extract_after_third_double_underscore)

# Combine Feature_Type and Feature_Type_two
rf_gini_importance_df <- rf_gini_importance_df %>%
  mutate(Combined_Feature = paste(Feature_Type, Feature_type_two, sep = "_")) %>%
  arrange(desc(Importance)) 

selected_genes_list <- c("ACTA2", "ADAMTS4", "ADAMTS5", "ANGPT1", "ANGPT2", "ARG1", "AURKA", "AURKB",
                         "B2M", "BLK", "BTLA", "BUB1", "CA9", "CCL1", "CCL11", "CCL15",
                         "CCL17", "CCL2", "CCL22", "CCL26", "CCL28", "CCL3", "CCL3", "CCL4",
                         "CCL5", "CCL7", "CCL8", "CCNB1", "CCND1", "CCNE1", "CCR10", "CCR2",
                         "CCR3", "CCR4", "CCR4", "CCR8", "CCR8", "CD160", "CD163", "CD177",
                         "CD19", "CD22", "CD226", "CD244", "CD248", "CD27", "CD274", "CD28",
                         "CD28", "CD3D", "CD3E", "CD3G", "CD40", "CD40LG", "CD40LG", "CD5",
                         "CD68", "CD70", "CD79A", "CD79B", "CD80", "CD83", "CD86", "CD8A",
                         "CD8B", "CDH2", "CDH5", "CDH5", "CDK2", "CETN3", "CIITA", "CLEC14A",
                         "CMKLR1", "COL11A1", "COL1A1", "COL1A1", "COL1A2", "COL1A2", "COL3A1", "COL4A1",
                         "COL5A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "CR2", "CSF1", "CSF1",
                         "CSF1R", "CSF1R", "CSF1R", "CSF2", "CSF2RA", "CSF3", "CSF3R", "CTLA4",
                         "CTLA4", "CTSG", "CX3CL1", "CX3CR1", "CXCL1", "CXCL10", "CXCL11", "CXCL12",
                         "CXCL12", "CXCL2", "CXCL5", "CXCL5", "CXCL5", "CXCL8", "CXCL8", "CXCL8",
                         "CXCL9", "CXCR1", "CXCR1", "CXCR2", "CXCR2", "CXCR2", "CXCR2", "CXCR3",
                         "CXCR4", "CYBB", "E2F1", "ELANE", "ELN", "ENG", "EOMES", "EOMES",
                         "ESCO2", "FAP", "FASLG", "FBLN1", "FCGR3B", "FCRL5", "FFAR2", "FGF2",
                         "FGFBP2", "FLT1", "FLT1", "FN1", "FOXP3", "GNLY", "GNLY", "GZMA",
                         "GZMB", "GZMB", "GZMH", "GZMK", "HAVCR2", "HLA-A", "HLA-B", "HLA-C",
                         "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
                         "ICOS", "ICOSLG", "IDO1", "IFNA2", "IFNB1", "IFNG", "IFNG", "IFNG",
                         "IKZF2", "IKZF4", "IL10", "IL10", "IL10", "IL10", "IL10", "IL12A",
                         "IL12B", "IL12RB2", "IL13", "IL1B", "IL2", "IL21", "IL21", "IL22",
                         "IL23A", "IL4", "IL4I1", "IL4I1", "IL5", "IL6", "IL6", "IL6",
                         "IL6R", "IRF5", "ITK", "KDR", "KDR", "KIR2DL4", "KITLG", "KLRC2",
                         "KLRF1", "KLRK1", "LAG3", "LAMA3", "LAMB3", "LAMC2", "LGALS7", "LGALS9",
                         "LOX", "LRP1", "LUM", "MCM2", "MCM6", "MFAP5", "MIF", "MKI67",
                         "MMP1", "MMP11", "MMP12", "MMP2", "MMP2", "MMP3", "MMP3", "MMP7",
                         "MMP9", "MMRN1", "MMRN2", "MPO", "MRC1", "MS4A1", "MSR1", "MYBL2",
                         "NCR1", "NCR3", "NKG7", "NLRC5", "NOS2", "NOS3", "PAX5", "PDCD1",
                         "PDCD1LG2", "PDGFC", "PDGFRA", "PDGFRB", "PGF", "PGLYRP1", "PLK1", "PLOD2",
                         "PRF1", "PRTN3", "PTGS2", "SH2D1B", "SIGLEC1", "SNAI1", "SNAI2", "SOCS3",
                         "STAP1", "STAT4", "TAP1", "TAP2", "TAPBP", "TBX21", "TBX21", "TBX21",
                         "TEK", "TGFB1", "TGFB2", "TGFB3", "TIGIT", "TNC", "TNF", "TNF",
                         "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF10", "TNFSF4",
                         "TNFSF9", "TRAC", "TRAT1", "TRBC1", "TRBC2", "TWIST1", "TWIST2", "VCAM1",
                         "VEGFA", "VEGFB", "VEGFC", "VSIR", "VTN", "VWF", "VWF", "XCL1",
                         "XCR1", "ZAP70", "ZEB1", "ZEB2")


filtered_df <- rf_gini_importance_df %>%
  filter(Feature_type_two %in% selected_genes_list)
# 
# 
# 
# Define the lists of terms for each feature type
tumor_antigen_terms <- c("B7H3", "EGFR", "GM2_GD2", "GPC2", "HER2", "NG2", "VISTA","Tumor")
immune_cell_terms <- c("Immune_unassigned", "Macrophage_CD206", "Macrophage_CD68",
                       "Macrophage_CD68_CD163", "Unassigned", "Myeloid_CD14",
                       "Myeloid_CD11b_HLADRminus", "Tcell_CD8", "Endothelial_cells",
                       "Myeloid_CD11b_HLADRplus", "DC_Mac_CD209", "Microglia",
                       "Tcell_CD4", "Myeloid_CD141", "Tcell_FoxP3", "APC",
                       "Microglia_CD163", "Myeloid_CD14_CD163", "Neurons",
                       "Neutrophils", "Bcells", "Mast_cells")

# Initialize the new column with NA
rf_gini_importance_df$Features_broad <- NA

rna_features_df <- rf_gini_importance_df %>%
  filter(Feature_Type == "RNA")

rna_features_df <- rna_features_df %>%
  mutate(Gene = Feature_type_two)

rna_features_filtered <- rna_features_df %>%
  filter(Importance > 0)

gene_list <- unique(rna_features_filtered$Gene)

#
# # # Assuming rna_features_filtered is your filtered dataframe
# gene_list <- rna_features_filtered$Importance
# names(gene_list) <- rna_features_filtered$Gene
# #
# # # Sort in decreasing order
# gene_list <- sort(gene_list, decreasing = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
# Assuming you have the Entrez IDs or need conversion, you may need an additional step to map gene symbols to Entrez IDs
# gene_list can be directly used in the GSEA analysis

#Assuming gene_list is a vector of gene symbols
ora_results <- enrichGO(gene         = gene_list,          # Your list of genes
                        OrgDb        = org.Hs.eg.db,       # Organism database (human)
                        keyType      = "SYMBOL",           # Use gene symbols as input
                        ont          = "ALL",              # Ontology: "BP", "CC", "MF", or "ALL"
                        pAdjustMethod = "BH",              # Adjust p-values using Benjamini-Hochberg
                        pvalueCutoff = 0.1,               # P-value cutoff
                        qvalueCutoff = 0.1)               # q-value cutoff

# View the results
head(ora_results)
#
# # Visualize the results
dotplot(ora_results)
# 
# 


# Assuming `selected_genes` is your dataframe
selected_genes <- selected_genes %>%
  rename_with(~ gsub("-", "_", .x)) %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  rename_with(~ gsub("[^A-Za-z0-9_]", "", .x))

# Display the modified column names
colnames(selected_genes)



# Flatten the selected_genes list while keeping track of the column names
selected_genes_mapping <- lapply(names(selected_genes), function(col) {
  genes <- selected_genes[[col]]
  data.frame(gene = genes, column = col, stringsAsFactors = FALSE)
})

# Combine all the data frames into one
selected_genes_mapping <- do.call(rbind, selected_genes_mapping)

# Remove rows with NA genes
selected_genes_mapping <- selected_genes_mapping[!is.na(selected_genes_mapping$gene), ]

# Create a list where each gene maps to a vector of columns
selected_genes_dict <- split(selected_genes_mapping$column, selected_genes_mapping$gene)

# Display the resulting dictionary
selected_genes_dict
# 
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

classes <- read_csv("../../tables/glycan_classes.csv")


# Capitalize the first letter of every column name
colnames(classes) <- str_to_title(colnames(classes))

# Rename specific columns
classes <- classes %>%
  dplyr::rename(
    `High Mannose` = `Highmannose`,
    `poly-LacNAc` = `Polylacnac`
  )

# Check the column names
colnames(classes)

# Create the mapping
glycan_classes_mapping <- classes %>%
  dplyr::select(-Mz) %>%
  pivot_longer(cols = -Composition, names_to = "class", values_to = "presence") %>%
  filter(presence == 1) %>%
  dplyr::select(Composition, class) %>%
  group_by(Composition) %>%
  summarise(classes = list(class)) %>%
  mutate(Composition = gsub("\\+", "_", Composition))

# Convert the tibble to a list for easier lookup
glycan_classes_dict <- setNames(glycan_classes_mapping$classes, glycan_classes_mapping$Composition)

# Function to assign Features_broad based on Combined_Feature
assign_features_broad <- function(combined_feature) {
  if (grepl("Endothelial_cells", combined_feature)) {
    return("Endothelial_cell_features")
  } else if (grepl("Glycan_", combined_feature)) {
    glycan <- sub("^Glycan_", "", combined_feature)
    return(paste(glycan_classes_dict[[glycan]], collapse = ", "))
  } else if (grepl("RNA_", combined_feature)) {
    gene <- sub("^RNA_", "", combined_feature)
    return(paste(selected_genes_dict[[gene]], collapse = ", "))
  } else if (grepl(paste(tumor_antigen_terms, collapse = "|"), combined_feature)) {
    return("Tumor_antigen_features")
  } else if (grepl(paste(immune_cell_terms, collapse = "|"), combined_feature) & !grepl("Endothelial_cells", combined_feature)) {
    return("Immune_cell_features")
  } else {
    return(NA)
  }
}

# Apply the function to create the Features_broad column
rf_gini_importance_df$Features_broad <- sapply(rf_gini_importance_df$Combined_Feature, assign_features_broad)


file_path <- sprintf("../../tables/rf_gini_importance_df_%s.csv", group_to_classify)

if (file.exists(file_path)) {
  response <- readline(prompt = sprintf("The file '%s' already exists. Do you want to overwrite it? (yes/no): ", file_path))
  if (tolower(response) == "no") {
    message("File not overwritten.")
  } else {
    write_csv(rf_gini_importance_df, file_path)
    message(sprintf("File '%s' has been overwritten.", file_path))
  }
} else {
  write_csv(rf_gini_importance_df, file_path)
  message(sprintf("File '%s' has been saved.", file_path))
}


# Run PCA ------------

rf_gini_importance_df <- read_csv("../../tables/rf_gini_importance_df_survival_status.csv")

# Prepare data for PCA by removing the survival_status column
set.seed(123)
if(run_pca_on_rf) {
  # Select top features for PCA
  top_features_for_pca <- rf_gini_importance_df[1:top_n_pca, ]
  # Exclude unwanted columns and include only the matching features
  data_pca <- data_for_model %>%
    dplyr::select(-survival_status) %>% 
    dplyr::select(all_of(c(top_features_for_pca$Feature))) 
    # dplyr::select(-Age.at.1st.Dx.in.Years)
  
  ncol(data_pca) 
} else {
  data_pca <- data_for_model %>%
    dplyr::select(-survival_status) %>% 
    dplyr::select(-age)
}


# Perform PCA
pca_result <- prcomp(data_pca, center = TRUE, scale. = FALSE)

# Summarize PCA results
summary(pca_result)

# Inspect the PCA loadings
pca_result$rotation

# Visualize the PCA results
# Create a data frame with the PCA results and survival_status
pca_df <- data.frame(pca_result$x, survival_status = data_for_model$survival_status)


# Calculate the variance explained by each principal component
explained_variance <- pca_result$sdev^2

# Calculate the proportion of variance explained
explained_variance_ratio <- explained_variance / sum(explained_variance)

# Create a data frame with the PC number and the variance explained
explained_variance_df <- data.frame(
  PC = seq_along(explained_variance_ratio),
  Variance_Explained = explained_variance_ratio
)

# Print the variance explained by each PC
explained_variance_df


# Plot the first two principal components with customized legend title
PCA_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = survival_status, fill = survival_status)) +
  geom_point(size = 5, alpha = 0.7, shape = 21, stroke = 0.8, color = "black")  +
  stat_ellipse(level = 0.95, linetype = 2, geom = "polygon", alpha = 0.1) +
  scale_color_manual(values = c("#2E86C1", "#E74C3C", "#28B463")) +
  scale_fill_manual(values = c("#2E86C1", "#E74C3C", "#28B463")) +
  labs(
    title = paste(toTitleCase(modality), "PCA of Survival Data"),
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Survival Status",
    fill = "Survival Status"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  theme(
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "gray")
  )

print(PCA_plot)

# Extract the required columns
pca_subset <- pca_df[, c("PC1", "PC2", "survival_status")]

# Save the subset as a CSV file
write.csv(pca_subset, "../../tables/pca_subset.csv", row.names = FALSE)




