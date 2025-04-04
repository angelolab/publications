library(glmnet)
library(caret)
library(pROC)
library(dplyr)
library(splitTools)
# setwd("~/Users/jolene/Documents/publications/2024-Ranek_etal_QUICHE")
source("utils.R")
datasets <- list(list(data = read.csv('data/output_files/normalized_cell_abundances_spain.csv'), outcome = "Relapse", positive = 1.0, name = "spain"),
            list(data = read.csv('data/output_files/normalized_cell_abundances_stanford.csv'), outcome = "Recurrence", positive = "POSITIVE", name = "stanford"),
            list(data = read.csv('data/output_files/normalized_cell_abundances_nt.csv'), outcome = "pCR", positive = "RD", name = "nt"),
            list(data = read.csv('data/output_files/normalized_cell_abundances_nt_metacluster.csv'), outcome = "pCR", positive = "RD", name = "nt_metacluster"))

compute_auc <- function(df, outcome_col, positive_class, seed_seq, nfolds = 5) {
  #preprocess
  feature_matrix <- df[, setdiff(names(df), c("Patient_ID", outcome_col))]
  response_label <- df[[outcome_col]]
  response_label_binary <- as.numeric(response_label == positive_class)
  feature_matrix <- as.matrix(feature_matrix)
  feature_matrix[is.na(feature_matrix)] <- 0 
  
  observed_auc_list <- c()
  permuted_auc_list <- c()
  
  #observed
  for (seed_i in seed_seq) {
    set.seed(seed_i)
    folds <- stratified_k_fold_indices(y = response_label_binary, k = nfolds)
    
    fit <- cv.glmnet(feature_matrix, 
                     response_label_binary, 
                     family = 'binomial', 
                     type.measure = 'auc', 
                     foldid = folds, 
                     standardize = FALSE)
    
    observed_auc <- max(fit$cvm)
    observed_auc_list <- c(observed_auc_list, observed_auc)
  }

  #permuted
  for (seed_i in seed_seq) {
    set.seed(seed_i)
    permuted_label <- sample(response_label_binary)
    folds <- stratified_k_fold_indices(y = permuted_label, k = nfolds)
    
    fit <- cv.glmnet(feature_matrix, 
                     permuted_label, 
                     family = 'binomial', 
                     type.measure = 'auc', 
                     foldid = folds, 
                     standardize = FALSE)
    
    permuted_auc <- max(fit$cvm)
    permuted_auc_list <- c(permuted_auc_list, permuted_auc)
  }
  
  auc_data <- data.frame(
    AUC = c(observed_auc_list, permuted_auc_list),
    Label = rep(c(d$outcome, "Permuted"), each = length(observed_auc_list))
  )
  
  return(auc_data)
}

seed_seq <- seq(81, 90)
nfolds <- 5

for (d in datasets) {
  auc_data <- compute_auc(d$data, d$outcome, d$positive, seed_seq, nfolds)
  file_name <- paste0("data/output_files/auc_data_", d$name, ".csv")
  write.csv(auc_data, file = file_name, row.names = FALSE)
  cat(paste("Saved AUC data for", d$name, "to", file_name, "\n"))
}
