library(glmnet)
library(caret)
library(xgboost)
library(pROC)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(dplyr)
library(MASS)
library(splitTools)
library(multiview)

source("utils_rna.R")

######################### ############################# 
######################### Baseline ####################
######################### DNA #########################
######################### #############################

df_X_matrix = read.csv('data/genomics/dna_df_matrix_baseline.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('data/genomics/processed_genomics_features_final.csv', header=TRUE)
df_label = df[,c("Patient_ID", "Clinical_benefit")]
aggregated_df = df_label %>%
  distinct()
seed_seq = seq(51,60)

merge_df_rna = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order_rna = merge_df_rna[, c(names(df_label), setdiff(names(merge_df_rna), names(df_label)))]

feature_matrix_rna = merged_df_order_rna[,3:dim(merged_df_order_rna)[2]]
feature_matrix_rna = scale(feature_matrix_rna)
rownames(feature_matrix_rna) = merged_df_order_rna$Patient_ID
response_label_raw = merged_df_order_rna[,2]
response_label = response_label_raw

response_label_baseline_rna = as.numeric(response_label == "Yes")
feature_baseline_rna = as.matrix(feature_matrix_rna)
#feature_baseline_rna = scale(feature_baseline_rna)
feature_baseline_rna[is.na(feature_baseline_rna)] = 0

nfolds = 3
auc_baseline_list_rna = c()

for (seed_i in seed_seq){
  set.seed(seed_i)
  folds =  stratified_k_fold_indices(y = response_label_baseline_rna, k = nfolds)
  
  fit_baseline = cv.glmnet(feature_baseline_rna, #feature_baseline, 
                           response_label_baseline_rna, 
                          family='binomial', 
                          type.measure='auc', 
                          foldid = folds,
                          standardize=FALSE)
  
  auc_baseline = max(fit_baseline$cvm) 
  auc_baseline_list_rna = c(auc_baseline_list_rna, auc_baseline)
}


