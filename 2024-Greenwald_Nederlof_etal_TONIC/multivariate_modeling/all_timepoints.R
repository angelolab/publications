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

source("utils.R")

######################### ######################### ######################### 
######################### On Nivo #########################
######################### ######################### ######################### 

df_X_matrix = read.csv('data/mibi/processed_data/df_matrix_on_nivo.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('data/mibi/combined_df.csv', header=TRUE)
df_label = df[,c("Patient_ID", "Clinical_benefit")]
aggregated_df = df_label %>%
  distinct()
seed_seq = seq(51,60)

merge_df = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order = merge_df[, c(names(df_label), setdiff(names(merge_df), names(df_label)))]

feature_matrix = merged_df_order[,3:dim(merged_df_order)[2]]
response_label_raw = merged_df_order[,2]
response_label = response_label_raw
response_label_on_nivo = as.numeric(response_label == "Yes")

feature_on_nivo = as.matrix(feature_matrix)
feature_on_nivo[is.na(feature_on_nivo)] = 0

nfolds = 3
auc_on_nivo_list = c()
feature_selected_on_nivo = list()
i = 1
for (seed_i in seed_seq){
  set.seed(seed_i)
  folds =  stratified_k_fold_indices(y = response_label_on_nivo, k = nfolds)
  
  fit_on_nivo = cv.glmnet(feature_on_nivo, 
                          response_label_on_nivo, 
                          family='binomial', 
                          type.measure='auc', 
                          foldid = folds,
                          standardize=FALSE)
  
  auc_on_nivo = max(fit_on_nivo$cvm) 
  auc_on_nivo_list = c(auc_on_nivo_list, auc_on_nivo)
  
  feature_selected_on_nivo[[i]] = get_ranked_feature(fit_on_nivo, feature_on_nivo)[1:20,]
  i = i + 1
}

mean(auc_on_nivo_list)
combined_results <- do.call(cbind, feature_selected_on_nivo)
write.csv(combined_results, file = "./results/mibi/top_features_results_on_nivo_MIBI.csv", row.names = FALSE)
#write.csv(auc_on_nivo_list, file = "./results_100_runs/results_on_nivo_MIBI_100_runs_previous.csv", row.names = FALSE)

######################### ######################### ######################### 
######################### Baseline #########################
######################### ######################### ######################### 

df_X_matrix = read.csv('data/mibi/processed_data/df_matrix_baseline.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

merge_df = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order = merge_df[, c(names(df_label), setdiff(names(merge_df), names(df_label)))]

feature_matrix = merged_df_order[,3:dim(merged_df_order)[2]]
response_label_raw = merged_df_order[,2]
response_label = response_label_raw
response_label_baseline = as.numeric(response_label == "Yes")

feature_baseline = as.matrix(feature_matrix)
feature_baseline[is.na(feature_baseline)] = 0

nfolds = 3
auc_baseline_list = c()
feature_selected_baseline = list()
i = 1

for (seed_i in seed_seq){
  set.seed(seed_i)
  folds =  stratified_k_fold_indices(y = response_label_baseline, k = nfolds)

  fit_baseline = cv.glmnet(feature_baseline, 
                           response_label_baseline, 
                           family='binomial', 
                           type.measure='auc', 
                           foldid = folds,
                           standardize=FALSE)
  
  auc_baseline = max(fit_baseline$cvm) 
  auc_baseline_list = c(auc_baseline_list, auc_baseline)
  
  feature_selected_baseline[[i]] = get_ranked_feature(fit_baseline, feature_baseline)[1:20,]
  i = i + 1
}

combined_results <- do.call(cbind, feature_selected_baseline)
write.csv(combined_results, file = "./results/mibi/top_features_results_baseline_MIBI.csv", row.names = FALSE)


######################### ######################### ######################### 
######################### Post Induction #########################
######################### ######################### ######################### 

df_X_matrix = read.csv('data/mibi/processed_data/df_matrix_post_induction.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

merge_df = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order = merge_df[, c(names(df_label), setdiff(names(merge_df), names(df_label)))]

feature_matrix = merged_df_order[,3:dim(merged_df_order)[2]]
response_label_raw = merged_df_order[,2]
response_label = response_label_raw
response_label_post_induction = as.numeric(response_label == "Yes")

feature_post_induction = as.matrix(feature_matrix)
feature_post_induction[is.na(feature_post_induction)] = 0

auc_post_induction_list = c()
nfolds = 3
feature_selected_post_induction = list()
i = 1

for (seed_i in seed_seq){
  set.seed(seed_i)
  folds = stratified_k_fold_indices(y = response_label_post_induction, k = nfolds)

  fit_post_induction = cv.glmnet(feature_post_induction, 
                                 response_label_post_induction, 
                                 family='binomial', 
                                 type.measure='auc', 
                                 foldid = folds,
                                 standardize=FALSE)
  
  auc_post_induction = max(fit_post_induction$cvm)
  auc_post_induction_list = c(auc_post_induction_list, auc_post_induction)
  feature_selected_post_induction[[i]] = get_ranked_feature(fit_post_induction, feature_post_induction)[1:20,]
  i = i + 1
}

combined_results <- do.call(cbind, feature_selected_post_induction)
write.csv(combined_results, file = "./results/mibi/top_features_results_induction_MIBI.csv", row.names = FALSE)

######################### ######################### ######################### 
######################### Primary #########################
######################### ######################### ######################### 

df_X_matrix = read.csv('data/mibi/processed_data/df_matrix_primary_untreated.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

merge_df = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order = merge_df[, c(names(df_label), setdiff(names(merge_df), names(df_label)))]

feature_matrix = merged_df_order[,3:dim(merged_df_order)[2]]
response_label_raw = merged_df_order[,2]
response_label = response_label_raw
response_label_primary = as.numeric(response_label == "Yes")

feature_primary = as.matrix(feature_matrix)
feature_primary[is.na(feature_primary)] = 0

nfolds = 3
auc_primary_list = c()
feature_selected_primary = list()
i = 1

for (seed_i in seed_seq){
  set.seed(seed_i)
  folds = stratified_k_fold_indices(y = response_label_primary, k = nfolds)
  fit_primary = cv.glmnet(feature_primary,
                                   response_label_primary,
                    family='binomial',
                    type.measure='auc',
                    foldid = folds,
                    standardize=FALSE)
  auc_primary = max(fit_primary$cvm) 
  auc_primary_list = c(auc_primary_list, auc_primary)
  feature_selected_primary[[i]] = get_ranked_feature(fit_primary, feature_primary)[1:20,]
  i = i + 1
}

combined_results <- do.call(cbind, feature_selected_primary)
write.csv(combined_results, file = "./results/mibi/top_features_results_primary_MIBI.csv", row.names = FALSE)


# Save all timepoint resutls in csv 
all_results = cbind(auc_on_nivo_list, auc_baseline_list, auc_post_induction_list, auc_primary_list)
write.csv(all_results, file = "./results/mibi/all_timepoints_results_MIBI.csv", row.names = FALSE)