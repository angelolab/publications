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


######################### #############################
######################### Baseline ####################
######################### RNA #########################
######################### #############################

df_X_matrix = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/multivariate_lasso/rna_df_matrix_baseline_signature.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/processed_genomics_features.csv', header=TRUE)
df_label = df[,c("Patient_ID", "Clinical_benefit")]
aggregated_df = df_label %>%
  distinct()
seed_seq = seq(51,60)

merge_df_rna = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order_rna = merge_df_rna[, c(names(df_label), setdiff(names(merge_df_rna), names(df_label)))]

feature_matrix_rna = merged_df_order_rna[,3:dim(merged_df_order_rna)[2]]
# Need to scale because not on the same unit
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

feature_selected_baseline_rna = list()
i = 1

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

  feature_selected_baseline_rna[[i]] = get_ranked_feature(fit_baseline, feature_baseline_rna)[1:20,]
  i = i + 1
}

combined_results <- do.call(cbind, feature_selected_baseline_rna)
write.csv(combined_results, file = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/top_features_results_baseline_RNA_signature.csv", row.names = FALSE)

######################### #############################
######################### Induction ####################
######################### RNA #########################
######################### #############################

df_X_matrix = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/multivariate_lasso/rna_df_matrix_post_induction_signature.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/processed_genomics_features.csv', header=TRUE)
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

response_label_induction_rna = as.numeric(response_label == "Yes")
feature_induction_rna = as.matrix(feature_matrix_rna)
#feature_baseline_rna = scale(feature_baseline_rna)
feature_induction_rna[is.na(feature_induction_rna)] = 0

nfolds = 3
auc_baseline_list_rna_induction = c()

feature_selected_induction_rna = list()
i = 1

for (seed_i in seed_seq){
  set.seed(seed_i)
  folds =  stratified_k_fold_indices(y = response_label_induction_rna, k = nfolds)

  fit_baseline = cv.glmnet(feature_induction_rna, #feature_baseline,
                           response_label_induction_rna,
                           family='binomial',
                           type.measure='auc',
                           foldid = folds,
                           standardize=FALSE)

  auc_baseline = max(fit_baseline$cvm)
  auc_baseline_list_rna_induction = c(auc_baseline_list_rna_induction, auc_baseline)

  feature_selected_induction_rna[[i]] = get_ranked_feature(fit_baseline, feature_induction_rna)[1:20,]
  i = i + 1
}

combined_results <- do.call(cbind, feature_selected_induction_rna)
write.csv(combined_results, file = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/top_features_results_inducion_RNA_signature.csv", row.names = FALSE)

######################### #############################
######################### On-nivo ####################
######################### RNA #########################
######################### #############################

df_X_matrix = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/multivariate_lasso/rna_df_matrix_on_nivo_signature.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/processed_genomics_features.csv', header=TRUE)
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

response_label_nivo_rna = as.numeric(response_label == "Yes")
feature_nivo_rna = as.matrix(feature_matrix_rna)
#feature_baseline_rna = scale(feature_baseline_rna)
feature_nivo_rna[is.na(feature_nivo_rna)] = 0

nfolds = 3
auc_baseline_list_rna_nivo = c()

feature_selected_nivo_rna = list()
i = 1

for (seed_i in seed_seq){
  set.seed(seed_i)
  folds =  stratified_k_fold_indices(y = response_label_nivo_rna, k = nfolds)

  fit_baseline = cv.glmnet(feature_nivo_rna, #feature_baseline,
                           response_label_nivo_rna,
                           family='binomial',
                           type.measure='auc',
                           foldid = folds,
                           standardize=FALSE)

  auc_baseline = max(fit_baseline$cvm)
  auc_baseline_list_rna_nivo = c(auc_baseline_list_rna_nivo, auc_baseline)

  feature_selected_nivo_rna[[i]] = get_ranked_feature(fit_baseline, feature_nivo_rna)[1:20,]
  i = i + 1
}

combined_results <- do.call(cbind, feature_selected_nivo_rna)
write.csv(combined_results, file = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/top_features_results_nivo_RNA_signature.csv", row.names = FALSE)

######################### #############################
######################### Baseline ####################
######################### DNA #########################
######################### #############################

df_X_matrix = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/multivariate_lasso/dna_df_matrix_baseline.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/processed_genomics_features.csv', header=TRUE)
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
auc_baseline_list_dna = c()

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
  auc_baseline_list_dna = c(auc_baseline_list_dna, auc_baseline)
}

# save all results in csv
all_results = cbind(auc_baseline_list_rna, auc_baseline_list_rna_induction, auc_baseline_list_rna_nivo, auc_baseline_list_dna)
colnames(all_results) = c("RNA Baseline", "RNA Induction", "RNA Nivo", "DNA Baseline")
write.csv(all_results, file = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data/prediction/multivariate_lasso/all_timepoints_results_genomics.csv", row.names = FALSE)




