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
seed_seq = seq(1,100)

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
hist(auc_on_nivo_list, main="Histogram of AUC on Nivo (Previous)", xlab="AUC", 
     col="lightblue", border="black")
text(mean(auc_on_nivo_list), 10, paste("Mean AUC: ", round(mean(auc_on_nivo_list), 2)), 
     adj=c(1.2, -10.9))

combined_results <- do.call(cbind, feature_selected_on_nivo)
write.csv(combined_results, file = "./top_features_results_on_nivo_MIBI_100_runs_previous.csv", row.names = FALSE)
write.csv(auc_on_nivo_list, file = "./results_on_nivo_MIBI_100_runs_previous.csv", row.names = FALSE)
