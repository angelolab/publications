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
library(precrec)
library(ROCR)
library(svglite)

source("utils.R")

###########################################################
######################### On Nivo #########################
##########################################################

df_X_matrix = read.csv('data/mibi/processed_data/df_matrix_on_nivo.csv', header=TRUE)
colnames(df_X_matrix)[1] = "Patient_ID"

df = read.csv('data/mibi/combined_df.csv', header=TRUE)
df_label = df[,c("Patient_ID", "Clinical_benefit")]
aggregated_df = df_label %>%
  distinct()

merge_df = merge(df_X_matrix, aggregated_df, by="Patient_ID")
merged_df_order = merge_df[, c(names(df_label), setdiff(names(merge_df), names(df_label)))]

feature_matrix = merged_df_order[,3:dim(merged_df_order)[2]]
response_label_raw = merged_df_order[,2]
response_label = response_label_raw
response_label_on_nivo = as.numeric(response_label == "Yes")

# Split into training and test sets with stratification
set.seed(6)
training_perc = 0.7

num_class_1 <- sum(response_label_on_nivo == 1)
num_class_0 <- sum(response_label_on_nivo == 0)

train_num_class_1 <- round(training_perc * num_class_1)
train_num_class_0 <- round(training_perc * num_class_0)

indices_class_1 <- which(response_label_on_nivo == 1)
indices_class_0 <- which(response_label_on_nivo == 0)

indices_class_1 <- sample(indices_class_1)
indices_class_0 <- sample(indices_class_0)

train_indices <- c(indices_class_1[1:train_num_class_1], 
                   indices_class_0[1:train_num_class_0])

test_indices <- setdiff(1:length(response_label_on_nivo), train_indices)

training_features <- feature_matrix[train_indices, ]
test_features <- feature_matrix[test_indices, ]

training_labels <- response_label_on_nivo[train_indices]
test_labels <- response_label_on_nivo[test_indices]

# CV split within training set and fit the model
# standardization in the correct way, only based on the training set
standardized_training_features <- scale(training_features)
train_mean <- attr(standardized_training_features, "scaled:center")
train_sd <- attr(standardized_training_features, "scaled:scale")
standardized_test_features <- scale(test_features, center = train_mean, scale = train_sd)

standardized_training_features = as.matrix(standardized_training_features)
standardized_training_features[is.na(standardized_training_features)] = 0

standardized_test_features = as.matrix(standardized_test_features)
standardized_test_features[is.na(standardized_test_features)] = 0

nfolds = 3
folds =  stratified_k_fold_indices(y = training_labels, k = nfolds)

fit_on_nivo = cv.glmnet(standardized_training_features, 
                        training_labels, 
                        family='binomial', 
                        type.measure='auc', 
                        foldid = folds,
                        standardize=FALSE)

coefs = coef(fit_on_nivo, s="lambda.min")
coefs_no_intercept = coefs[2:length(coefs)]
feature_nonzero = colnames(standardized_training_features)
feature_selected = data.frame('feature'=feature_nonzero, 'coef'=coefs_no_intercept)
feature_selected_order = feature_selected[rev(order(abs(feature_selected[,'coef']))), ]
feature_selected_order_sel = feature_selected_order[abs(feature_selected_order[,'coef']) > 0,]

predictions = predict(fit_on_nivo, newx = standardized_test_features, 
                      s = "lambda.min", type = "response")

# Evaluate the model with AUROC and AUPRC
pred <- prediction(predictions, test_labels)
perf <- performance(pred, "tpr", "fpr")
auroc <- performance(pred, measure = "auc")@y.values[[1]]

pr_perf <- performance(pred, "prec", "rec")
auprc <- sum((pr_perf@x.values[[1]][-1] - pr_perf@x.values[[1]][-length(pr_perf@x.values[[1]])]) * pr_perf@y.values[[1]][-1])

roc_obj = roc(test_labels, predictions)
auprc_obj = evalmod(scores = predictions, labels = test_labels)
auc(auprc_obj)

svglite("results/mibi/on_nivo_auroc.svg")
plot(auprc_obj, "ROC", main="AUROC")
mtext("AUROC = 0.875", side=3, line=0.5, cex=1)
dev.off()

svglite("results/mibi/on_nivo_auprc.svg")
plot(auprc_obj, "PRC", main="AURPC")
mtext("AURPC = 0.822", side=3, line=0.5, cex=1)
dev.off()




