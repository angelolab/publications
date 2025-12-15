# Clear workspace
rm(list = ls(all = TRUE))

# Restart R
rstudioapi::restartSession()

#################### Load necessary packages ########################
require(ranger)
require(caret)
require(survival)
require(pheatmap)
require(pROC)
require(ggpubr)
require(ggbeeswarm)
require(survminer)
require(tidyverse)
require(data.table)

library(glmnet)
library(caret)
library(ggplot2)
library(reshape2)
library(smotefamily)

#################### functions ########################

#### Impute nan values #########
imputeValues <- function(tr, te, ti) {
  for(marker in ti) {
    tr[is.na(eval(parse(text=marker))), eval(marker):=mean(tr[, eval(parse(text=marker))], na.rm=T)]
    te[is.na(eval(parse(text=marker))), eval(marker):=mean(te[, eval(parse(text=marker))], na.rm=T)]
  }
}

#### Classifier #########
lasso_classifier_model <- function(train_set, test_set,alpha_fac, discard_tr) {
  
  # Separate feature matrix and response vector in training set
  x <- model.matrix(Status~., train_set)[,-1]
  
  center <- colMeans(x)
  scale <- apply(x, 2, sd)
  
  # Scale training data
  x_normalized <- t((t(x) - center) / scale)
  y <- train_set$Status
  
  # Scale test data
  x.test <- model.matrix(Status~., test_set)[,-1]
  x_test_normalized <- t((t(x.test) - center) / scale)
  
  #################### remove non predictive variables using training set ########################
  predictors_to_discard <- c()

  # For each predictor, perform a Kruskal-Wallis test
  for(predictor in colnames(x_normalized)) {
    tryCatch({
      # "Try" part
      # message(paste("Processing predictor:", predictor))
      data_for_test <- data.frame(Predictor = x_normalized[, predictor],
                                  Status = train_set$Status)
      # Check if all predictor values in test_set are NA
      if(any(is.na(test_set[[predictor]]))) {
        predictors_to_discard <- c(predictors_to_discard, predictor)
        message(paste("discarding based on test:", predictor))
      } else {
        kruskal_test_result <- kruskal.test(Predictor ~ Status, data = data_for_test)
        if(kruskal_test_result$p.value > discard_tr) {
          predictors_to_discard <- c(predictors_to_discard, predictor)
          # message(paste("predictor ", predictor ))
          # message(paste("p-val ", kruskal_test_result$p.value))
        }
      }
    },
    error = function(cond) {
      # "Error" part
      message(paste("Failed to process predictor:", predictor))
      message("Here's the original error message:")
      message(conditionMessage(cond))
      predictors_to_discard <- c(predictors_to_discard, predictor)
    },
    finally = {
      # "Finally" part
      # message(paste("Processed predictor:", predictor))
    }
    )
  }

  # print(predictors_to_discard)
  
  # Exclude predictors based on Krukal-Wallis test result
  x_normalized_filtered <- x_normalized[, !(colnames(x_normalized) %in% predictors_to_discard)]
  x_test_normalized_filtered <- x_test_normalized[, colnames(x_test_normalized) %in% colnames(x_normalized_filtered)]
  
  #################### train model ########################
  
  # Train the Lasso model
  cv.lasso <- cv.glmnet(x_normalized_filtered, y, family="binomial", alpha=alpha_fac, nfolds=10)
  
  # Use the lambda that gives the minimum mean cross-validated error
  best.lambda <- cv.lasso$lambda.min
  
  # Fit final model on the training data
  weights <- ifelse(y == 1, 5, 1)
  lasso.model <- glmnet(x_normalized_filtered, y, family="binomial", alpha=alpha_fac, lambda=best.lambda, weights=weights)
  
  #################### predict test set ########################
  predictions <- predict(lasso.model, newx =  x_test_normalized_filtered, type = "response")
  # Create the response vector for the test set
  y.test <- test_set$Status
  
  #################### evaluate and plot model performance ########################
  
  # Calculate metrics for different thresholds
  thresholds <- seq(0, 1, by = 0.01)
  metrics_results <- sapply(thresholds, function(thresh) {
    predicted.classes <- ifelse(predictions > thresh, 1, 0)

    confusionMatrix <- table(
      Predicted = factor(predicted.classes, levels = c(0,1)),
      Actual = factor(y.test, levels = c(0,1))
    )
    #print(confusionMatrix)  # see the full confusion matrix
    
    # Calculate NPV
    if(sum(confusionMatrix[1,]) > 0) {
      NPV <- confusionMatrix[1,1] / sum(confusionMatrix[1,])
    } else {
      NPV <- 0
    }
    
    # Calculate Specificity (TNR)
    if((confusionMatrix[1,1] + confusionMatrix[2,1]) > 0) {
      Specificity <- confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[2,1])
    } else {
      Specificity <- 0
    }
    
    return(c(NPV = NPV, Specificity = Specificity))
  })
  
  # Convert results to data frame
  metrics_df <- data.frame(
    threshold = thresholds,
    NPV = metrics_results[1,],
    Specificity = metrics_results[2,]
  )
  
  # Find best NPV and its specificity
  best_npv <- max(metrics_df$NPV)
  
  # Get all rows that have the best NPV
  best_npv_rows <- metrics_df[metrics_df$NPV == best_npv, ]
  
  # Among those rows, find the best Specificity
  best_npv_Specificity <- max(best_npv_rows$Specificity)
  #best_npv_Specificity <- metrics_df$Specificity[which.max(metrics_df$NPV)]
  
  
  # Calculate metrics
  predictions <- as.vector(predictions)
  roc_curve <- roc(y.test, predictions)
  auc <- auc(roc_curve)
  
  cat("AUC-ROC: ", auc, "\n")
  cat("Best NVP: ", best_npv, "\n")
  cat("Specificity for best NVP: ",  best_npv_Specificity, "\n")
  
  
  #################### analyze model coefficients ########################
  
  # Extract coefficients
  coef_lasso <- coef(lasso.model, s = best.lambda)
  coef_lasso_vector <- as.vector(coef_lasso)
  sorted_coefs <- sort(coef_lasso_vector, decreasing = TRUE, index.return = TRUE)
  
  # Get names of the predictors corresponding to the sorted coefficients
  predictor_names <- rownames(coef_lasso)[sorted_coefs$ix]
  var_names <- rownames(coef_lasso)
  coefficients <- as.numeric(coef_lasso)
  
  # Combine in a dataframe
  coef_df <- data.frame(Predictor = var_names, Coefficient = coefficients)
  coef_df <- coef_df[order(abs(coef_df$Coefficient)), ]
  
  # Add a new column with absolute values of coefficients
  coef_df$Abs_Coefficient <- abs(coef_df$Coefficient)
  
  # Exclude the intercept
  coef_df <- coef_df[coef_df$Predictor != "(Intercept)", ]
  
  return(list(auc = auc, coef_df = coef_df, roc_obj = roc_curve, metrics_by_threshold = metrics_df ))
}



################################################################################################
################################################################################################
#################### read input ########################
################################################################################################
# directory for saving results
path <- "/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/Lasso classifier tables/"    

#run_name <- 'with_compartment_score_st'
# tables subdirectory
#tables.path <- paste0(path,run_name, "/")
#table.file <- paste0(tables.path, "feature_table_" ,run_name, ".csv")

# read feature table
table.file <-'/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/__Tables to upload/Lasso classifier tables/classifier_input_table_features.csv'
drop.cols <- c("diag_short","event_timing","fov", "event_recur_type", "Mastectomy","Radiation","Tamoxifen","Years_Since_First_CIS","Days_To_First_Event","PointNumber")
data <- fread(table.file) %>%
  .[order(PointNumber)] %>%
  .[, setdiff(colnames(.), drop.cols), with=F]

# Replace placeholder with NA
data[data == -999999] <- NA

#################### arrange data ########################

# Convert Status to a factor 
data$Status <- factor(data$Status, levels = c("case", "ctrl"))
# Check the levels 
levels(data$Status)
id.cols <- c( "Status")
unfiltered.predictors <- setdiff(colnames(data), c(id.cols))
# Convert status to binary
data$Status <- ifelse(data$Status=="case", 1, 0)

################################################################################################
#########################              LOOP OVER SEEDS             #############################
################################################################################################
# Define a vector with all different seeds
seeds <- c(123, 124,125, 457, 852, 147, 753, 258,130,963)
use_alpha <- 1
use_discrard_p_tr <- 0.99

# Determine model type for file names 
if (use_alpha == 1) {
  model_type <- 'Lasso'
} else if (use_alpha > 0 && use_alpha < 1) {
  model_type <- 'Enet'
} else if (use_alpha == 0) {
  model_type <- 'Ridge'
}

# Initialize a table to store the auc_cur values
auc_table <- data.frame(seed = numeric(), 
                        auc = numeric(), 
                        auc_permuted = numeric())
# Initialize the list before the loop
coef_dfs <- list()
coef_dfs_permuted <- list()
roc_objs <- list()
metrics_by_threshold <- list()
###########################################################################
#################### Loop over each seed starts here ######################
###########################################################################
for(i in seq_along(seeds)) {
  # Set the random seed
  set.seed(seeds[i])
  # Split the data into training and testing sets - ensuring enough cases in training set (~80-20)
  N_ctrl_test=28 
  N_case_test=4 
  testing.id <- c(sample(which(data$Status==0), N_ctrl_test),
                  sample(which(data$Status==1), N_case_test))
  training.id <- setdiff(1:nrow(data), testing.id)
  train_cur <- copy(data[training.id])
  test_cur <-  copy(data[testing.id])
  # Imputing values
  to.impute <-data[, lapply(.SD, function(x) any(is.na(x))), .SDcols=unfiltered.predictors] %>%
    t() %>%
    as.data.table(keep.rownames="predictor") %>%
    .[V1==T, predictor]
  imputeValues(tr=train_cur, te=test_cur, ti=to.impute)

  # create synthetic samples
  set.seed(123) # For reproducibility
  smote_result <- SMOTE(X = train_cur, 
                        target = train_cur$Status, 
                        K = 5, 
                        dup_size = 5)
  # Combine the new synthetic samples with the original data
  training_set_balanced <- cbind(smote_result$data, Status = smote_result$target)
  # Remove last column
  training_set_balanced <- training_set_balanced[,1:(length(training_set_balanced)-1)]
  
  
  # Run classifier
  result <- lasso_classifier_model(train_set=training_set_balanced,test_set = test_cur,alpha_fac=use_alpha,discard_tr = use_discrard_p_tr)
  # Access the auc
  auc_cur <- result$auc
  # Add the updated coef_df to the list
  coef_dfs[[i]] <- result$coef_df
  roc_objs[[i]] <- result$roc_obj
  metrics_by_threshold[[i]] <- result$metrics_by_threshold
  
  # ### Permuted run ##########
  train_permuted <- training_set_balanced
  # Permute the 'Status' labels
  train_permuted$Status <- sample(training_set_balanced$Status)
  
  # run classifier on permuted data
  result_permuted <- lasso_classifier_model(train_set=train_permuted,test_set = test_cur,alpha_fac=use_alpha,discard_tr = use_discrard_p_tr)
  
  # Access the auc
  auc_permuted <- result_permuted$auc
  # Add the updated coef_df to the list
  coef_dfs_permuted[[i]] <- result_permuted$coef_df
  
  # Save AUCs
  auc_table <- rbind(auc_table, data.frame(seed = seeds[i], auc = auc_cur, auc_permuted = auc_permuted))
  
}

# After the loop, merge all data frames in the list by the common column, which is 'variable'
merged_coef_dfs <- Reduce(function(x, y) merge(x, y, by = "Predictor", all = TRUE), coef_dfs)
merged_coef_dfs[is.na(merged_coef_dfs)] <- 0
merged_coef_dfs_permuted <- Reduce(function(x, y) merge(x, y, by = "Predictor", all = TRUE), coef_dfs_permuted)
merged_coef_dfs_permuted[is.na(merged_coef_dfs_permuted)] <- 0

# Change the data form
long_form_data <- melt(auc_table, id.vars = "seed")

# Perform paired Wilcoxon test
paired_wilcox_test_result <- wilcox.test(auc_table$auc, auc_table$auc_permuted, paired = TRUE)
p.value <- paired_wilcox_test_result$p.value
median_auc <- median(auc_table$auc)
min.y <- floor(min(long_form_data$value) * 10)/10
# Plot the data
ggplot(long_form_data, aes(variable, value, fill = variable)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.05, seed=123)) +
  labs(x="Permuted", y="AUC") +
  scale_y_continuous(limits=c(min.y, 1)) +
  annotate("text", x = 1.5, y = 0.97, label = paste0("p = ", round(p.value, 2)), hjust = 0.5, vjust = -0.2, size = 6) +
  annotate("segment", x = 1, xend = 2, y = 0.97, yend = 0.97, colour = "black") +
  ggtitle(paste("Median AUC: ", round(median_auc, digits=3))) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="none")



#################################################################################################

#################################################################################################


#################################################################################################
#######################       Save results      #################################################
ggsave(paste0(tables.path,model_type,"_permutation_analysis_", run_name,".png"), height=5, width=2.5)
fwrite(merged_coef_dfs, paste0(tables.path,model_type, "_merged_coef_dfs_", run_name,".csv"))
fwrite(auc_table, paste0(tables.path,model_type, "_auc_table_", run_name,".csv"))
fwrite(merged_coef_dfs_permuted, paste0(tables.path, model_type, "_merged_coef_dfs_permuted_", run_name,".csv"))

