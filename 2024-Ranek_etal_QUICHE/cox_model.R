library(survival)
library(glmnet)
library(caret)
library(survminer)
library(MASS)
options(warn = -1)
setwd("~/Documents/Angelo_lab/quiche")
run_cox_analysis <- function(train_path, test_path, feature_exclude, output_pdf) {
  #preprocess
  df_train <- read.csv(train_path, header = TRUE)
  df_test <- read.csv(test_path, header = TRUE)
  feature_names <- setdiff(names(df_test), feature_exclude)
  X_train <- as.matrix(df_train[, feature_names])
  X_train <- X_train[, colSums(X_train != 0) > 0]
  X_test <- as.matrix(df_test[, feature_names])
  X_test <- X_test[, colSums(X_test != 0) > 0]

  time_train <- as.numeric(df_train$Time_Relapse.days._capped)
  event_train <- as.numeric(df_train$RECURRENCE_LABEL)

  time_test <- as.numeric(df_test$Time_Relapse.days._capped)
  event_test <- ifelse(df_test$RECURRENCE_LABEL == "POSITIVE", 1, 0)
  
  train_data <- data.frame(Time = time_train, Relapse = event_train, X_train)
  test_data <<- data.frame(Time = time_test, Relapse = event_test, X_test)
  
  #fit cox model on Spain cohort
  cox_model <- coxph(Surv(Time, Relapse) ~ ., data = train_data)
  
  #predict risk scores on Stanford cohort
  risk_scores <- predict(cox_model, newdata = test_data, type = "risk")
  
  risk_groups <- cut(risk_scores, 
                     breaks = quantile(risk_scores, probs = c(0, 0.5, 1)), ## med
                     labels = c("Low", "High"))
  
  test_data$risk_groups <- risk_groups
  surv_object <<- Surv(test_data$Time, test_data$Relapse)
  
  #KM+log rank test
  km_fit <<- survfit(surv_object ~ test_data$risk_groups)
  log_rank_test <- survdiff(Surv(time_test, event_test) ~ risk_groups, data = test_data)
  chi_square_stat <- log_rank_test$chisq
  df <- length(log_rank_test$n) - 1
  p_value <- 1 - pchisq(chi_square_stat, df)
  
  pdf(output_pdf, width = 10, height = 8)
  p <- ggsurvplot(km_fit,
    data = test_data,
    risk.table = TRUE,
    pval = paste0("P-value: ", signif(p_value, digits = 2)),
    ggtheme = theme_classic(),
    xlab = "Time (days)",
    tables.height = 0.2,
    ylab = "Relapse Free Survival",
    palette = c("#8FBAE0", "#F97068"),
    xlim = c(0, 5000),
    legend.title = "Risk",
    legend.labs = c("Low", "High"),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.3)
  
  print(p, newpage = FALSE)
  dev.off()
}

run_cox_analysis(train_path = "data/output_files/spain_niche_count.csv",
                test_path = "data/output_files/stanford_niche_count.csv",
                feature_exclude = c("Patient_ID", "Time_Relapse.days._capped", "RECURRENCE_LABEL"),
                output_pdf = "publications/figures/figure5/figure5n_quiche.pdf")

run_cox_analysis(train_path = "data/output_files/spain_celltype_count.csv",
                test_path = "data/output_files/stanford_celltype_count.csv",
                feature_exclude = c("Patient_ID", "Time_Relapse.days._capped", "RECURRENCE_LABEL"),
                output_pdf = "publications/figures/figure5/figure5n_abundance.pdf")

run_cox_analysis(train_path = "data/output_files/spain_kmeans_count.csv",
                test_path = "data/output_files/stanford_kmeans_count.csv",
                feature_exclude = c("Patient_ID", "Time_Relapse.days._capped", "RECURRENCE_LABEL"),
                output_pdf = "publications/figures/figure5/figure5n_kmeans.pdf")