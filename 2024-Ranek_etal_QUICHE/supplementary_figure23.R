library(survival)
library(survminer)
library(MASS)
library(ggplot2)
options(warn = -1)

preprocess_features <- function(data, outcome_vars = c("Time", "Relapse"), min_nonzero = 1) {
  feature_data <- data[, !(names(data) %in% outcome_vars)]
  valid_cols <- sapply(feature_data, function(col) {
    nonzero_count <- sum(col != 0)
    unique_vals <- length(unique(col))
    return(nonzero_count >= min_nonzero && unique_vals > 1)
  })
  clean_data <- data[, c(outcome_vars, names(feature_data)[valid_cols])]
  return(clean_data)
}

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
  train_data <- preprocess_features(train_data)
  
  test_data <- data.frame(Time = time_test, Relapse = event_test, X_test)
  test_data <<- test_data[, colnames(test_data) %in% colnames(train_data)]  # Match train features
  
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
  log_rank_test <- survdiff(Surv(time_test, event_test) ~ risk_groups, data = test_data, rho = 1)
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
                  palette = c("Low" = "#8FBAE0", "High" = "#F97068"),
                  xlim = c(0, 5000),
                  legend.labs = c("Low", "High"),
                  legend.title = "Risk",
                  risk.table.y.text.col = TRUE,
                  risk.table.height = 0.3,
                  risk.table.y.text.fontsize = 12)
  
  print(p, newpage = FALSE)

  #permutation test
  set.seed(42)
  n_perm <- 10000
  perm_stats <- numeric(n_perm)
  for (i in 1:n_perm) {
    perm_data <- test_data
    perm_data$Relapse <- sample(perm_data$Relapse)
    perm_test <- survdiff(Surv(perm_data$Time, perm_data$Relapse) ~ perm_data$risk_group, rho = 1)
    perm_stats[i] <- perm_test$chisq
  }
  empirical_p <- mean(perm_stats >= chi_square_stat)
  
  #permutation histogram
  hist(perm_stats, breaks = 30,
       main = "Permutation Test (Chi-squared Distribution)",
       xlab = "Chi-squared under null", col = "#AAAAAA", border = "white")
  abline(v = chi_square_stat, col = "red", lwd = 2)
  legend("topright",
         legend = c(
           bquote("Observed " ~ chi^2 ~ "=" ~ .(round(chi_square_stat, 2))),
           bquote("Empirical p = " ~ .(signif(empirical_p, 2)))
         ),
         bty = "n",
         text.col = "black")
  
  #bootstrap
  n_boot <- 10000
  boot_pvals <- numeric(n_boot)
  for (i in 1:n_boot) {
    boot_idx <- sample(1:nrow(test_data), replace = TRUE)
    boot_data <- test_data[boot_idx, ]
    boot_surv <- Surv(boot_data$Time, boot_data$Relapse)
    boot_test <- survdiff(boot_surv ~ boot_data$risk_group, rho = 1)
    boot_chi <- boot_test$chisq
    boot_df <- length(boot_test$n) - 1
    boot_pvals[i] <- 1 - pchisq(boot_chi, boot_df)
  }
  
  #bootstrap histogram
  hist(boot_pvals, breaks = 30,
       main = "Bootstrap Distribution of Log-rank P-values",
       xlab = "Log-rank p-value", col = "#AAAAAA", border = "white")
  abline(v = p_value, col = "red", lwd = 2)
  legend("topright",
         legend = c(
           bquote("Observed p = " ~ .(signif(p_value, 2))),
           bquote("Bootstrap median = " ~ .(signif(median(boot_pvals), 2)))
         ),
         bty = "n", text.col = "black")
  dev.off()
}

run_cox_analysis(train_path = "data/output_files/spain_niche_count.csv",
                 test_path = "data/output_files/stanford_niche_count.csv",
                 feature_exclude = c("Patient_ID", "Time_Relapse.days._capped", "RECURRENCE_LABEL"),
                 output_pdf = "publications/supplementary_figures/supplementary_figure23_quiche.pdf")

run_cox_analysis(train_path = "data/output_files/spain_celltype_count.csv",
                 test_path = "data/output_files/stanford_celltype_count.csv",
                 feature_exclude = c("X", "Patient_ID", "Time_Relapse.days._capped", "RECURRENCE_LABEL"),
                 output_pdf = "publications/supplementary_figures/supplementary_figure23_abundance.pdf")

run_cox_analysis(train_path = "data/output_files/spain_kmeans_count.csv",
                 test_path = "data/output_files/stanford_kmeans_count.csv",
                 feature_exclude = c("X", "Patient_ID", "Time_Relapse.days._capped", "RECURRENCE_LABEL"),
                 output_pdf = "publications/supplementary_figures/supplementary_figure23_kmeans.pdf)