library(dplyr)
library("survival")
library("survminer")
library(RColorBrewer)
library(ggplot2)
library(splines)
library(car)
library(nnet)

############################################################### read data  ###################################################################################
rahbt <- read.csv("/Users/innaa/Library/CloudStorage/GoogleDrive-innaa@stanford.edu/My Drive/DCIS 2.0 Masking/Tables/hazard_ratio_input_table.csv", header = T, sep = ",",stringsAsFactors = F ,fill=T)
rahbt$Treatment_type[rahbt$Treatment_type == 'Lumpectomy_RT_unknown'] <- NA


# Create a binary event indicator
rahbt$recurrence = ifelse(rahbt$event_recur_type == "Invasive_Ipsilateral", 1, 
                          ifelse(rahbt$event_recur_type == "Non-progressor", 0, NA))

############################################################################################################################################################
############################################################### Testing recurrence  ########################################################################
############################################################################################################################################################
t4 = coxph(Surv(Mo_FU, recurrence) ~ 
             Race + Age_at_diagnosis + Overall_tumor_grade + 
             Treatment_type + ER_RNA + Her2_RNA + ER_clin + PR_clin, 
           data = rahbt)

ggforest(t4)

# Check model summary
summary(t4)

# Check if model converged and coefficients
print(t4)

########################## full model does not converge because of missing ER/PR/Her2 values - split into simpler models: #####################################
# model with fewer variables - no ER/PR/Her2
t4_simple = coxph(Surv(Mo_FU, recurrence) ~ 
                    Age_at_diagnosis + Race + Overall_tumor_grade + Treatment_type, 
                  data = rahbt)
summary(t4_simple)
ggforest(t4_simple)

# model with fewer variables - RNA based tumor annotation
t4_simple = coxph(Surv(Mo_FU, recurrence) ~ 
                    ER_RNA + Her2_RNA , 
                  data = rahbt)
summary(t4_simple)
ggforest(t4_simple)

# models with fewer variables - clinical tumor annotation
t4_simple = coxph(Surv(Mo_FU, recurrence) ~ 
                    ER_clin + PR_clin , 
                  data = rahbt)
summary(t4_simple)
ggforest(t4_simple)

#################################################### Logistic regression for association of recurrence with FOV properties ################################
model_simple <- glm(recurrence ~ 
                      N_cells_fov + N_cells_emask + size_emask + 
                      N_cells_smask + size_smask + size_fov,
                    family = binomial(link = "logit"),
                    data = rahbt)

# Summary of results
summary(model_simple)
library(broom)

# Create forest plot for logistic regression
tidy_results <- tidy(model_simple, conf.int = TRUE, exponentiate = TRUE)
ggplot(tidy_results, aes(y = term, x = estimate)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Odds Ratio", y = "Variable")

############################################################################################################################################################
################################### Test meta data association with  score  ################################################################################
############################################################################################################################################################
# using clinical ER
lm_result = lm(score ~ Age_at_diagnosis + Race + 
                 ER_clin + PR_clin + Her2_RNA+ Treatment_type + Overall_tumor_grade , data=rahbt)
summary(lm_result)
anova(lm_result)

######### Check model assumptions ########
par(mfrow=c(2,2))
plot(lm_result) 

# check for collinearity
vif(lm_result)

######### using RNA based ER ############
lm_result = lm(score ~ Age_at_diagnosis + Race + 
                 PR_clin + Her2_RNA + ER_RNA+ Treatment_type + Overall_tumor_grade , data=rahbt)
summary(lm_result)
anova(lm_result)

# Check model assumptions
par(mfrow=c(2,2))
plot(lm_result) 

# check for collinearity
vif(lm_result)


# Look at metadata variables that associate with score to quantify how much of the variance in score they explain 
lm_result = lm(score ~ Age_at_diagnosis  , data=rahbt)
summary(lm_result)
anova(lm_result)


lm_result = lm(score ~ Overall_tumor_grade  , data=rahbt)
summary(lm_result)
anova(lm_result)

lm_result = lm(score ~ Overall_tumor_grade +Age_at_diagnosis  , data=rahbt)
summary(lm_result)
anova(lm_result)

####################################### Testing association between score and FOV properties ##########################################
lm_result = lm(score ~  N_cells_fov + N_cells_emask + size_emask + 
                 N_cells_smask + size_smask + size_fov , data=rahbt)
summary(lm_result)
anova(lm_result)

# Evaluate how much of the variance in score is explained by each FOV property
evaluate_individual_predictors <- function(data, dependent_var, predictors, plot_title = "Individual Predictor Contributions") {
  # Initialize a list to store results
  results_list <- list()
  
  # Calculate R-squared for each predictor individually
  for(predictor in predictors) {
    formula <- as.formula(paste(dependent_var, "~", predictor))
    model <- lm(formula, data = data)
    summary_stats <- summary(model)
    
    # Append results to the list
    results_list[[predictor]] <- data.frame(
      Predictor = predictor,
      R_squared = summary_stats$r.squared,
      P_value = anova(model)$`Pr(>F)`[1]
    )
  }
  
  # Combine list into a single data frame
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL # Clean up row names
  
  # Add percentage column
  results$R_squared_percent <- results$R_squared * 100
  
  # Create a factor for coloring
  results$Significance <- ifelse(results$P_value < 0.05, "p < 0.05", "p >= 0.05")
  
  # Sort data for the table output
  results_sorted_for_table <- results[order(-results$R_squared), ]
  
  # Create bar plot
  plot <- ggplot(results, 
                 aes(x = reorder(Predictor, R_squared_percent), 
                     y = R_squared_percent,
                     fill = Significance)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(name = "P-value", # Legend title
                      values = c("p < 0.05" = "darkblue", "p >= 0.05" = "lightblue")) +
    coord_flip() +
    labs(x = "Predictor", 
         y = "Variance Explained (%)",
         title = plot_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = "bottom" 
    )
  
  # Print numerical results
  cat("\nIndividual Predictor Analysis Results:\n")
  cat("----------------------------------------\n")
  print(results_sorted_for_table[, c("Predictor", "R_squared_percent", "P_value")])
  
  # Return both results and plot
  return(list(
    results = results_sorted_for_table,
    plot = plot
  ))
}

predictors <- c("N_cells_fov", "N_cells_emask", "size_emask",
                "N_cells_smask", "size_smask", "size_fov")

results <- evaluate_individual_predictors(
  data = rahbt,
  dependent_var = "score",
  predictors = predictors,
  plot_title = "Variance Explained by Individual Predictors"
)

print(results$plot)

####################################### Testing association between tumor grade and treatment type ##########################################

# First create ordinal grade
rahbt$grade_ordinal <- as.numeric(factor(rahbt$Overall_tumor_grade, 
                                       levels=c("Well differentiated - Grade I - 3-5 points",
                                              "Moderately differentiated - Grade II - 6-7 points",
                                              "Poorly differentiated - Grade III - 8-9 points")))

# Run GLM
glm1 <- multinom(Treatment_type ~ grade_ordinal, data=rahbt)
summary(glm1)
# Calculate z-scores
z <- summary(glm1)$coefficients/summary(glm1)$standard.errors

# Calculate p-values
p <- (1 - pnorm(abs(z), 0, 1)) * 2

# Print p-values
print(p)

# Overall model significance - likelihood ratio test
Anova(glm1, type="II")


# or
glm1 <- glm(grade_ordinal ~ Treatment_type, family=gaussian, data=rahbt)
summary(glm1)
anova(glm1)


#### try grouping grade and treatment type ####
# Create binary grade
rahbt$grade_binary <- ifelse(rahbt$Overall_tumor_grade == "Well differentiated - Grade I - 3-5 points", 
                            "Low", "High")

# Create binary treatment
rahbt$treatment_binary <- ifelse(rahbt$Treatment_type == "Lumpectomy_no_RT", 
                                "No_RT", "RT_or_Mastectomy")

# Create contingency table
cont_table <- table(rahbt$grade_binary, rahbt$treatment_binary)

# Chi-square test
chisq.test(cont_table)

# Fisher's exact test 
fisher.test(cont_table)
