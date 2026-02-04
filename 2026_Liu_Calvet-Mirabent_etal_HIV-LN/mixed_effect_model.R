# Mixed effect model to compare features
# Author: Candace Liu

library(data.table)
library(nlme)

feature_tab_path = "../data/tables/feature_tab.csv"
metadata_path = "../data/tables/metadata.csv"

feature_tab = fread(feature_tab_path)
metadata = fread(metadata_path)
feature_tab = metadata[,c("fov","sample_id")][feature_tab, on=c("fov","sample_id")]

feature_tab[,fov:=NULL]
all_features = colnames(feature_tab)
all_features = setdiff(all_features, c("status","sample_id","status_with_viremia","status_with_viremia","status_with_p24"))

# Mixed effect model
fit_model <- function(feature) {
  print(feature)
  formula = as.formula(paste(feature,"~status"))
  model = lme(formula, random= ~1|sample_id, data = feature_tab, na.action=na.omit)
  p_value = summary(model)$tTable['statushiv_pos','p-value']
  fixeff_est = summary(model)$tTable['statushiv_pos','Value']
  return(c(estimate=fixeff_est, p_value = p_value))
}

results_matrix = sapply(all_features, fit_model)
results_dt = data.table(t(results_matrix))
results_dt[,feature_name:=colnames(results_matrix)]
results_dt[, p_adj_bh := p.adjust(p_value, method = "BH")]

fwrite(results_dt, "../data/tables/mixed_effect_model_output.csv")
