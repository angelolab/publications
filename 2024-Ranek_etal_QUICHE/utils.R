library(glmnet)

stratified_k_fold_indices <- function(y, k) {
  # Initialize a vector to hold fold assignments for each data point
  fold_assignments <- numeric(length(y))
  
  # Get the unique class labels
  unique_classes <- unique(y)
  
  # Loop through each class for stratified splitting
  for (cl in unique_classes) {
    indices <- which(y == cl)
    len <- length(indices)
    size_per_fold <- floor(len / k)
    remaining <- len - (k * size_per_fold)
    
    # Shuffle indices for randomness
    shuffled_indices <- sample(indices)
    
    start_idx <- 1
    for (i in 1:k) {
      end_idx <- start_idx + size_per_fold - 1
      if (remaining > 0) {
        end_idx <- end_idx + 1
        remaining <- remaining - 1
      }
      
      # Assign the fold number to the corresponding indices
      fold_assignments[shuffled_indices[start_idx:end_idx]] <- i
      start_idx <- end_idx + 1
    }
  }
  
  return(fold_assignments)
}

get_ranked_feature <- function(fit, df, s='lambda.min'){
  coefs = coef(fit, s=s)
  coefs_no_intercept = coefs[2:length(coefs)]
  feature_nonzero = colnames(df) #[coef_fit_nonzero_ind]
  feature_selected = data.frame('feature'=feature_nonzero, 'coef'=coefs_no_intercept)
  feature_selected_order = feature_selected[rev(order(abs(feature_selected[,'coef']))), ]
  feature_selected_order_sel = feature_selected_order[abs(feature_selected_order[,'coef']) > 0,]
  return(feature_selected_order_sel)
}