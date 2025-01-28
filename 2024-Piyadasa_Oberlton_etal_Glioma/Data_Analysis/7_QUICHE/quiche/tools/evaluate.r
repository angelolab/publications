library('limma')
evaluate_graphcompass <- function(condition, interaction_mat, ref){
    cell_type_interaction_mat = data.frame(interaction_mat)
    cell_type_interaction_mat$subject_id <- as.factor(cell_type_interaction_mat$subject_id)
    cell_type_interaction_mat$condition <- as.factor(cell_type_interaction_mat$condition)
    #contrasts(cell_type_interaction_mat$condition) <- contr.treatment(levels(cell_type_interaction_mat$condition), base = which(levels(cell_type_interaction_mat$condition) == ref))
    #cell_type_interaction_mat$condition <- relevel(cell_type_interaction_mat$condition, ref = ref)
    # model <- lm(value ~ subject_id + interaction_id*condition, data=cell_type_interaction_mat)
    #model <- lm(value ~ subject_id + interaction_id*condition, data=cell_type_interaction_mat)
    model <- lm(value ~ interaction_id*condition, data=cell_type_interaction_mat)
    model_df <- as.data.frame(coef(summary(model)))

    return(model_df)
}