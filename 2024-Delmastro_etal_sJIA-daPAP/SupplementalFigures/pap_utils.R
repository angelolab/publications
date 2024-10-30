
# pagebreak
# description: 
#
# input: 
# 
# output: 
#
pagebreak <- function() {
  if(knitr::is_latex_output())
    return("\\newpage")
  else
    return('<div style="page-break-before: always;" />')
}


# pap_read_tif
# description: A function that takes a tif file representing where
# surfactant is present in an FOV and reads it into a tidy tibble
# a tibble using tidy format.
#
# input:
#
# output:
pap_read_tif <- function(filepath) {
  
  read_tif <- quietly(tiff::readTIFF)
  
  result <- 
    filepath %>% 
    read_tif(source = .) %>% 
    pluck("result") %>% 
    as_tibble() %>% 
    mutate(y = n():1) %>% 
    pivot_longer(
      cols = -y, 
      names_to = "x", 
      values_to = "values", 
      names_prefix = "V", 
    ) %>% 
    mutate(
      x = as.integer(x), 
      values = if_else(values == 1, "surfactant", "no surfactant") 
    ) %>% 
    select(x, y, values) %>% 
    arrange(x, y)
  
  return(result)
}

#pap_perform_daa()
# description: 
#
# input:
#
# output:
pap_perform_daa <- 
  function(
    data_tibble, 
    sample_col, 
    cluster_col,
    fixed_effect_cols,
    random_effect_cols,
    include_observation_level_random_effects = TRUE, 
    min_cells = 3, 
    min_samples = 5,
    ...
  ) { 
    
    # extract names from selected columns
    sample_colname <- 
      data_tibble %>% 
      select({{sample_col}}) %>% 
      colnames()
    
    fixed_effect_colnames <- 
      data_tibble %>% 
      select({{fixed_effect_cols}}) %>% 
      colnames()
    
    random_effect_colnames <- 
      data_tibble %>% 
      select({{random_effect_cols}}) %>% 
      colnames()
    
    marker_names <- 
      data_tibble %>% 
      select(
        -c({{sample_col}}, {{cluster_col}}),
        -any_of(fixed_effect_colnames), 
        -any_of(random_effect_colnames)
      ) %>% 
      colnames()
    
    # create experiment_info
    experiment_info <- 
      data_tibble %>% 
      select({{sample_col}}, {{fixed_effect_cols}}, {{random_effect_cols}}) %>% 
      distinct() %>% 
      rename(sample_id = {{sample_col}}) %>% 
      arrange(sample_id)
  
    # create marker_info
    marker_info <- 
      tibble(marker_name = marker_names) %>% 
      mutate(marker_class = "state")
    
    # create formula
    if (include_observation_level_random_effects) { 
      random_effect_colnames <- 
        c("sample_id", random_effect_colnames)
    }
    
    my_formula <- 
      createFormula(
        experiment_info = experiment_info,
        cols_fixed = fixed_effect_colnames,
        cols_random = random_effect_colnames
      )
    
    # create design matrix
    my_design <- 
      createDesignMatrix(
        experiment_info = experiment_info, 
        cols_design = fixed_effect_colnames
      )
    
    # make contrast matrix 
    my_contrast <- 
      createContrast(attr(my_design, "assign"))
    
    # configure data into the format diffcyt likes
    data_nested <- 
      data_tibble %>% 
      group_by({{sample_col}}) %>% 
      nest() %>% 
      arrange({{sample_col}}) %>% 
      pull(data)
    
    data_diff <- 
      prepareData(
        d_input = data_nested,
        experiment_info = as.data.frame(experiment_info),
        marker_info = as.data.frame(marker_info), 
        cols_to_include = 
          (
            colnames(select(data_tibble, -c({{sample_col}}))) %in% 
              marker_info$marker_name
          )
      )
    
    # add clusters to diffcyt object 
    temp <- 
      data_diff %>% 
      SummarizedExperiment::rowData()
    
    temp[,"cluster_id"] <- 
      data_tibble %>% 
      pull({{cluster_col}}) %>% 
      as.factor()
    
    SummarizedExperiment::rowData(data_diff) <- temp
    
    # fix some typing issues in the exprs component of the SummarizedExperiment
    data_exprs <- 
      data_diff %>% 
      SummarizedExperiment::assays() %>% 
      `[[`("exprs")
    
    data_colnames <- colnames(data_exprs)
    
    data_exprs <- 
      data_exprs %>% 
      apply(MARGIN = 2, FUN = as.numeric)
    
    colnames(data_exprs) <- data_colnames
    
    SummarizedExperiment::assays(data_diff)[["exprs"]] <- data_exprs
    
    # perform differential abundance testing using GLMMs
    cell_counts <- calcCounts(data_diff)
    
    da_results <- 
      testDA_GLMM(
        d_counts = cell_counts, 
        formula = my_formula, 
        contrast = my_contrast, 
        min_cells = min_cells, 
        min_samples = min_samples
      )
    
    da_results_limma <- NULL
      # testDA_voom(
      #   d_counts = cell_counts, 
      #   design = my_design, 
      #   contrast = my_contrast, 
      #   block_id = 
      #     experiment_info %>% 
      #     pull({{random_effect_cols}}) %>% 
      #     as.factor(), 
      #   min_cells = min_cells, 
      #   min_samples = min_samples
      # )
    
    # return result
    return(
      list(
        fixed_effect_colnames = fixed_effect_colnames, 
        random_effect_colnames = random_effect_colnames, 
        marker_names = marker_names, 
        experiment_info = experiment_info, 
        data_nested = data_nested,
        marker_info = marker_info,
        my_formula = my_formula, 
        my_design = my_design,
        my_contrast = my_contrast,
        data_diff = data_diff, 
        cell_counts = cell_counts,
        da_results = da_results, 
        da_results_limma = da_results_limma
      )
    )
    
  }

pap_perform_daa_manual <- 
  function(
    data_tibble, 
    sample_col, 
    cluster_col, 
    fixed_effect_cols, 
    random_effect_cols, 
    include_observation_level_random_effects = TRUE, 
    min_cells = 20, 
    min_samples = 5
  ) {
    NULL
  }

pap_perform_dea <- 
  function(
    data_tibble, 
    sample_col, 
    cluster_col, 
    fixed_effect_cols, 
    random_effect_col,
    min_cells = 20, 
    min_samples = 5, 
    ...
  ) { 
    
    # extract names from selected columns
    sample_colname <- 
      data_tibble %>% 
      select({{sample_col}}) %>% 
      colnames() 
    
    fixed_effect_colnames <- 
      data_tibble %>% 
      select({{fixed_effect_cols}}) %>% 
      colnames()
    
    random_effect_colnames <- 
      data_tibble %>% 
      select({{random_effect_col}}) %>% 
      colnames()
    
    marker_names <- 
      data_tibble %>% 
      select(
        -c({{sample_col}}, {{cluster_col}}),
        -any_of(fixed_effect_colnames), 
        -any_of(random_effect_colnames)
      ) %>% 
      colnames()
    
    # create experiment_info
    experiment_info <- 
      data_tibble %>% 
      select({{sample_col}}, {{fixed_effect_cols}}, {{random_effect_col}}) %>% 
      distinct() %>% 
      rename(sample_id = {{sample_col}}) %>% 
      arrange(sample_id)
    
    # create marker_info
    marker_info <- 
      tibble(marker_name = marker_names) %>% 
      mutate(marker_class = "state")
    
    # create formula
    my_formula <- 
      createFormula(
        experiment_info = experiment_info,
        cols_fixed = fixed_effect_colnames,
        cols_random = random_effect_colnames
      )
    
    # create design matrix
    my_design <- 
      createDesignMatrix(
        experiment_info = experiment_info, 
        cols_design = fixed_effect_colnames
      )
    
    # make contrast matrix 
    my_contrast <- 
      createContrast(attr(my_design, "assign"))
    
    # configure data into the format diffcyt likes
    data_nested <- 
      data_tibble %>% 
      group_by({{sample_col}}) %>% 
      nest() %>% 
      arrange({{sample_col}}) %>% 
      pull(data)
    
    data_diff <- 
      prepareData(
        d_input = data_nested,
        experiment_info = as.data.frame(experiment_info),
        marker_info = as.data.frame(marker_info), 
        cols_to_include = 
          (
            colnames(select(data_tibble, -c({{sample_col}}))) %in% 
              marker_info$marker_name
          )
      )
    
    # add clusters to diffcyt object 
    temp <- 
      data_diff %>% 
      SummarizedExperiment::rowData()
    
    temp[,"cluster_id"] <- 
      data_tibble %>% 
      pull({{cluster_col}}) %>% 
      as.factor()
    
    SummarizedExperiment::rowData(data_diff) <- temp
    
    # fix some typing issues in the exprs component of the SummarizedExperiment
    data_exprs <- 
      data_diff %>% 
      SummarizedExperiment::assays() %>% 
      `[[`("exprs")
    
    data_colnames <- colnames(data_exprs)
    
    data_exprs <- 
      data_exprs %>% 
      apply(MARGIN = 2, FUN = as.numeric)
    
    colnames(data_exprs) <- data_colnames
    
    SummarizedExperiment::assays(data_diff)[["exprs"]] <- data_exprs
    
    # perform differential expression testing using GLMMs
    cell_counts <- calcCounts(data_diff)
    cell_medians <- calcMedians(data_diff)
    
    de_results <- 
      testDS_limma(
        d_counts = cell_counts, 
        d_medians = cell_medians,
        design = my_design, 
        contrast = my_contrast, 
        block_id = 
          experiment_info %>% 
          pull({{random_effect_col}}) %>% 
          as.factor(),
        min_cells = min_cells, 
        min_samples = min_samples, 
        markers_to_test = rep(TRUE, nrow(marker_info)) 
        #weights = FALSE
      )
    
    # return result
    return(
      list(
        fixed_effect_colnames = fixed_effect_colnames, 
        random_effect_colnames = random_effect_colnames, 
        marker_names = marker_names, 
        experiment_info = experiment_info, 
        data_nested = data_nested,
        marker_info = marker_info,
        my_formula = my_formula, 
        my_design = my_design,
        my_contrast = my_contrast,
        data_diff = data_diff, 
        cell_medians = cell_medians,
        de_results = de_results 
      )
    )
  }

