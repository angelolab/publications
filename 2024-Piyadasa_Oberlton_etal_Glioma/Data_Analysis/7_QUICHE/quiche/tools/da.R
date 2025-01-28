suppressPackageStartupMessages(
    {
        library(SingleCellExperiment)
        library(scran)
        library(DAseq)
        library(tibble)
        library(dplyr)
        library(igraph)
        library(cydar)
        library(pdist)
        library(reshape2)
    }
)
run_daseq <- function(
    sce,
    k.vec,
    condition_col,
    sample_col="synth_samples",
    reduced.dim = "PCA",
    d=30
){
  # prepare the params needed by func: getDAcells
  condition_vec <- colData(sce)[[condition_col]]
  conds <- levels(factor(condition_vec))
  cell.labels <- sce[[sample_col]]
  if (length(conds) != 2){
    stop(str_c("daseq can only handle binary labels, but got ", length(conds)))
  }
  # main routine to get DA cells
  daseq_res <- getDAcells(X=reducedDim(sce, reduced.dim)[,1:d],
                          cell.labels=cell.labels,
                          labels.1=unique(cell.labels[condition_vec == conds[1]]),
                          labels.2=unique(cell.labels[condition_vec == conds[2]]),
                          k.vector=k.vec,
                          size=1,
                          do.plot=FALSE)
  return(daseq_res)
}

run_cydar <- function(
    sce,
    condition_col="synth_labels",
    sample_col="synth_samples",
    reduced.dim="pca.corrected",
    d=30,
    batch_col=NULL,
    tol=0.5,
    downsample=10,
    returnCd=TRUE
){
  ## Make the design matrix as Milo
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Make list for each sample
  sample_ls <- split(1:ncol(sce), sce[[sample_col]])
  processed.exprs <- lapply(sample_ls, function(s) reducedDim(sce[,s], reduced.dim)[,1:d])
  cd <- prepareCellData(processed.exprs)
  ## Count cells in hyperspheres
  cd <- cydar::countCells(cd, tol=tol, filter=1, downsample=downsample)
  # do DA testing with edgeR
  cd.dge <- DGEList(assay(cd), lib.size=cd$totals)
  
  sim.design <- model.matrix(design, data=design_df)[colnames(cd.dge),]
  sim.dge <- estimateDisp(cd.dge, sim.design)
  sim.fit <- glmQLFit(sim.dge, sim.design)
  sim.res <- glmQLFTest(sim.fit, coef=2)
  
  # control the spatial FDR
  cydar.res <- sim.res$table
  cydar.res$SpatialFDR <- spatialFDR(intensities(cd), sim.res$table$PValue)
  if (returnCd) {
    list(Cd=cd, DAres=cydar.res)
  } else {
    cydar.res
  }
}