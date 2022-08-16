# Test reproducibility of clustering using CyTOF data
# Data downloaded from https://doi.org/10.5281/zenodo.3951613
# Author: Candace Liu
# Date: 8/15/22

library(data.table)
library(FlowSOM)
library(flowCore)
library(Seurat)
library(ConsensusClusterPlus)

markers = c("CD45","CD3","CD4","CD8","CD45RA","CD66","CD14","CD19","CD20","HLADR","CD56","CD57","CD11c","CD123","FceRI","CD235ab")
antibody_mapping = fread("antibody_mapping.csv") #obtained from publication (Hartmann et al, 2020)
antibody_to_metal = antibody_mapping$marker
names(antibody_to_metal) = antibody_mapping$metal
metals = antibody_mapping$metal

fcs = read.FCS("Experiment_1/e0b3_wb_exp1_CD45.fcs", transformation=FALSE, emptyValue=FALSE)
dt = data.table(exprs(fcs))
dt_markers = dt[,..metals]
colnames(dt_markers) = antibody_to_metal[colnames(dt_markers)]

# Randomly subset 5000 cells from dataset
set.seed(329)
rand_samp = sample(1:nrow(dt_markers), 5000)
dt_sub = dt_markers[rand_samp]


## Test reproducibilty of clustering
all_seeds = c(2022,32,115,6,22)
k = 15 #total number of clusters

one_seed <- function(seed) {
  set.seed(seed)
  output = SOM(as.matrix(dt_sub))
  dt_sub$cluster = output$mapping[,1]
  mean_dat = dt_sub[,lapply(.SD,mean), by=cluster]
  mat_dat = data.frame(mean_dat[,..markers])
  # Z-score the columns
  mat_dat = scale(mat_dat)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k, seed=seed)
  clust_to_hclust = consensus[[k]]$consensusClass
  names(clust_to_hclust) = mean_dat$cluster
  dt_sub$hCluster = clust_to_hclust[as.character(dt_sub$cluster)]
  return(dt_sub$hCluster)
}
all_seed_output = lapply(all_seeds, one_seed)
dt = data.table(rep1=as.integer(all_seed_output[[1]]),
                rep2=as.integer(all_seed_output[[2]]),
                rep3=as.integer(all_seed_output[[3]]),
                rep4=as.integer(all_seed_output[[4]]),
                rep5=as.integer(all_seed_output[[5]]))
fwrite(dt,"cytof_flowsom.csv")


