# Cell clustering using integrated expression from cell table
# Author: Candace Liu
# Date: 8/15/22

library(RColorBrewer)
library(grid)
source("0_pixelClustering_functions.R")

## Variables to change
markers = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")
cap_cell = 3 #z-score cap for hierarchical clustering of cell clusters
k_cell = 14 #number of cell metaclusters
cell_filen = "cellClustering_integrated" #name to append to output
seeds = c(4,13,17,55,32)

# Colors for plotting
all_colors = c(brewer.pal(name="Set3", n=12), brewer.pal(name="Paired", n=12), brewer.pal(name="Pastel1", n=9), brewer.pal(name="Dark2", n=8), brewer.pal(name="Set1", n=9), brewer.pal(name="Pastel2", n=8))
all_cell_colors = all_colors[1:k_cell]
names(all_cell_colors) = paste0("cell_h",1:k_cell)

# Get data from cell segmentation (output from Mesmer)
seg_dat = fread("cell_table_size_normalized.csv")
# Remove cells with zero expression for clustering markers
seg_dat = seg_dat[rowSums(seg_dat[,..markers]) != 0]
# 99.9% normalization
dat_trans = seg_dat[,..markers]
dat_trans = dat_trans[, lapply(.SD, function(x) quantile_norm(x))]

# Function for one replicate
one_rep <- function(file_suffix, seed) {
  # Do cell clustering
  set.seed(seed)
  som_out = SOM(as.matrix(dat_trans[,..markers]))
  dat_trans$cluster = som_out$mapping[,1]

  ## Get consensus hierarchical clusters with no z-score capping
  # Get mean of each cluster
  mean_dat = dat_trans[, lapply(.SD, mean), by = cluster]
  mat_dat = data.frame(mean_dat[, ..markers])
  # Z-score the columns
  mat_dat = scale(mat_dat)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k_cell, seed=seed)
  clust_to_hclust = consensus[[k_cell]]$consensusClass
  names(clust_to_hclust) = mean_dat$cluster
  dat_trans$hCluster = clust_to_hclust[as.character(dat_trans$cluster)]

  ## Get hierarchical clusters with z-score capping
  mat_dat = pmin(mat_dat, cap_cell)
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k_cell, seed=seed)
  clust_to_hclust = consensus[[k_cell]]$consensusClass
  names(clust_to_hclust) = mean_dat$cluster
  dat_trans$hCluster_cap = clust_to_hclust[as.character(dat_trans$cluster)]

  out_dat = cbind(seg_dat[,c("sample","label","cell_size")], dat_trans[,c("cluster","hCluster","hCluster_cap")])
  fwrite(out_dat, paste0(file_suffix,".csv"))

  # Make cluster heatmap
  mean_dat = dat_trans[, lapply(.SD, mean), by=cluster]
  mean_dat = mean_dat[order(hCluster_cap)]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$cluster)
  # Z-score
  mat_dat = scale(mat_dat)
  mat_dat = pmin(mat_dat,cap_cell)
  # Annotations
  mat_col = data.frame(hclust_cap = paste0("cell_h",mean_dat$hCluster_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$cluster)
  mat_colors = list(hclust_cap = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$hCluster_cap)))])
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p1 = pheatmap(t(mat_dat),
                color = viridis(length(breaks)-1),
                breaks = breaks,
                cluster_cols = FALSE,
                show_colnames = FALSE,
                annotation_col = mat_col,
                annotation_colors = mat_colors,
                main = "Average across cell clusters")

  # Make hCluster heatmap
  mean_dat = dat_trans[, lapply(.SD, mean), .SDcols=markers, by=hCluster_cap]
  mean_dat = mean_dat[order(hCluster_cap)]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("cell_h",mean_dat$hCluster_cap)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust_cap = paste0("cell_h",mean_dat$hCluster_cap))
  rownames(mat_col) = paste0("cell_h",mean_dat$hCluster_cap)
  mat_colors = list(hclust_cap = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$hCluster_cap)))])
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p2 = pheatmap(t(mat_dat),
                color = viridis(length(breaks)-1),
                breaks = breaks,
                show_colnames = FALSE,
                annotation_col = mat_col,
                annotation_colors = mat_colors,
                main = "Average across cell hierarchical clusters")

  return(list(p1[[4]],p2[[4]]))
}

# Rep 1
cell_plots = one_rep(file_suffix=paste0(cell_filen,"_cellRep1"), seed=seeds[1])

# Rep 2
cell_plots = one_rep(file_suffix=paste0(cell_filen,"_cellRep2"), seed=seeds[2])

# Rep 3
cell_plots = one_rep(file_suffix=paste0(cell_filen,"_cellRep3"), seed=seeds[3])

# Rep 4
cell_plots = one_rep(file_suffix=paste0(cell_filen,"_cellRep4"), seed=seeds[4])

# Rep 5
cell_plots = one_rep(file_suffix=paste0(cell_filen,"_cellRep5"), seed=seeds[5])

