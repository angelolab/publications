# Cell clustering using pixel composition, run 5 replicates for cluster consistency score
# Author: Candace Liu
# Date: 8/15/22

library(RColorBrewer)
library(grid)
source("0_pixelClustering_functions.R")

## Variables to change
markers = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")
cap_cell = 3 #z-score cap for hierarchical clustering of cell clusters
k_cell = 14 #number of cell metaclusters
passes = 10 #passes for SOM
px_filen = "pixelClustering_sigma2_passes10_rep1_clusters.csv" #output of pixel clustering
px_mean_dat_path = "pixelClustering_sigma2_passes10_rep1_hCluster_mean.csv" #file that has mean expression of each pixel cluster
cell_filen = "passes10" #name to append to output
seeds = c(4,13,17,55,32)

# Get data from cell segmentation (output from Mesmer)
seg_dat = fread("cell_table_size_normalized.csv")
allPoints = unique(seg_dat$sample)

# Colors for plotting
all_colors = c(brewer.pal(name="Set3", n=12), brewer.pal(name="Paired", n=12), brewer.pal(name="Pastel1", n=9), brewer.pal(name="Dark2", n=8), brewer.pal(name="Set1", n=9), brewer.pal(name="Pastel2", n=8))


# Count number of pixels per cell
clusts = fread(px_filen)
clust_freq = pixelClusters_count(pixel_dat=clusts,
                                 coln="hCluster_cap",
                                 cell_dat=seg_dat[,c("sample","label","cell_size")],
                                 file_suffix=cell_filen,
                                 save_dat=FALSE)
# Get mean of each pixel cluster
pixel_mean_dat = fread(px_mean_dat_path)

# Function for one replicate
one_rep <- function(file_suffix, seed) {
  # Do cell clustering
  cell_clusts = cellClustering_flowsom_consensus(dat=clust_freq,
                                                 k_cell_hclust=k_cell,
                                                 cap=cap_cell,
                                                 quant_norm=0.999,
                                                 npasses=passes,
                                                 seed=seed,
                                                 file_suffix=file_suffix,
                                                 save_som=FALSE)
  # Make heatmaps
  cell_plots = cellClustering_plots(clust_freq=clust_freq,
                                    cell_clusts=cell_clusts,
                                    cell_seg_dat=seg_dat,
                                    pixel_mean_dat=pixel_mean_dat,
                                    pixel_dat_coln="hCluster_cap",
                                    markers=markers,
                                    cap=cap_cell,
                                    mat_colors=all_colors,
                                    cluster_rows=TRUE,
                                    show_rownames=TRUE,
                                    file_suffix=file_suffix)
  return(cell_plots)
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

