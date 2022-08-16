# Run pixel clustering
# 5 replicates for cluster consistency score calculation
# Author: Candace Liu
# Date: 8/15/22 

library(RColorBrewer)
library(grid)
library(gridExtra)
source("0_pixelClustering_functions.R")

## Variables to change
markers = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")
norm_path = "avg999_sigma2.csv" #path to save normalization values
pixel_mat_dir = "pixel_mats/" #where output of imageToMatrix.py is stored
sigma = 2 #sigma for Gaussian blur (must be same as used in imageToMatrix.py)
cap = 3 #z-score cap for hierarhcical clustering
k = 15 #number of pixel metaclusters
npasses = 10 #number of passes through dataset for FlowSOM
xdim = 10 #xdim for SOM grid
ydim = 10 #ydim for SOM grid
seeds = c(22,59,10,33,329) #seeds for replicates
name = "passes10" #appended to output name
save_dat = FALSE #whether to write full data table to file
save_som = FALSE #whether to save SOM object

# Get list of points
point_tab = fread("cell_table_size_normalized.csv") #cell table output from Mesmer
allPoints = unique(point_tab$sample)

# Get 99.9% normalization values
one_point <- function(point_num) {
  # Read data
  one_tab = fread(file.path(pixel_mat_dir,paste0("Point", point_num, "_sigma", sigma,".csv")))
  one_tab = one_tab[rowSums(one_tab[,..markers])>0, ] # remove empty pixels

  # Pixel normalization
  tab_markers = one_tab[,..markers]
  dat = tab_markers[,lapply(.SD, function(x) x/rowSums(tab_markers))]

  # Get 99.9th percentile
  all_quants = dat[,lapply(.SD, function(x) get_quantile_norm(x))]

  return(all_quants)
}
all_point_quant = rbindlist(lapply(allPoints, function(x) one_point(x)))
norm_vals = all_point_quant[,lapply(.SD,mean), .SDcols=markers]
fwrite(norm_vals, norm_path)


# Colors for plotting
all_colors = c(brewer.pal(name="Set3", n=12), brewer.pal(name="Paired", n=12), brewer.pal(name="Pastel1", n=9), brewer.pal(name="Dark2", n=8), brewer.pal(name="Set1", n=9), brewer.pal(name="Pastel2", n=8))


# Replicate 1
pixel_dat = pixelClustering_flowsom_consensus(points=allPoints,
                                              markers=markers,
                                              norm_vals=norm_vals,
                                              pixel_mat_dir=pixel_mat_dir,
                                              sigma=sigma,
                                              k=k,
                                              cap=cap,
                                              npasses=npasses,
                                              xdim=xdim,
                                              ydim=ydim,
                                              seed=seeds[1],
                                              file_suffix=paste0(name,"_rep1"),
                                              save_dat=save_dat,
                                              save_som=save_som)
all_plots = pixelClustering_plots(pixel_dat=pixel_dat,
                                  hclust_coln = "hCluster_cap",
                                  colors=all_colors,
                                  markers=markers,
                                  cap=cap)

# Replicate 2
pixel_dat = pixelClustering_flowsom_consensus(points=allPoints,
                                              markers=markers,
                                              norm_vals=norm_vals,
                                              pixel_mat_dir=pixel_mat_dir,
                                              sigma=sigma,
                                              k=k,
                                              cap=cap,
                                              npasses=npasses,
                                              xdim=xdim,
                                              ydim=ydim,
                                              seed=seeds[2],
                                              file_suffix=paste0(name,"_rep2"),
                                              save_dat=save_dat,
                                              save_som=save_som)
all_plots = pixelClustering_plots(pixel_dat=pixel_dat,
                                  hclust_coln = "hCluster_cap",
                                  colors=all_colors,
                                  markers=markers,
                                  cap=cap)

# Replicate 3
pixel_dat = pixelClustering_flowsom_consensus(points=allPoints,
                                              markers=markers,
                                              norm_vals=norm_vals,
                                              pixel_mat_dir=pixel_mat_dir,
                                              sigma=sigma,
                                              k=k,
                                              cap=cap,
                                              npasses=npasses,
                                              xdim=xdim,
                                              ydim=ydim,
                                              seed=seeds[3],
                                              file_suffix=paste0(name,"_rep3"),
                                              save_dat=save_dat,
                                              save_som=save_som)
all_plots = pixelClustering_plots(pixel_dat=pixel_dat,
                                  hclust_coln = "hCluster_cap",
                                  colors=all_colors,
                                  markers=markers,
                                  cap=cap)

# Replicate 4
pixel_dat = pixelClustering_flowsom_consensus(points=allPoints,
                                              markers=markers,
                                              norm_vals=norm_vals,
                                              pixel_mat_dir=pixel_mat_dir,
                                              sigma=sigma,
                                              k=k,
                                              cap=cap,
                                              npasses=npasses,
                                              xdim=xdim,
                                              ydim=ydim,
                                              seed=seeds[4],
                                              file_suffix=paste0(name,"_rep4"),
                                              save_dat=save_dat,
                                              save_som=save_som)
all_plots = pixelClustering_plots(pixel_dat=pixel_dat,
                                  hclust_coln = "hCluster_cap",
                                  colors=all_colors,
                                  markers=markers,
                                  cap=cap)

# Replicate 5
pixel_dat = pixelClustering_flowsom_consensus(points=allPoints,
                                              markers=markers,
                                              norm_vals=norm_vals,
                                              pixel_mat_dir=pixel_mat_dir,
                                              sigma=sigma,
                                              k=k,
                                              cap=cap,
                                              npasses=npasses,
                                              xdim=xdim,
                                              ydim=ydim,
                                              seed=seeds[5],
                                              file_suffix=paste0(name,"_rep5"),
                                              save_dat=save_dat,
                                              save_som=save_som)
all_plots = pixelClustering_plots(pixel_dat=pixel_dat,
                                  hclust_coln = "hCluster_cap",
                                  colors=all_colors,
                                  markers=markers,
                                  cap=cap)


