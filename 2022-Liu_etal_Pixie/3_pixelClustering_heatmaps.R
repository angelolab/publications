# Make heatmaps of mean cluster expression
# Author: Candace Liu
# Date: 8/15/22

library(data.table)
library(pheatmap)
library(viridis)
library(RColorBrewer)

pixel_mat_dir = "pixel_mats" #where pixel matrices from imageToMatrix.py are saved
sigma = 2 #Gaussian blur sigma used in imageToMatrix.py
cap = 3 #hierarchical clustering cap
hclust_coln = "hCluster_cap" #column name of cluster ids

name = "sigma2_passes10_rep1" #name of file generated from pixelClustering.R
clust_path = paste0("pixelClustering_",name,"_clusters.csv")
clust_to_pheno_path = paste0("pixelClustering_",name,"_mapping.csv") #manual mapping file of each cluster to its phenotype
norm_path = "avg999_sigma2.csv" #99.9% normalization values, generated in pixelClustering.R 
colors_path = "px_colors.csv" #colors for each pixel cluster, one row for each cluster

# Markers
cluster_markers = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")

# Phenotype to color mapping
clust_to_pheno = fread(clust_to_pheno_path)
colors_tab = fread(colors_path)
clust_to_color = clust_to_pheno[colors_tab, on=.(phenotype)]
clust_to_color = clust_to_color[order(hCluster_cap)]

mat_colors = clust_to_color$color
names(mat_colors) = paste0("pixel_h",clust_to_color$hCluster_cap)
mat_colors = list(hclust = mat_colors)

# Blue-white-red colors
rwb_cols = colorRampPalette(c("royalblue4","white","red4"))(99)

# Pixel clustering clusters
clusters = fread(clust_path)

# Get all data
points = 1:12
dat_all_marks = rbindlist(lapply(points, function(x) {print(x)
                                                      one_tab = fread(file.path(pixel_mat_dir,paste0("Point", x, "_sigma", sigma,".csv")))
                                                      one_tab = one_tab[rowSums(one_tab[,..cluster_markers])>0, ] # remove empty pixels
                                                      return(one_tab)}))
# Combine with clusters
dat_all = dat_all_marks[clusters, on=.(sample,x,y)]
# Pixel normalization
dat = dat_all[,..cluster_markers]
dat = dat[,lapply(.SD, function(x) x/rowSums(dat))]
# 99.9% normalization
norm_vals = fread(norm_path)
norm_vals = norm_vals[,..cluster_markers]
dat = dat[,Map(`/`,.SD,norm_vals)]

dat$cluster = dat_all$cluster
dat$hCluster_cap = dat_all$hCluster_cap

## Clustering markers
pdf(paste0("pixelClustering_heatmaps_",name,".pdf"),height=8,width=8)
## Heatmap of pixel clusters x markers, average across pixel clusters
mean_dat = dat[, lapply(.SD, mean), by = cluster]
mean_dat = mean_dat[order(get(hclust_coln))]
mat_dat = data.frame(mean_dat[, ..cluster_markers])
rownames(mat_dat) = paste0("clust_", mean_dat$cluster)
# Z-score and cap
mat_dat = scale(mat_dat)
mat_dat = pmin(mat_dat, cap)
# Annotations
mat_col = data.frame(hclust = paste0("pixel_h", mean_dat[,get(hclust_coln)]))
rownames(mat_col) = paste0("clust_", mean_dat$cluster)
# Determine breaks
range = max(abs(mat_dat))
breaks = seq(-range, range, length.out=100)
# Make heatmap
pheatmap(mat_dat,
         color = rwb_cols,
         breaks = breaks,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         annotation_row = mat_col,
         annotation_colors = mat_colors,
         main = "Average across pixel clusters for markers used in clustering")

## Heatmap of pixel hierarchical cluster x markers, average across hierarchical clusters
mean_dat = dat[, lapply(.SD, mean), by = eval(hclust_coln), .SDcols = cluster_markers]
mat_dat = data.frame(mean_dat[,..cluster_markers])
rownames(mat_dat) = paste0("clust_",mean_dat[,get(hclust_coln)])
# Z-score the columns
mat_dat = scale(mat_dat)
# Make annotations
mat_col = data.frame(hclust = paste0("pixel_h",mean_dat[,get(hclust_coln)]))
rownames(mat_col) = paste0("clust_",mean_dat[,get(hclust_coln)])
# Determine breaks
range = 3
breaks = seq(-range, range, length.out=100)
# Make heatmap
pheatmap(mat_dat,
         color = rwb_cols,
         breaks = breaks,
         show_rownames = FALSE,
         annotation_row = mat_col,
         annotation_colors = mat_colors,
         main = "Average across pixel hierarchical clusters for markers used in clustering")
dev.off()

