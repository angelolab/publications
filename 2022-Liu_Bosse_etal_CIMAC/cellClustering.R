# Liu and Bosse et al., Reproducible, high-dimensional imaging in archival human tissue by Multiplexed Ion Beam Imaging by Time-of-Flight (MIBI-TOF)
# Cell clustering using mean expression of each cell, make heatmap in Figure 4

library(data.table)
library(FlowSOM)
library(pheatmap)
library(RColorBrewer)

set.seed(2021)

# Cell table is output from Mesmer, table with the data (and cluster assignments already appended) is available on Zenodo
seg_dat_path = "cell_table_size_normalized_clusters.csv"
# Markers to include in clustering
markers = c("CD3","CD4","CD8","CD11c","CD20","CD31","CD45","CD56","CD68","PANCK","PAX5","Vimentin")
# Z-score cap to use when hierarchical clustering
cap = 3
# Number of metaclusters
k = 20
# Blue-white-red colors for heatmaps
rwb_cols = colorRampPalette(c("royalblue4","white","red4"))(99)


## 1. Use a SOM to cluster cells into 100 clusters
# Read in size-normalized single cell data (output from Mesmer)
dat = fread(seg_dat_path)
# Remove rows with all 0's
dat = dat[rowSums(dat[,..markers])>0,]
# Frequency normalization
norm_dat = dat[,..markers]
norm_dat = norm_dat[,lapply(.SD, function(x) x/rowSums(norm_dat))]
# 99.9% normalization
norm_dat = norm_dat[,lapply(.SD,function(x) x/quantile(x,0.999))]
fSOM = SOM(as.matrix(norm_dat))
dat$cluster = fSOM$mapping[,1]
norm_dat$cluster = fSOM$mapping[,1]

# Get mean of each cluster
mean_dat = norm_dat[, lapply(.SD, mean), by = cluster, .SDcols=markers]
mat_dat = data.frame(mean_dat[, ..markers])
rownames(mat_dat) = paste0("clust_",mean_dat$cluster)
# Z-score the columns and cap
mat_dat = pmin(scale(mat_dat), cap)
# Determine breaks
breaks = seq(-cap, cap, length.out=100)
# Make heatmap
pdf("cellClustering_som_heatmaps.pdf",height=10,width=20)
p = pheatmap(t(mat_dat),
             color = rwb_cols,
             breaks = breaks)

# Get hierarchical clusters
clust_to_hclust = cutree(p$tree_col, k=k)
names(clust_to_hclust) = mean_dat$cluster
dat$hCluster_cap = clust_to_hclust[as.character(dat$cluster)]

# Get mean of hierarchical clusters
mean_dat = dat[, lapply(.SD, mean), by = hCluster_cap, .SDcols = markers]
mat_dat = data.frame(mean_dat[,..markers])
rownames(mat_dat) = paste0("clust_",mean_dat$hCluster_cap)
# Z-score the columns
mat_dat = scale(mat_dat)
# Determine breaks
breaks = seq(-cap, cap, length.out=100)
# Make heatmap
p = pheatmap(t(mat_dat),
             color = rwb_cols,
             breaks = breaks)
dev.off()


## 2. Using the heatmaps, manually annotate all 100 clusters (outside of R). Save clusters to file named "clust_to_pheno.csv" with 2 columns, 'cluster' with cluster numbers 1-100, and 'phenotype' with cell phenotypes for each cluster. Save data and make final heatmap.
annot = fread("clust_to_pheno.csv")
save_dat = annot[,c("cluster","phenotype")][dat, on=.(cluster)]
# Colors for heatmap (these are the phenotypes we identified in our clustering)
colors = data.table(phenotype_num = 1:9, 
                    phenotype=c("B_cell","Epithelial","Endothelial","T_cell","Macrophage","Fibroblast","DC","NK_cell","Other"),
                    color=c("#A0CBE8","#F1CE63","#8CD17D","#FF9D9A","#B07AA1","#499894","#9D7660","#FFBE7D","#808080"))
save_dat = colors[,c("phenotype_num","phenotype")][save_dat, on=.(phenotype)]
fwrite(save_dat, "cell_table_size_normalized_clusters.csv")

## Once data is annotated, make final heatmap
all_cell_colors = colors$color
names(all_cell_colors) = colors$phenotype
mat_colors = list(phenotype = all_cell_colors)

pdf("cellClustering_pheno_heatmap.pdf",height=8,width=12)
mean_dat = save_dat[,lapply(.SD,mean),by=phenotype,.SDcols=markers]
mat_dat = data.frame(mean_dat[,..markers])
rownames(mat_dat) = mean_dat$phenotype
# Z-score
mat_dat = scale(mat_dat)
# Annotations
mat_col = data.frame(phenotype = mean_dat$phenotype)
rownames(mat_col) = mean_dat$phenotype
# Make heatmap
breaks = seq(-3,3,length.out=100)
pheatmap(mat_dat,
         color = rwb_cols,
         breaks = breaks,
         annotation_row = mat_col,
         annotation_colors = mat_colors)
dev.off()

