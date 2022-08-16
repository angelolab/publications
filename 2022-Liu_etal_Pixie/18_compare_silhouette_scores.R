# Compare Silhouette scores of cell clustering using pixel composition or integrated expression
# Author: Candace Liu
# Dat: 8/16/22

library(data.table)
library(cluster)
library(ggplot2)
library(ggpubr)

## Pixel cluster composition
reps = 1:5 #all replicates
hclust_coln = "pixelfreq_hclust_cap" #name of column with cluster ids
coln = paste0("pixel_",1:14) #column names of feature table 
feature_dat = fread("freq_px_passes10_rep5.csv") #feature table for all cells
one_rep <- function(r) {
  dat = fread(paste0("cellClustering_kcell14_passes10_pixelComposition_cellRep",r,".csv"))
  dat = feature_dat[dat, on=.(sample,label)]
  ss = silhouette(dat[,get(hclust_coln)], dist(dat[,..coln]/dat$cell_size))
  return(mean(ss[,c("sil_width")]))
}
all_reps = lapply(reps,one_rep)
pxclustfreq_dt = data.table(group="pxclustfreq",ss=unlist(all_reps))

## Integrated expression
reps = 1:5 #all replicates
hclust_coln = "hCluster_cap" #name of column with cluster ids
coln = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")
feature_dat = fread("cell_table_size_normalized.csv") #cell table output from Mesmer
one_rep <- function(r) {
  dat = fread(paste0("cellClustering_integrated_cellRep",r,".csv"))
  dat = feature_dat[dat, on=.(sample,label)]
  ss = silhouette(dat[,get(hclust_coln)], dist(dat[,..coln]))
  return(mean(ss[,c("sil_width")]))
}
all_reps = lapply(reps,one_rep)
seg_dt = data.table(group="seg",ss=unlist(all_reps))

# All data
all_dt = rbindlist(list(pxclustfreq_dt,seg_dt))

# Plot
all_dt$group = factor(all_dt$group, levels=c("seg","pxclustfreq"))
ggplot(all_dt, aes(x=group,y=ss,fill="lightgray")) +
  geom_jitter(width=0.05, height=0) +
  stat_summary(fun=mean, geom="point", color="red", shape="-", size=10) +
  ylim(c(-0.2,0.5)) +
  xlab("") +
  ylab("Silhouette score") +
  theme_light() +
  theme(legend.position="none")

# Get p-value
compare_means(ss ~ group, data=all_dt, method="wilcox.test", paired=FALSE)

