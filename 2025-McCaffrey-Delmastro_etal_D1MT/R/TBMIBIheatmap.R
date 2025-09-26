# TBMIBIheatmap.R
# Author: Erin McCaffrey 
# Date created: 190204
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it produces a heat map ordered by cell-type and with the normalized expression of all clustering channels plus
# others that were not included in clustering, but important. Columns are arranged with heirarchical clustering. It can
# also produce a heatmap where each cell is a row. 

require(dplyr)
library(gplots)
library(viridis)

##..Import data..##
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data_gran<-read.csv('NHP_cohort_data_norm_annotated.csv')

##..Drop excluded samples..##
drop_samples <- c(17,22,26)
data_gran <- data_gran[!data_gran$SampleID %in% drop_samples,]

##..Specify clustering channels and row order..##
major<-c('CD45','SMA','CD31','Keratin',"VIM")
immune<-c("CD3","CD20","CD14","CD16","CD11c","HLA.DR","MastChyTry","Calprotectin")
myeloid<-c("CD14","CD16","CD11c","CD11b","CD68","CD163","CD206")
lymph<-c("CD3","CD20","CD4","CD8","Foxp3")
functional<-c("IDO","PDL1","GrzB","IFNg","H3K9Ac","H3K27me3","pS6","CD40","Ki67")

heatmap_channels<-unique(c(major,immune,myeloid,lymph))

##..Make heatmap matrix..##
row_order = unique(data_gran$cell_type)                                    
hm_allclusters <- matrix(, nrow = length(unique(data_gran$cell_type)), ncol = length(heatmap_channels))
for(i in 1:length(row_order)) {
  temp_mat <- data_gran[data_gran$cell_type==row_order[i], heatmap_channels]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

hm_allclusters
rownames(hm_allclusters) <- paste(row_order)
colnames(hm_allclusters) <- heatmap_channels

##..Make heatmap..##

heatmap.2(hm_allclusters,
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "column",
          trace = "none",
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="black",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          cexRow=0.4,
          cexCol=0.4,
          breaks=seq(0, 1, length.out=257))




