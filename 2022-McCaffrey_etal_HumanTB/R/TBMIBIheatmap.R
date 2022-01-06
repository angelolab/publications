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
data_gran<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Specify clustering channels and row order..##

major<-c('CD45','SMA','E.cadherin','CD31','Keratin.pan',"Vimentin")
immune<-c("CD3","CD20","CD14","CD16","CD11c","HLA.DR.DQ.DP","MastChyTry","MPO")
myeloid<-c("CD14","CD16","CD11c","CD11b","HLA.DR.DQ.DP","CD68","CD163","CD206","CD209")
lymph<-c("CD3","CD20","CD4","CD8","Foxp3","gdTCR")

heatmap_channels<-unique(c(major,immune,myeloid,lymph))

row_order<-c("CD4_T","B_cell","CD11b/c_CD206_Mac/Mono","CD8_T","CD14_Mono","CD11c_DC/Mono","imm_other","CD16_CD14_Mono","fibroblast","endothelial",
             "CD68_Mac","CD163_Mac","CD206_Mac","neutrophil","epithelial","CD209_DC","Treg","mast","gdT_cell","giant_cell")

##..Make heatmap matrix..##
                                                    
hm_allclusters <- matrix(, nrow = length(unique(data_gran$cell_type)), ncol = length(heatmap_channels))
for(i in 1:length(row_order)) {
  temp_mat <- data_gran[data_gran$cell_type==row_order[i], heatmap_channels]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

hm_allclusters
rownames(hm_allclusters) <- paste(row_order)
colnames(hm_allclusters) <- heatmap_channels

##..Get vector of color codes for annotation column..##

colorkey<-read.csv('data/colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% row_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = row_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color<-as.vector(colorkey_imm$code)

##..Make heatmap..##

heatmap.2(hm_allclusters,
          scale = "none",
          Colv = T, Rowv = F,
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
          RowSideColors = color,
          breaks=seq(0, 1, length.out=257))




