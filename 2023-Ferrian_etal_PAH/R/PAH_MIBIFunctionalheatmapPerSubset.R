# PAH_MIBIFunctionalheatmapPerSubset.R
# Author: Erin McCaffrey 
# Date created: 200310
# Overview: This script reads in the csv for all PAH cell data and produces a heatmap of the mean expression
# for all phenotypic markers across the cell populations and a second heatmap with expression of key functional
# markers across all populations in healthy versus PAH to identify any expression patterns that vary between 
# the groups. 

require(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(viridis)
library(colorspace)
library(pals)


##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("allpatients_freq-of-positivecells-perregion.csv")

##..Make separate healthy and PAH matrices for functional markers..##

hlt<-droplevels(data[data$Group=='hlt',])
pah<-droplevels(data[data$Group =='pah',])
ipah<-droplevels(data[data$Subgroup=='ipah',])
hpah<-droplevels(data[data$Subgroup=='hpah',])


##..Define markers and lineages..##

functional_markers<-colnames(data[,c(4,7,9,12,19,23,24,25,26,27,28,29,30,34,35,37)])

regions<-c('vessel_associated', 'vessel_proximal', 'non_vascular')

color<-c('#651E38','#C57050','#C8C8C7')

##..Generate heatmap of mean marker expression per subset in PAH and Hlt..##

hm_hlt <- matrix(, nrow = length(regions), ncol = length(functional_markers))
hm_pah <- hm_hlt
hm_ipah <- hm_hlt
hm_hpah <- hm_hlt

for(i in 1:length(regions)) {
  #hlt
  temp_mat_hlt <- hlt[hlt$region_modified==regions[i], functional_markers]
  hm_hlt[i,] <- as.matrix(apply(temp_mat_hlt, 2, function (x) mean(x, na.rm = T)))
  #pah
  temp_mat_pah <- pah[pah$region_modified==regions[i], functional_markers]
  hm_pah[i,] <- as.matrix(apply(temp_mat_pah, 2, function (x) mean(x, na.rm = T)))
  #ipah
  temp_mat_ipah <- ipah[ipah$region_modified==regions[i], functional_markers]
  hm_ipah[i,] <- as.matrix(apply(temp_mat_ipah, 2, function (x) mean(x, na.rm = T)))
  #hpah
  temp_mat_hpah <- hpah[hpah$region_modified==regions[i], functional_markers]
  hm_hpah[i,] <- as.matrix(apply(temp_mat_hpah, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_hlt) <- paste(regions)
colnames(hm_hlt) <- functional_markers

rownames(hm_pah) <- rownames(hm_hlt)
colnames(hm_pah) <- colnames(hm_hlt)

rownames(hm_ipah) <- rownames(hm_hlt)
colnames(hm_ipah) <- colnames(hm_hlt)

rownames(hm_hpah) <- rownames(hm_hlt)
colnames(hm_hpah) <- colnames(hm_hlt)


hm.hlt<-heatmap.2(hm_hlt,
                  scale = "none",
                  Colv = F, Rowv = F,
                  dendrogram = "none",
                  trace = "none",
                  col = magma(256),
                  density.info = 'none',
                  key.title = '',
                  sepwidth=c(0.01,0.01),
                  sepcolor="grey",
                  colsep=0:ncol(hm_hlt),
                  rowsep=0:nrow(hm_hlt),
                  cexRow=0.4,
                  cexCol=0.4,
                  RowSideColors = color,
                  breaks = seq(0, 0.5, length.out = 257))

hm.pah<-heatmap.2(hm_pah,
                  scale = "none",
                  Colv = F, Rowv = F,
                  dendrogram = "none",
                  trace = "none",
                  col = magma(256),
                  density.info = 'none',
                  key.title = '',
                  sepwidth=c(0.01,0.01),
                  sepcolor="grey",
                  colsep=0:ncol(hm_pah),
                  rowsep=0:nrow(hm_pah),
                  cexRow=0.4,
                  cexCol=0.4,
                  RowSideColors = color,
                  breaks = seq(0, 0.5, length.out = 257))

hm.ipah<-heatmap.2(hm_ipah,
                   scale = "none",
                   Colv = F, Rowv = F,
                   dendrogram = "none",
                   trace = "none",
                   col = magma(256),
                   density.info = 'none',
                   key.title = '',
                   sepwidth=c(0.01,0.01),
                   sepcolor="grey",
                   colsep=0:ncol(hm_ipah),
                   rowsep=0:nrow(hm_ipah),
                   cexRow=0.4,
                   cexCol=0.4,
                   RowSideColors = color,
                   breaks = seq(0, 0.5, length.out = 257))

hm.hpah<-heatmap.2(hm_hpah,
                   scale = "none",
                   Colv = F, Rowv = F,
                   dendrogram = "none",
                   trace = "none",
                   col = magma(256),
                   density.info = 'none',
                   key.title = '',
                   sepwidth=c(0.01,0.01),
                   sepcolor="grey",
                   colsep=0:ncol(hm_hpah),
                   rowsep=0:nrow(hm_hpah),
                   cexRow=0.4,
                   cexCol=0.4,
                   RowSideColors = color,
                   breaks = seq(0, 0.5, length.out = 257))