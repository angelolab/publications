# PAH_MIBIheatmapNeighborhood.R
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
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")
data_neighborhood<-read.csv("PAH_AllCells_Neighborhood_K=10_PatientAnnotated.csv")

##..Produce a heatmap of the cell frequencies..##

# go through all clusters and get mean frequency of each cell type

clusters<-c(1,2,3,4,5,6,7,8,9,10)
cell_types<-colnames(data_neighborhood[,3:14])
hm_allclusters <- matrix(, nrow = length(clusters), ncol = length(cell_types))
for(i in 1:length(clusters)) {
  temp_mat <- data_neighborhood[data_neighborhood[,"cluster"] == clusters[i], cell_types]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- clusters
colnames(hm_allclusters) <- cell_types  
hm_allclusters

# plot heatmap of all clusters

colorkey_neighborhood<-read.csv('neighborhood_colorkey_1.csv')
color<-as.vector(colorkey_neighborhood$code)

heatmap.2(hm_allclusters, 
          Colv = T, Rowv = F,
          hclustfun = hclust,
          scale = "column",
          dendrogram = c("column"),
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          cexRow=1,
          cexCol=1,
          RowSideColors = color)


##..Merge the neighborhood data with the expression data..##

data_neighborhood_markers<-merge(data_neighborhood, data, by=c("Point_num","label",'PID'))

##..Define markers and clusters..##

functional_markers<-c('CD45RO','GranzB','CD163','IDO1','IFNg','iNOS','KI67',
                      'PDL1b','TIM.3','MPO')

clusters<-c(1,2,3,4,5,6,7,8,9,10)

##..Create heatmap for all PAH..##

hm_neighborhood <- matrix(, nrow = length(clusters), ncol = length(functional_markers))

for(i in 1:length(clusters)) {
  temp_mat <- data_neighborhood_markers[data_neighborhood_markers$cluster == clusters[i], functional_markers]
  hm_neighborhood[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_neighborhood)<-clusters
colnames(hm_neighborhood)<-functional_markers

heatmap.2(hm_neighborhood,
          scale = "column",
          Colv = T, Rowv = F,
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "column",
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_neighborhood),
          rowsep=0:nrow(hm_neighborhood),
          cexRow=1,
          cexCol=1,
          RowSideColors = color)


##..Create heatmap for IPAH..##

data_neighborhood_markers_ipah<-droplevels(data_neighborhood_markers[data_neighborhood_markers$Subgroup =='IPAH',])

hm_neighborhood <- matrix(, nrow = length(clusters), ncol = length(functional_markers))

for(i in 1:length(clusters)) {
  temp_mat <- data_neighborhood_markers_ipah[data_neighborhood_markers_ipah$cluster == clusters[i], functional_markers]
  hm_neighborhood[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_neighborhood)<-clusters
colnames(hm_neighborhood)<-functional_markers

colorkey_neighborhood<-read.csv('neighborhood_colorkey_1.csv')
color<-as.vector(colorkey_neighborhood$code)

heatmap.2(hm_neighborhood,
          scale = "none",
          Colv = T, Rowv = F,
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "column",
          trace = "none",
          # col = diverging_hcl(100, palette = 'Blue-Red 3'),
          col = as.vector((ocean.thermal(100))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_neighborhood),
          rowsep=0:nrow(hm_neighborhood),
          cexRow=1,
          cexCol=1,
          RowSideColors = color)


##..Create heatmap for HPAH..##

data_neighborhood_markers_hpah<-droplevels(data_neighborhood_markers[data_neighborhood_markers$Subgroup =='HPAH',])

hm_neighborhood <- matrix(, nrow = length(clusters), ncol = length(functional_markers))

for(i in 1:length(clusters)) {
  temp_mat <- data_neighborhood_markers_hpah[data_neighborhood_markers_hpah$cluster == clusters[i], functional_markers]
  hm_neighborhood[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_neighborhood)<-clusters
colnames(hm_neighborhood)<-functional_markers

colorkey_neighborhood<-read.csv('neighborhood_colorkey_1.csv')
color<-as.vector(colorkey_neighborhood$code)

heatmap.2(hm_neighborhood,
          scale = "column",
          Colv = T, Rowv = F,
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "column",
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_neighborhood),
          rowsep=0:nrow(hm_neighborhood),
          cexRow=1,
          cexCol=1,
          RowSideColors = color)
