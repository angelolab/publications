# PAH_MIBIheatmapDCbyCluster.R
# Author: Erin McCaffrey 
# Date created: 200731

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

##..Join the neighborhood and functional data..##

data_neighborhood_markers<-merge(data_neighborhood, data, by=c("Point_num","label",'PID'))

##..Filter to include just DCs..##

data_DCs<-droplevels(data_neighborhood_markers[data_neighborhood_markers$cell_lineage=='DC',])

##..Define clusters and markers to assess..##

functional_markers<-c('CD163','IDO1','iNOS','KI67','PDL1b','TIM.3','CD141')

clusters<-c(1,2,3,4,5,6,7,8,9,10)

##..Filter dataframe to only include clusters..##

data_DC_clusters<-droplevels(data_DCs[data_DCs$cluster %in% clusters, ])

##..Heatmap with mean expression..##

hm_neighborhood <- matrix(, nrow = length(clusters), ncol = length(functional_markers))

for(i in 1:length(clusters)) {
  temp_mat <- data_DC_clusters[data_DC_clusters$cluster == clusters[i], functional_markers]
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
          # col = diverging_hcl(100, palette = 'Blue-Red 3'),
          col = as.vector((cividis(100))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_neighborhood),
          rowsep=0:nrow(hm_neighborhood),
          cexRow=1,
          cexCol=1)

