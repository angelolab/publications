# MIBI_plot_pixel_cluster_data.R
# Created by: Erin McCaffrey
#
# Overview: This script takes the output of the pixel clustering data and 
# visualizes the resulting clusters as a heatmap.

library(pals)
library(gplots)

##..Step 1: Import data..##

meta_som_data <-read.csv('pixel_channel_avg_meta_cluster.csv')
som_data <-read.csv('pixel_channel_avg_som_cluster.csv')

##..Step 2: Plot a heatmap of the 100 clusters..##

som.m<-as.matrix(som_data[,2:18])
rownames(som.m)<-som_data$pixel_som_cluster

heatmap.2(t(som.m),
          scale = "row",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          col = rev(as.vector((brewer.rdbu(100)))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          breaks=seq(-3,3, length.out=101))

##..Step 3: Plot a heatmap of the metaclusters..##

metasom.m<-as.matrix(meta_som_data[,2:18])
rownames(metasom.m)<-meta_som_data$pixel_meta_cluster_rename

heatmap.2(t(metasom.m),
          scale = "row",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          col = rev(as.vector((brewer.rdbu(100)))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          breaks=seq(-3,3, length.out=101))

