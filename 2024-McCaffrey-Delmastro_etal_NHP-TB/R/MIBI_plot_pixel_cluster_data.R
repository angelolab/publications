# MIBI_plot_pixel_cluster_data.R

library(pals)
library(gplots)

##..Import data..##

# Read in necrotic and non-necrotic separately
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/NHP_TB/Cohort/downstream_analysis/phenotyping/Panel2/panel2_pixel_output_dir")
meta_som_data <-read.csv('pixel_channel_avg_meta_cluster.csv')
som_data <-read.csv('pixel_channel_avg_som_cluster.csv')

# Heat map of the 100 clusters

som.m<-as.matrix(som_data[,2:18])
rownames(som.m)<-som_data$pixel_som_cluster

heatmap.2(t(som.m),
          scale = "row",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          # col = magma(256),
          col = rev(as.vector((brewer.rdbu(100)))),
          # col = (as.vector((ocean.delta(100)))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          breaks=seq(-3,3, length.out=101))

# Heat map of metaclust

metasom.m<-as.matrix(meta_som_data[,2:18])
rownames(metasom.m)<-meta_som_data$pixel_meta_cluster_rename

heatmap.2(t(metasom.m),
          scale = "row",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          # col = magma(256),
          col = rev(as.vector((brewer.rdbu(100)))),
          # col = (as.vector((ocean.delta(100)))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          breaks=seq(-3,3, length.out=101))

