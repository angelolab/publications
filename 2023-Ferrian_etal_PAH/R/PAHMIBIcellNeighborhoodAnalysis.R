# PAHMIBIcellNeighborhoodAnalysis.R
# Author: Erin McCaffrey
# Date created: 191201
# Overview: This script reads in the csv for the cell neighborhood frequency data. 
# It separates PAH and healthy. It summarizes neighborhood clusters for immune cells.

library(dplyr)
library(flowCore)           
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(colorspace)

# set working directory and read in data
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
neighbor_data<-read.csv('cell_neighbor_counts_50px.csv')
cell_data<-read.csv("celldata_region_annotated.csv")

# convert to freqs
neighbor_freqs <- neighbor_data
neighbor_freqs[,3:14] <- neighbor_freqs[,3:14]/rowSums(neighbor_freqs[,3:14])
neighbor_freqs <- replace(neighbor_freqs, is.na(neighbor_freqs), 0)

# get just the immune cell neighborhoods

neighbor_freqs$cell_type<-cell_data$cell_lineage
immune_neighbor_freqs<-droplevels(neighbor_freqs[!neighbor_freqs$cell_type %in% 
                                                   c('fibroblast','endothelial','mesenchymal','epithelial'),])

# separate healthy and ipah+hpah
ipah<-unique(cell_data[cell_data$Subgroup=='IPAH',]$Point_num)
hpah<-unique(cell_data[cell_data$Subgroup=='HPAH',]$Point_num)
hlt<-unique(cell_data[cell_data$Subgroup=='Healthy Controls',]$Point_num)

neighbor_pah<-droplevels(neighbor_freqs[neighbor_freqs$Point_num %in% c(ipah, hpah),])
neighbor_hlt<-droplevels(neighbor_freqs[neighbor_freqs$Point_num %in% hlt,])

immune_neighbor_pah<-droplevels(immune_neighbor_freqs[immune_neighbor_freqs$Point_num %in% c(ipah, hpah),])
immune_neighbor_hlt<-droplevels(immune_neighbor_freqs[immune_neighbor_freqs$Point_num %in% hlt,])

######......PAH......######

## Visualize and hierarchically cluster to assess structure/ good k for k-means

hm<-as.matrix(neighbor_pah[,3:14])
hm_immune<-as.matrix(immune_neighbor_pah[,3:14])

heatmap.2(hm, 
          Colv = T, Rowv = T,
          hclustfun = hclust,
          scale = "row",
          dendrogram = c("none"),
          trace = "none",
          #col = viridis(256),
          col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
          density.info = 'none',
          key.title = '',
          cexRow = 0.6, cexCol = 0.9)

##..K means cluster the data..##

## k-means cluster the TB data
pah_clusters<-kmeans(immune_neighbor_pah[,3:14], 10)

# append cluster annotation to TB data
immune_neighbor_pah$cluster <- pah_clusters$cluster

# go through all clusters and get mean frequency of each cell type
clusters<-unique(immune_neighbor_pah$cluster)
cell_types<-colnames(immune_neighbor_pah[,3:14])
hm_allclusters <- matrix(, nrow = length(clusters), ncol = length(cell_types))
for(i in 1:length(clusters)) {
  temp_mat <- immune_neighbor_pah[immune_neighbor_pah[,"cluster"] == clusters[i], cell_types]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- clusters
colnames(hm_allclusters) <- cell_types  
hm_allclusters

# plot heatmap of all clusters
dev.off()
heatmap.2(hm_allclusters, 
          Colv = T, Rowv = T,
          hclustfun = hclust,
          scale = "column",
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          # col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
          density.info = 'none',
          key.title = '',
          cexRow = 0.6, cexCol = 0.9)

##...Save results..##
write.csv(neighbor_pah,"PAH_AllCells_Neighborhood_K=10.csv", row.names = FALSE)

## See how many clusters phenograph pulls out

library(Rphenograph)
phenograph<-Rphenograph(as.matrix(neighbor_pah[,3:14]))
pheno_clusters<-as.data.frame(factor(membership(phenograph[[2]])))
colnames(pheno_clusters)<-'phenograph'

neighbor_pah$pheno_cluster <- pheno_clusters$phenograph

# produce heatmap
clusters<-unique(pheno_clusters$phenograph)
cell_types<-colnames(neighbor_pah[,3:14])
hm_allclusters <- matrix(, nrow = length(clusters), ncol = length(cell_types))
for(i in 1:length(clusters)) {
  temp_mat <- neighbor_pah[neighbor_pah[,"pheno_cluster"] == clusters[i], cell_types]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- clusters
colnames(hm_allclusters) <- cell_types  
hm_allclusters

# plot heatmap of all clusters
dev.off()
heatmap.2(hm_allclusters, 
          Colv = T, Rowv = T,
          hclustfun = hclust,
          scale = "row",
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          # col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
          density.info = 'none',
          key.title = '',
          cexRow = 0.6, cexCol = 0.9)

# #meta-cluster the phenograph clusters based on mean frequency
# 
# pheno_clusters<-as.data.frame(hm_allclusters)
# 
# # determine K
# library(factoextra)
# library(NbClust)
# 
# #source: https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
# 
df<-neighbor_pah[,3:14]
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 10, linetype = 2)+
  labs(subtitle = "Elbow method")
# 
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# 
# set.seed(123)
# fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50, k.max = 10)+
#   labs(subtitle = "Gap statistic method")
# 
# 
# # K means cluster the phenograph clusters
# pheno_metaclusters<-kmeans(pheno_clusters, 10)
# 
# # append cluster annotation to data
# pheno_clusters$cluster <- pheno_metaclusters$cluster
# 
# # append cluster assignment to all immune cells in neighbor freq matrix
# 
# pheno_clusters$phenograph<-row.names(pheno_clusters)
# pheno_clusters$phenograph<- lapply(pheno_clusters$phenograph, gsub, 
#                                    pattern = "cluster", replacement = "", fixed = TRUE)
# kclust1<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==1,]$phenograph))
# kclust2<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==2,]$phenograph))
# kclust3<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==3,]$phenograph))
# kclust4<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==4,]$phenograph))
# kclust5<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==5,]$phenograph))
# kclust6<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==6,]$phenograph))
# kclust7<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==7,]$phenograph))
# kclust8<-as.numeric(unique(pheno_clusters[pheno_clusters$cluster==8,]$phenograph))
# 
# immune_neighbor_pah$pheno_meta<-immune_neighbor_pah$pheno_cluster
# immune_neighbor_pah<- immune_neighbor_pah %>% 
#   mutate(pheno_meta=case_when(immune_neighbor_pah$pheno_meta %in% kclust1~ 1,
#                               immune_neighbor_pah$pheno_meta %in% kclust2~ 2,
#                               immune_neighbor_pah$pheno_meta %in% kclust3~ 3,
#                               immune_neighbor_pah$pheno_meta %in% kclust4~ 4,
#                               immune_neighbor_pah$pheno_meta %in% kclust5~ 5,
#                               immune_neighbor_pah$pheno_meta %in% kclust6~ 6,
#                               immune_neighbor_pah$pheno_meta %in% kclust7~ 7,
#                               immune_neighbor_pah$pheno_meta %in% kclust8~ 8))
#                               
# 
# # go through all clusters and get mean frequency of each cell type
# cell_types <- colnames(immune_neighbor_pah[,3:14])
# hm_allclusters <- matrix(, nrow = length(unique(immune_neighbor_pah$pheno_meta)), ncol = length(cell_types))
# for(i in 1:nrow(hm_allclusters)) {
#   temp_mat <-immune_neighbor_pah[immune_neighbor_pah[,"pheno_meta"] == i, cell_types]
#   hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
# }
# hm_allclusters[is.na(hm_allclusters)] <- 0
# 
# # add names to rows and cols
# rownames(hm_allclusters) <- paste("cluster", 1: nrow(hm_allclusters), sep = "")
# colnames(hm_allclusters) <- cell_types  
# hm_allclusters
# 
# # plot heatmap of all clusters
# 
# dev.off()
# heatmap.2(hm_allclusters, 
#           Colv = T, Rowv = T,
#           hclustfun = hclust,
#           scale = "column",
#           dendrogram = c("both"),
#           trace = "none",
#           col = diverging_hcl(100, palette = 'Blue-Red 3'),
#           #col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
#           sepcolor="grey35",
#           colsep=0:ncol(hm_allclusters),
#           rowsep=0:nrow(hm_allclusters),
#           sepwidth=c(0.01,0.01),
#           density.info = 'none',
#           key = TRUE,
#           key.title = '',
#           cexRow = 0.6, cexCol = 0.9)


