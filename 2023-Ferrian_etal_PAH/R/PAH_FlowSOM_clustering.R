# PAH_MIBIrunFlowSOM.R
# Author: Erin McCaffrey - some modifications by Selena 
# Date created: 02/21/2020
# Overview: Script imports concatenated data set of asinh transformed sc data. Runs FlowSOM on the transformed and
# normalized data.


###########################################
##..Install packages and open libraries..##
###########################################

# install.packages("BiocManager")
# BiocManager::install("BiocUpgrade")  # upgrading Bioconductor
# BiocManager::install("flowCore")     # for handling of fcs files in R
# BiocManager::install("ggplot2")      # for advanced data plotting
# BiocManager::install("gplots")       # for heatmaps
# BiocManager::install("RColorBrewer") # additional color schemes
# BiocManager::install("reshape2")     # reorganizing data
# BiocManager::install("FlowSOM")      # FlowSOM clustering
# BiocManager::install("plyr")         # for data handling
# install.packages("viridis")          # for badass color palettes
# install.packages("naniar")           # for visualize missing values

##..Open necessary libraries..##

library(flowCore)           
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(FlowSOM)
library(plyr)
library(dplyr)
library(viridis)


################################################################################
##..Data Preparation of data , CS-normalized, scaled, and asinh trans..##
################################################################################

##..Figure1A..##

options(stringsAsFactors = F)
rm(list =ls())
setwd("~/Desktop/FinalData/PAHanalyses")
library(ggplot2)
library(ggmosaic)
dir()
mydata <- read.table("MosaicPlot.txt", header=T, sep="\t")
dim(mydata)
head(mydata)
sapply(mydata,function(x) sum(is.na(x)))
classes <- sapply(mydata, class)
classes
str(mydata)
data=table(mydata$ORDER.1, mydata$Score, mydata$PID)
mosaicplot(data,main="PID",xlab="ORDER.1",ylab="Score",col=c(3,4,6))
data


##..Read in the concatenated expression matrix..##
options(stringsAsFactors = F)
rm(list =ls())
setwd("~/Desktop/FinalData/PAHanalyses/single_cell_data")
dir()
data_CSscaled<-read.csv("combined_dataT.csv") #cell size normalized, linearly scaled
head(data_CSscaled)
dim(data_CSscaled)

classes <- sapply(data_CSscaled, class)
classes


##..Percent normalization and scale data 0-99th percentile..##
v <- 1:1000
v
quantile_value <- 0.999
quantile(v, quantile_value)
data_trans<- data_CSscaled


# calculating percentile for 1 vector
percentile.vector <- apply(data_trans[,c(4:21,23:32,34:40)], 2, function(x) quantile(x, quantile_value, names = F))
percentile.vector
data_norm<-data_trans


data_norm[,c(4:21,23:32,34:40)] <- data.frame(t(t(data_trans[,c(4:21,23:32,34:40)]) / as.numeric(percentile.vector)))

#max normalize Foxp3 due to low positive cell count affecting 99.9th percentile
foxp3_max<-max(data_trans$FoxP3)
PD1_max<-max(data_trans$PD1)
data_norm[,22] <- data_norm[,22]/foxp3_max
data_norm[,33] <- data_norm[,33]/PD1_max
head(data_norm[,c(22,33)])

# filter out cells based on H3 signal to exlcude artifacts
hist(data_norm$H3, breaks = 20)
max(data_norm$H3)
min(data_norm$H3)
sum(data_norm$H3[data_norm$H3<min(data_norm$H3)])
sum(data_norm$H3)
sum(data_norm$Na)
sum(data_norm$H3[data_norm$H3<0.1])
sum(data_norm$H3[data_norm$H3>0.1])

H3_cutoff <- 0.1
data_norm_cells <- droplevels(data_norm[data_norm$H3 > H3_cutoff,])
sum(data_norm_cells$H3)
dim(data_norm_cells)
hist(data_norm_cells$H3)
summary(data_norm_cells)


write.csv(data_norm_cells,"~/Desktop/FinalData/PAHanalyses/FlowSom_Data_Norm.csv", row.names = FALSE)


#######################################################################

        ###..FlowSOM Cluster 1: Separate major lineages..### 

#######################################################################
#Subset PAH patients
#PAH<-data_norm_cells[which(data_norm_cells$Group=='PAH'),]
#head(PAH)
#data_norm_cells<-PAH


dev.off()
dim(data_norm_cells)
head(data_norm_cells)
ff_new <- flowFrame(exprs = data.matrix(data_norm_cells[,-c(47,48,49,50,51,52,53)]), desc = list(FIL = 1)) #exclude point name column 46

clusterChannels_1<-c('SMA','PanCK','CD31','CD14','Vimentin','CD45')

##..Run FlowSOM random seed for reproducibility..##

set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusterChannels_1)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
fs_clusters <- out_fSOM$map$mapping[,1]

out_fSOM <- UpdateNodeSize(out_fSOM, reset = TRUE)
FlowSOM::PlotStars(out_fSOM, view = "grid", markers = clusterChannels_1)
out_fSOM <- UpdateNodeSize(out_fSOM)
FlowSOM::PlotStars(out_fSOM, view = "MST", markers = clusterChannels_1)
FlowSOM::PlotStars(out_fSOM, view="tSNE",markers = clusterChannels_1)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
data_fs_clusters <- cbind(data_norm_cells, fs_clusters)
head(data_fs_clusters)

write.csv(data_fs_clusters,"~/Desktop/FinalData/PAHanalyses/LineageClusters.csv", row.names = FALSE)


# check the amount of cells in each cluster
table(data_fs_clusters$fs_clusters)

# go through all clusters and calculate mean for every channel
hm_allclusters <- matrix(, nrow = length(unique(fs_clusters)), ncol = length(clusterChannels_1))
for(i in 1:length(unique(fs_clusters))) {
  clust = fs_clusters[i]
  temp_mat <- data_fs_clusters[data_fs_clusters[,"fs_clusters"] == clust, clusterChannels_1]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- paste("cluster", unique(fs_clusters), sep = "")
colnames(hm_allclusters) <- clusterChannels_1  
hm_allclusters

heatmap.2(hm_allclusters, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))

#write.csv(hm_allclusters,"~/Desktop/FinalData/FinalDataDenoise_adjusted/segmentation_output/single_cell_output_expansion_3/allclusters.csv", row.names = FALSE)


##..Meta-cluster..##

# try the suggested automatic metaclustering method for a hint for k
auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 20)
max(auto_meta)

# do a manual metaclustering
chosen_k=5
set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = chosen_k)
meta_results <- out_meta[fs_clusters]

##..Visualize FlowSOM metaclusters output on heatmap..##

# make combined expression matrix
data_meta_clusters <- cbind(data_norm_cells, meta_results)

# check the amount of cells in each cluster
table(data_meta_clusters$meta_results)

# go through all clusters and calculate mean for every channel
hm_metaclusters <- matrix(, nrow = chosen_k, ncol = length(clusterChannels_1))
for(i in 1:chosen_k) {
  temp_mat <- data_meta_clusters[data_meta_clusters[,"meta_results"] == i, clusterChannels_1]
  hm_metaclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
rownames(hm_metaclusters) <- paste("cluster", 1:chosen_k, sep = "")
colnames(hm_metaclusters) <- clusterChannels_1
hm_metaclusters

dev.off() #if the plot shows errors

# make a metacluster heatmap

#heatmap.2(hm_metaclusters,
          #scale = "none",
          #Colv = T, Rowv = T,
          #hclustfun = hclust,
          #dendrogram = c("both","row","column","none"),
          #trace = "none",
          #col = colorRampPalette(brewer.pal(9,"Blues"))(100),
          #density.info = 'none',
          #key.title = '',
          #lhei = c(1,7),
          #cexRow = 0.7, cexCol = 0.7, 
          #breaks=seq(0, 1, length.out=101))

heatmap.2(hm_metaclusters, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = cividis(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))

#write.csv(hm_metaclusters,"~/Desktop/FinalData/FinalDataDenoise_adjusted/segmentation_output/single_cell_output_expansion_3/hm_metaclusters.csv", row.names = TRUE)

##..Label metaclusters to annotate data..##

fibro_clust<-4
endo_clust<-5
imm_clust<-3
mesen_clust<-2
epi_clust<-1
cell_type<-meta_results

cell_type<-replace(cell_type,cell_type %in% imm_clust,"immune")
cell_type<-replace(cell_type,cell_type %in% epi_clust,"epithelial")
cell_type<-replace(cell_type,cell_type %in% endo_clust,"endothelial")
cell_type<-replace(cell_type,cell_type %in% fibro_clust,"fibroblast")
cell_type<-replace(cell_type,cell_type %in% mesen_clust,"mesenchymal")
head(data_norm_cells)

##..Add cell type to each event..##
data_norm_cells$cell_lineage<-cell_type

# make combined expression matrix
data_meta_clusters <- cbind(data_norm_cells, meta_results)
head(data_meta_clusters)

write.csv(data_meta_clusters,"~/Desktop/FinalData/PAHanalyses/LineageMetaclust_k=5.csv", row.names = FALSE)

# check the amount of cells in each cluster and point
table(data_meta_clusters$meta_results)
table(data_meta_clusters$point)

cell_counts_clusters<-table(data_meta_clusters$point)
head(cell_counts_clusters)
write.csv(cell_counts_clusters,"~/Desktop/FinalData/PAHanalyses/LineageMetaclust_K=5CellCountsPerPoint.csv", row.names = FALSE)

#Group Healthy Controls/PAH
cluster_group_summary<-as.data.frame(table(data_meta_clusters$Group, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary)<-c('Group','cell_type','Count')
totals<-as.data.frame(table(data_meta_clusters$Group))$Freq
head(cluster_group_summary)
p <- cbind(cluster_group_summary, totals)
p$Cluster_Freq<-(p$Count/p$totals)*100
Cluster_Freq<-p$Cluster_Freq
p2 <- cbind(p, Cluster_Freq)
p2
write.csv(p2,"~/Desktop/FinalData/PAHanalyses/LineageMetaclust_k=5CellFreqPerGroup.csv", row.names = FALSE)

#Subgroup HPAH/IPAH
cluster_group_summary1<-as.data.frame(table(data_meta_clusters$Subgroup, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary1)<-c('Subgroup','cell_type','Count')
totals1<-as.data.frame(table(data_meta_clusters$Subgroup))$Freq
head(cluster_group_summary1)
s <- cbind(cluster_group_summary1, totals1)
s$Cluster_Freq<-(s$Count/s$totals1)*100
Cluster_Freq<-s$Cluster_Freq
s2 <- cbind(s, Cluster_Freq)
s2
write.csv(s2,"~/Desktop/FinalData/PAHanalyses/LineageMetaclust_k=5CellFreqPerSubgroup.csv", row.names = FALSE)

#Patient
cluster_group_summary2<-as.data.frame(table(data_meta_clusters$PID, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary2)<-c('PID','cell_type','Count')
totals2<-as.data.frame(table(data_meta_clusters$PID))$Freq
head(cluster_group_summary2)
f <- cbind(cluster_group_summary2, totals2)
f$Cluster_Freq<-(f$Count/f$totals2)*100
Cluster_Freq<-f$Cluster_Freq
f2 <- cbind(f, Cluster_Freq)
f2
write.csv(f2,"~/Desktop/FinalData/PAHanalyses/LineageMetaclust_K=5CellFreqPerPID.csv", row.names = FALSE)


#######################################################################
        ###..FlowSOM Cluster 1: Immune Phenotype..###
#######################################################################

# filter out cells based on CD45 signal to exclude artifacts
#hist(data_norm_cells$CD45, breaks = 10)
#max(data_norm_cells$CD45)
#min(data_norm_cells$CD45)
#sum(data_norm_cells$CD45[data_norm_cells$CD45>min(data_norm_cells$CD45)])
#sum(data_norm_cells$CD45)
#sum(data_norm_cells$CD45[data_norm_cells$CD45<0.025])

#CD45_cutoff <- 0
#data_immune <- droplevels(data_immune[data_immune$CD45 > CD45_cutoff,])


data_immune<-droplevels(data_norm_cells[data_norm_cells$cell_lineage == 'immune',])
dim(data_immune)
head(data_immune)
#sum(data_immune$CD45[data_immune$CD45>0]) #how many CD45+ cells have the immune cluster

ff_new <- flowFrame(exprs = data.matrix(data_immune[,-c(46,47,48,49,50,51,52,53,54)]), desc = list(FIL = 1)) #exclude point name column 46

clusterChannels_1<-c('CD14','CD15','CD4','CD11c','CD68','CD11b','CD8','CD20','CD3e','CD16','CD209','CD163')
#clusterChannels_1<-c('CD14','CD15','CD4','CD11c','CD68','CD11b','CD8','CD20','CD3e','CD16')
#clusterChannels_1<-c('CD14','CD15','CD11c','CD68','CD11b','CD20','CD3e','CD16')
##..Run FlowSOM random seed for reproducibility..##

set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusterChannels_1)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
fs_clusters <- out_fSOM$map$mapping[,1]

out_fSOM <- UpdateNodeSize(out_fSOM, reset = TRUE)
FlowSOM::PlotStars(out_fSOM, view = "grid", markers = clusterChannels_1)
out_fSOM <- UpdateNodeSize(out_fSOM)
FlowSOM::PlotStars(out_fSOM, view = "MST", markers = clusterChannels_1)
FlowSOM::PlotStars(out_fSOM, view="tSNE",markers = clusterChannels_1)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
data_fs_clusters <- cbind(data_immune, fs_clusters)

write.csv(data_fs_clusters,"~/Desktop/FinalData/PAHanalyses/ImmuneClusters.csv", row.names = FALSE)

# check the amount of cells in each immune cluster
table(data_fs_clusters$fs_clusters)

# go through all clusters and calculate mean for every channel
hm_allclusters <- matrix(, nrow = length(unique(fs_clusters)), ncol = length(clusterChannels_1))
for(i in 1:length(unique(fs_clusters))) {
  clust = fs_clusters[i]
  temp_mat <- data_fs_clusters[data_fs_clusters[,"fs_clusters"] == clust, clusterChannels_1]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- paste("cluster", unique(fs_clusters), sep = "")
colnames(hm_allclusters) <- clusterChannels_1  
hm_allclusters

dev.off()
heatmap.2(hm_allclusters, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))

# run if you get a plot margin error


##..Meta-cluster..##

# try the suggested automatic metaclustering method for a hint for k
#auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 90)
#max(auto_meta)

# do a manual metaclustering
chosen_k=30
set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = chosen_k)
meta_results <- out_meta[fs_clusters]


##..Visualize FlowSOM metaclusters output on heatmap..##

# make combined expression matrix
data_meta_clusters <- cbind(data_immune, meta_results)

# check the amount of cells in each cluster
table(data_meta_clusters$meta_results)

# go through all clusters and calculate mean for every channel
hm_metaclusters <- matrix(, nrow = chosen_k, ncol = length(clusterChannels_1))
for(i in 1:chosen_k) {
  temp_mat <- data_meta_clusters[data_meta_clusters[,"meta_results"] == i, clusterChannels_1]
  hm_metaclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
rownames(hm_metaclusters) <- paste("cluster", 1:chosen_k, sep = "")
colnames(hm_metaclusters) <- clusterChannels_1
hm_metaclusters

# make a metacluster heatmap

dev.off()

heatmap.2(hm_metaclusters,
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(brewer.pal(7, "RdPu"))(100),
          #col = colorRampPalette(brewer.pal(7, "Set1"))(100),
          #col = magma(256),
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))

##..Label metaclusters to annotate data..##

CD3_clust<-9
CD3CD8_clust<-c(5,2,1,7,3)
CD3CD4_clust<-c(8,24,15,29)

DC_clust<-c(21,22,20,11,13,25,27,28,22,26,21,30)

CD68_clust<-c(23,18,16,17)

CD14_clust<-c(19)
CD16_clust<-c(14)

CD15_clust<-c(10,12)
CD20_clust<-c(4,6)

cell_type<-meta_results

cell_type<-replace(cell_type,cell_type %in% CD3_clust,"Th")

cell_type<-replace(cell_type,cell_type %in% CD3CD8_clust,"Tc")
cell_type<-replace(cell_type,cell_type %in% CD3CD4_clust,"Th")

cell_type<-replace(cell_type,cell_type %in% DC_clust,"DC")

cell_type<-replace(cell_type,cell_type %in% CD68_clust,"Macro")
cell_type<-replace(cell_type,cell_type %in% CD14_clust,"Mono")
cell_type<-replace(cell_type,cell_type %in% CD16_clust,"NK")

cell_type<-replace(cell_type,cell_type %in% CD15_clust,"Neutro")
cell_type<-replace(cell_type,cell_type %in% CD20_clust,"Bcell")


head(data_immune)

##..Add cell type to each event..##
data_immune$cell_lineage<-cell_type
head(data_immune)


# make combined expression matrix
data_meta_clusters <- cbind(data_immune, meta_results)
head(data_meta_clusters)

write.csv(data_meta_clusters,"~/Desktop/FinalData/PAHanalyses/ImmuneMetaclust_K=30.csv", row.names = FALSE)


#check the amount of cells in each cluster and point
table(data_meta_clusters$meta_results)
table(data_meta_clusters$point)

cell_counts_clusters<-table(data_meta_clusters$point)
head(cell_counts_clusters)
write.csv(cell_counts_clusters,"~/Desktop/FinalData/PAHanalyses/ImmuneMetaclust_K=30CellCountsPerPoint.csv", row.names = FALSE)


# Visualize the counts for each cluster by condition
#setwd("~/Desktop/FinalData/PAHanalyses")
#data_meta_clusters<-read.csv("ImmuneMetaclust_K=30_FINAL.csv")

#Group
cluster_group_summary<-as.data.frame(table(data_meta_clusters$Group, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary)<-c('Group','cell_type','Count')
totals<-as.data.frame(table(data_meta_clusters$Group))$Freq
head(cluster_group_summary)
p <- cbind(cluster_group_summary, totals)
p$Cluster_Freq<-(p$Count/p$totals)*100
Cluster_Freq<-p$Cluster_Freq
p2 <- cbind(p, Cluster_Freq)
p2
write.csv(p2,"~/Desktop/FinalData/PAHanalyses/ImmuneMetaclust_K=30CellFreqPerGroup.csv", row.names = FALSE)

#Patient
cluster_group_summary2<-as.data.frame(table(data_meta_clusters$PID_name, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary2)<-c('PID_name','cell_type','Count')
totals2<-as.data.frame(table(data_meta_clusters$PID_name))$Freq
head(cluster_group_summary2)
f <- cbind(cluster_group_summary2, totals2)
f$Cluster_Freq<-(f$Count/f$totals2)*100
Cluster_Freq<-f$Cluster_Freq
f2 <- cbind(f, Cluster_Freq)
f2
write.csv(f2,"~/Desktop/FinalData/PAHanalyses/ImmuneMetaclust_K=30CellFreqPerPID.csv", row.names = FALSE)
