
####DC analysis###

  
# Visualize the counts for each cluster by condition
options(stringsAsFactors = F)
rm(list =ls())
setwd("~/Desktop/FinalData/PAHanalyses/single_cell_data")
dir()
data_norm_cells<-read.csv("ImmuneMetaclust_K=30.csv")
head(data_norm_cells)

PAH<-data_norm_cells[which(data_norm_cells$Group=='PAH'),]
head(PAH)
data_norm_cells<-PAH

data_immune<-droplevels(data_norm_cells[data_norm_cells$cell_lineage == 'DC',])
dim(data_immune)
head(data_immune)

ff_new <- flowFrame(exprs = data.matrix(data_immune[,-c(1,2,3,39,40,41,42,43,44,45,46,47,48,49,50,51,52)]), desc = list(FIL = 1)) #exclude point name column 46

#clusterChannels_1<-c('CD14','CD11c','CD141','CD209') # also works
clusterChannels_1<-c('CD14','CD163','CD11c','CD11b','CD141','CD209','HLA.DR')
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

#write.csv(data_fs_clusters,"~/Desktop/FinalData/PAHanalyses/DCClusters.csv", row.names = FALSE)

# check the amount of cells in each immune cluster
table(data_fs_clusters$fs_clusters)

# go through all clusters and calculate mean for every channel
hm_allclusters <- matrix(, nrow = length(unique(fs_clusters)), ncol = length(clusterChannels_1))
for(i in 1:length(unique(fs_clusters))) {
  clust = fs_clusters[i]
  temp_mat <- data_fs_clusters[data_fs_clusters[,"fs_clusters"] == clust, clusterChannels_1]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

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
auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 90)
max(auto_meta)

# do a manual metaclustering
chosen_k=8
set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = chosen_k)
meta_results2 <- out_meta[fs_clusters]


##..Visualize FlowSOM metaclusters output on heatmap..##

# make combined expression matrix
data_meta_clusters <- cbind(data_immune, meta_results2)

# check the amount of cells in each cluster
table(data_meta_clusters$meta_results2)

# go through all clusters and calculate mean for every channel
hm_metaclusters <- matrix(, nrow = chosen_k, ncol = length(clusterChannels_1))
for(i in 1:chosen_k) {
  temp_mat <- data_meta_clusters[data_meta_clusters[,"meta_results2"] == i, clusterChannels_1]
  hm_metaclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

#hm_allclusters[is.na(hm_allclusters)] <- 0
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
          cexRow = 0.8, cexCol = 0.8, 
          breaks=seq(0, 1, length.out=257))

##..Label metaclusters to annotate data..##

CD209_clust<-c(1,5)
cDC1_clust<-c(7)
cDC2_clust<-c(3,4)
MDSC_clust<-c(2,6,8)

#CD209_clust<-c(10,6,13,11,9)  --with option 1
#cDC1_clust<-c(8,3,7,1) --with option 1
#cDC2_clust<-c(2,5,4) --with option 1
#MDSC_clust<-c(12) --with option 1


cell_type<-meta_results2

cell_type<-replace(cell_type,cell_type %in% CD209_clust,"moDC")
cell_type<-replace(cell_type,cell_type %in% cDC1_clust,"cDC1")
cell_type<-replace(cell_type,cell_type %in% cDC2_clust,"cDC2")
cell_type<-replace(cell_type,cell_type %in% MDSC_clust,"pDC/MDSC")


head(data_immune)


##..Add cell type to each event..##
data_meta_clusters$cell_lineage<-cell_type
head(data_meta_clusters)


# make combined expression matrix
#data_meta_clusters <- cbind(data_immune, meta_results2)
#head(data_meta_clusters)

write.csv(data_meta_clusters,"~/Desktop/FinalData/PAHanalyses/DCMetaclust_K=20_corrected.csv", row.names = FALSE)


data_meta_clusters<-read.csv("~/Desktop/FinalData/PAHanalyses/DCMetaclust_K=20_corrected.csv")
head(data_meta_clusters)
dim(data_meta_clusters)

#check the amount of cells in each cluster and point
table(data_meta_clusters$meta_results2)
table(data_meta_clusters$Point_num)

cell_counts_clusters<-table(data_meta_clusters$Point_num)
head(cell_counts_clusters)
write.csv(cell_counts_clusters,"~/Desktop/FinalData/PAHanalyses/DC_CountsPerPoint_corrected.csv", row.names = FALSE)


# Visualize the counts for each cluster by condition

#Group
cluster_group_summary<-as.data.frame(table(data_meta_clusters$Group, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary)<-c('Group','cell_type','Count')
Tot_Immune<-as.data.frame(table(data_meta_clusters$Group))$Freq
head(cluster_group_summary)
p <- cbind(cluster_group_summary, Tot_Immune)
p$Cluster_Freq<-(p$Count/p$Tot_Immune)*100
p
write.csv(p,"~/Desktop/FinalData/PAHanalyses/DC_CellFreqPerGroup_corrected.csv", row.names = FALSE)


#Patient
cluster_group_summary1<-as.data.frame(table(data_meta_clusters$PID, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary1)<-c('PID','cell_type','Count')
Tot_Immune<-as.data.frame(table(data_meta_clusters$PID))$Freq
head(cluster_group_summary1)
f <- cbind(cluster_group_summary1, Tot_Immune)
f$Cluster_Freq<-(f$Count/f$Tot_Immune)*100
f
write.csv(f,"~/Desktop/FinalData/PAHanalyses/DC_K=20CellFreqPerPID_corrected.csv", row.names = FALSE)


#Subgroup/Point
cluster_group_summary2<-as.data.frame(table(data_meta_clusters$Subgroup, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary2)<-c('Subgroup','cell_type','Count')
Tot_Immune<-as.data.frame(table(data_meta_clusters$Subgroup))$Freq
head(cluster_group_summary2)
s <- cbind(cluster_group_summary2, Tot_Immune)
s$Cluster_Freq<-(s$Count/s$Tot_Immune)*100
s
write.csv(s,"~/Desktop/FinalData/PAHanalyses/DC_CellFreqPerSubgroup_Point_corrected.csv", row.names = FALSE)

#Point
cluster_group_summary3<-as.data.frame(table(data_meta_clusters$Point_num, data_meta_clusters$cell_lineage))
colnames(cluster_group_summary3)<-c('Point_num','cell_type','Count')
Tot_Immune<-as.data.frame(table(data_meta_clusters$Point_num))$Freq
head(cluster_group_summary3)
q <- cbind(cluster_group_summary3, Tot_Immune)
q$Cluster_Freq<-(q$Count/q$Tot_Immune)*100
q
write.csv(q,"~/Desktop/FinalData/PAHanalyses/DC_CellFreqPerPoint_corrected.csv", row.names = FALSE)



