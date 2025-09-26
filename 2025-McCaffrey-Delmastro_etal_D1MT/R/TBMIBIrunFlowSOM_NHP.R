# TBMIBIrunFlowSOM_val.R
# Author: Erin McCaffrey (with many sections adapted from Felix's R Demo)
# Date created: 190128
# Overview: Script imports concatenated data set of asinh transformed sc data. Runs FlowSOM on the transformed and
# normalized data.

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
##..Data Preparation of A granulomas, CS-normalized, scaled, and asinh trans..##
################################################################################

##..Read in the concatenated expression matrix..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data_CSscaled<-read.csv("allsamples_dataCS_annotated.csv") #cell size normalized and untransformed

##..Drop controls and the granulomas with poor image quality issues..##

drop_samples<-c(3,4,5,6,7,17,26)

data_gran<-data_CSscaled[!data_CSscaled$SampleID %in% drop_samples,] 
# linear transform the data by 100 (only the expression data and not sample ID, cell label, cell size, or tissue)

data_gran[,4:46]<- data_gran[,4:46]*100

# arcsinh transform data (only the expression data and not sample ID, cell label, cell size, or tissue)

asinh_scale <- 5
data_trans<-data_gran
data_trans[,4:46]<- asinh(data_gran[,4:46]/ asinh_scale)

##..Percent normalization and scale data 0-99th percentile..##

v <- 1:1000
v
quantile_value <- 0.999
quantile(v, quantile_value)

# calculating percentile for 1 vector
percentile.vector <- apply(data_trans[,4:46], 2, function(x) quantile(x, quantile_value, names = F))
percentile.vector
data_gran_norm<-data_trans
data_gran_norm[,4:46] <- data.frame(t(t(data_trans[,4:46]) / as.numeric(percentile.vector)))
data_gran_norm<-data_gran_norm[,-c(4,19,23,24,27)] #remove channels with no expression or staining issues (C, PD1,)

# write.csv(data_gran_norm, 'gran_cohort_CSnorm-scaled.csv',row.names = F)

################################################################
###..FlowSOM Cluster round 1: Immune, Endo, and Fibroblast..###
###############################################################

ff_new <- flowFrame(exprs = data.matrix(data_gran_norm[,-c(42,43)]), desc = list(FIL = 1)) #exclude tissue type column

clusterChannels_1<-c('CD45','SMA','CD31','VIM',"CD3","CD20","CD14",
                     "HLA.DR","MastChyTry","Calprotectin","Keratin")

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
data_fs_clusters <- cbind(data_gran_norm, fs_clusters)

# go through all clusters and calculate mean for every channel
hm_allclusters <- matrix(, nrow = length(unique(fs_clusters)), ncol = length(clusterChannels_1))
for(i in 1:length(unique(fs_clusters))) {
  temp_mat <- data_fs_clusters[data_fs_clusters[,"fs_clusters"] == i, clusterChannels_1]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- paste("cluster", 1:length(unique(fs_clusters)), sep = "")
colnames(hm_allclusters) <- clusterChannels_1  
hm_allclusters

# plot heatmap of all clusters

heatmap.2(hm_allclusters, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
          #col = colorRampPalette(brewer.pal(9,"Blues"))(100),
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))


##..Meta-cluster..##

# try the suggested automatic metaclustering method for a hint for k
auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 20)
max(auto_meta)

# do a manual metaclustering
chosen_k=10
set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = chosen_k)
meta_results <- out_meta[fs_clusters]

##..Visualize FlowSOM metaclusters output on heatmap..##

# make combined expression matrix
data_meta_clusters <- cbind(data_gran_norm, meta_results)

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

heatmap.2(hm_metaclusters,
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both"),
          trace = "none",
          col = viridis(256),
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          density.info = 'none',
          key.title = '',
          cexRow = 0.7, cexCol = 0.7, 
          breaks=seq(0, 1, length.out=257))

# Assign clusters as being endothelial, fibroblast, neutrophil, or immune..##

imm_clust<-c(1,3,6,7,10)
fibro_clust<-c(2)
endo_clust<-c(4)
epi_clust<-c(5,8)
mast_clust<-c(9)
cell_type<-meta_results

cell_type<-replace(cell_type,cell_type %in% imm_clust,"immune")
cell_type<-replace(cell_type,cell_type %in% mast_clust,"mast")
cell_type<-replace(cell_type,cell_type %in% endo_clust,"endothelial")
cell_type<-replace(cell_type,cell_type %in% fibro_clust,"fibroblast")
cell_type<-replace(cell_type,cell_type %in% epi_clust,"epithelial")

data_gran_norm$cell_type<-cell_type


###################################################
####..FlowSOM Cluster round 2: Immune Cells..#####
###################################################

data_immune<-data_gran_norm[data_gran_norm$cell_type == 'immune',]
ff_immune<- flowFrame(exprs = data.matrix(data_immune[,-c(42,43,44)]), desc = list(FIL = 1))

clusterChannels_2=c("CD3","CD20","CD4","CD8","Foxp3","CD14","CD16","CD11c","CD11b","CD206","Calprotectin",
                    "CD68","CD163",'HLA.DR')

##..Run FlowSOM random seed for reproducibility..##

set.seed(123)
out_fSOM_imm <- FlowSOM::ReadInput(ff_immune, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_imm <- FlowSOM::BuildSOM(out_fSOM_imm, colsToUse = clusterChannels_2)
out_fSOM_imm <- FlowSOM::BuildMST(out_fSOM_imm)
labels_imm <- out_fSOM_imm$map$mapping[,1]

out_fSOM_imm <- UpdateNodeSize(out_fSOM_imm, reset = TRUE)
FlowSOM::PlotStars(out_fSOM_imm, view = "grid", markers = clusterChannels_2)
out_fSOM_imm <- UpdateNodeSize(out_fSOM_imm)
FlowSOM::PlotStars(out_fSOM_imm, view = "MST", markers = clusterChannels_2)
FlowSOM::PlotStars(out_fSOM_imm, view="tSNE",markers = clusterChannels_2)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_imm<-out_fSOM_imm[["map"]][["mapping"]]
fs_clusters_imm<-fs_clusters_imm[,1]
data_fs_clusters_imm <- cbind(data_immune, fs_clusters_imm)

# go through all clusters and calculate mean for every channel
hm_allclusters_imm <- matrix(, nrow = length(unique(fs_clusters_imm)), ncol = length(clusterChannels_2))
for(i in 1:length(unique(fs_clusters_imm))) {
  temp_mat <- data_fs_clusters_imm[data_fs_clusters_imm[,"fs_clusters_imm"] == i, clusterChannels_2]
  hm_allclusters_imm[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_imm) <- paste("cluster", 1:length(unique(fs_clusters_imm)), sep = "")
colnames(hm_allclusters_imm) <- clusterChannels_2  
hm_allclusters_imm

hmap_imm<-heatmap.2(hm_allclusters_imm, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both"),
          trace = "none",
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          cexRow = 0.3, cexCol = 0.6, 
          breaks=seq(0, 1, length.out=257))

hmap_imm<-heatmap.2(hm_allclusters_imm[c(42,51,73,75,83),], 
                    scale = "none",
                    Colv = T, Rowv = T,
                    hclustfun = hclust,
                    dendrogram = c("both"),
                    trace = "none",
                    col = viridis(256),
                    density.info = 'none',
                    key.title = '',
                    cexRow = 0.3, cexCol = 0.6, 
                    breaks=seq(0, 1, length.out=257))

# Assign cell types
#lymph
B_clust<-c(61,71,72,81,82,91,92)
CD8_clust <-c(50,59,60,69,70,79,80,90)
CD4_clust <-c(64,67,74,76,77,78,84,85,86,87,88,94,95,96,97,98,99,100)
treg_clust<-c(89)
tother_clust<-c(49,57,68,93)
#gran
neut_clust<-c(1,2,3,4,5,11,12,13,14,15,21,22,31)
#mac/mono
mac_mono_clust<-c(6,7,8,9,10,16,17,18,19,20,23,24,25,26,27,28,29,30,
                  32,33,34,35,36,37,38,39,40,41,43,44,45,46,47,48,52,53,54,55,56,58,62,63,65,66)
#other
other_clust<-c(42,51,73,75,83)


# Create vector of immune phenotypes

imm_pheno<-fs_clusters_imm
imm_pheno<-replace(imm_pheno, imm_pheno %in% B_clust,"B_cell")
imm_pheno<-replace(imm_pheno, imm_pheno %in% mac_mono_clust,"Mac_Mono")
imm_pheno<-replace(imm_pheno, imm_pheno %in% neut_clust,"neutrophil")
imm_pheno<-replace(imm_pheno, imm_pheno %in% CD4_clust,"CD4_T")
imm_pheno<-replace(imm_pheno, imm_pheno %in% CD8_clust,"CD8_T")
imm_pheno<-replace(imm_pheno, imm_pheno %in% other_clust,"imm_other")
imm_pheno<-replace(imm_pheno, imm_pheno %in% treg_clust,"Treg")
imm_pheno<-replace(imm_pheno, imm_pheno %in% tother_clust,"T_other")
data_immune$cell_type <- imm_pheno
  
# Create single dataset

data_gran_norm[rownames(data_immune),]$cell_type<-imm_pheno

# Add lineage data and codes

myeloid<-c("Mac_Mono")
lymphocyte<-c("CD8_T","CD4_T","B_cell","Treg","T_other")
granulocyte<-c("mast","neutrophil")
imm_other<-c("imm_other")
nonimmune<-c("fibroblast","endothelial","epithelial")


data_gran_norm<-data_gran_norm %>% mutate(cell_lin=case_when(data_gran_norm$cell_type %in% myeloid ~ "myeloid",
                                                             data_gran_norm$cell_type %in% lymphocyte ~ "lymphocyte",
                                                             data_gran_norm$cell_type %in% granulocyte ~ "granulocyte",
                                                             data_gran_norm$cell_type %in% imm_other ~ "other",
                                                             data_gran_norm$cell_type %in% nonimmune ~ "nonimmune")) 

data_gran_norm<-data_gran_norm %>% mutate(lintype_num=case_when(data_gran_norm$cell_lin == "myeloid" ~ 1,
                                                                data_gran_norm$cell_lin == "lymphocyte" ~ 2,
                                                                data_gran_norm$cell_lin == "granulocyte" ~ 3,
                                                                data_gran_norm$cell_lin == "other" ~ 4,
                                                                data_gran_norm$cell_lin == "nonimmune" ~ 5))

# Annotate giant cells

data_gran_norm[data_gran_norm$SampleID ==1 & data_gran_norm$cellLabelInImage %in% c(980,1218),]$cell_type <-'giant_cell'
data_gran_norm[data_gran_norm$SampleID ==14 & data_gran_norm$cellLabelInImage %in% c(1218,594,943),]$cell_type <-'giant_cell'
data_gran_norm[data_gran_norm$SampleID ==29 & data_gran_norm$cellLabelInImage %in% c(616,1328),]$cell_type <-'giant_cell'
data_gran_norm[data_gran_norm$SampleID ==41 & data_gran_norm$cellLabelInImage %in% c(820,718,1180),]$cell_type <-'giant_cell'

# Export
write.csv(data_gran_norm, 'NHP_cohort_data_norm_annotated.csv', row.names = F)

