# PAH_MIBIheatmap.R
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


##..Define heatmap markers and lineages..##

lineage_markers<-c('CD11b','CD11c','CD14','CD15','CD16','CD163','CD20','CD209','CD31','CD4',
                   'CD3e','CD45','CD68','CD8','HLA.DR','PanCK','MPO','SMA','Vimentin') 

cell_lineages<-unique(data$cell_lineage)

hm_allclusters <- matrix(, nrow = length(cell_lineages), ncol = length(lineage_markers))
for(i in 1:length(cell_lineages)) {
  temp_mat <- data[data$cell_lineage==cell_lineages[i], lineage_markers]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_allclusters) <- paste(cell_lineages)
colnames(hm_allclusters) <- lineage_markers
hm_allclusters

##..Get vector of color codes for annotation column..##

row_order<-c('DC','Mono','Macro','Th','Tc','NK','Bcell','endothelial','epithelial',
             'mesenchymal','fibroblast','Neutro')
colorkey<-read.csv('colorkey.csv')
colorkey$cell_lineage<-factor(colorkey$cell_lineage, levels = cell_lineages)
colorkey<-colorkey[order(colorkey$cell_lineage),]
color<-as.vector(colorkey$codes)

##..Plot heatmap..##

# pdf("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Figures/heatmaps/all_cells.pdf",width=15, height=11)
heatmap.2(hm_allclusters,
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "both",
          trace = "none",
          # col = magma(256),
          col = as.vector((ocean.thermal(100))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          cexRow=0.4,
          cexCol=0.4,
          RowSideColors = color,
          breaks=seq(0, 1, length.out=101))
# dev.off()


##..Make separate healthy and PAH and separate immune and non-immune matrices for functional markers..##

imm<-c('Th','Tc','Bcell','NK','Macro','Mono','Neutro','DC')
nonimm<-c('epithelial','endothelial','mesenchymal','fibroblast')

hlt_imm<-droplevels(data[data$Group=='Healthy Controls'& data$cell_lineage %in% imm,])
hlt_nonimm<-droplevels(data[data$Group=='Healthy Controls'& !data$cell_lineage %in% imm,])

pah_imm<-droplevels(data[data$Group=='PAH' & data$cell_lineage %in% imm,])
pah_nonimm<-droplevels(data[data$Group=='PAH' & !data$cell_lineage %in% imm,])


functional_imm_markers<-c('CD45RO','GranzB','CD163','HLA.Class.I','HLA.DR','IDO1','IFNg','iNOS','KI67',
                      'MPO','PDL1b','SAMHD1','TIM.3')
functional_nonimm_markers<-c('bCatenin','SAMHD1','HLA.Class.I','HLA.DR','CD141','iNOS','KI67','IFNg',
                             'NaK.ATPase')


hm_hlt_imm <- matrix(, nrow = length(imm), ncol = length(functional_imm_markers))
hm_pah_imm <- hm_hlt_imm
hm_imm <- hm_hlt_imm

for(i in 1:length(imm)) {
  #hlt
  temp_mat_hlt <- hlt_imm[hlt_imm$cell_lineage==imm[i], functional_imm_markers]
  hm_hlt_imm[i,] <- as.matrix(apply(temp_mat_hlt, 2, function (x) mean(x, na.rm = T)))
  #pah
  temp_mat_pah <- pah_imm[pah_imm$cell_lineage==imm[i], functional_imm_markers]
  hm_pah_imm[i,] <- as.matrix(apply(temp_mat_pah, 2, function (x) mean(x, na.rm = T)))
  #merged
  temp_mat <- data[data$cell_lineage==imm[i], functional_imm_markers]
  hm_imm[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

hm_hlt_nonimm <- matrix(, nrow = length(nonimm), ncol = length(functional_nonimm_markers))
hm_pah_nonimm <- hm_hlt_nonimm
hm_nonimm <- hm_hlt_nonimm

for(i in 1:length(nonimm)) {
  #hlt
  temp_mat_hlt <- hlt_nonimm[hlt_nonimm$cell_lineage==nonimm[i], functional_nonimm_markers]
  hm_hlt_nonimm[i,] <- as.matrix(apply(temp_mat_hlt, 2, function (x) mean(x, na.rm = T)))
  #pah
  temp_mat_pah <- pah_nonimm[pah_nonimm$cell_lineage==nonimm[i], functional_nonimm_markers]
  hm_pah_nonimm[i,] <- as.matrix(apply(temp_mat_pah, 2, function (x) mean(x, na.rm = T)))
  #merged
  temp_mat <- data[data$cell_lineage==nonimm[i], functional_nonimm_markers]
  hm_nonimm[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

#immune
rownames(hm_pah_imm) <- paste(imm)
colnames(hm_pah_imm) <- functional_imm_markers

rownames(hm_hlt_imm) <- paste(imm)
colnames(hm_hlt_imm) <- functional_imm_markers

rownames(hm_imm) <-paste(imm)
colnames(hm_imm) <-functional_imm_markers

#non-immune
rownames(hm_hlt_nonimm) <- paste(nonimm)
colnames(hm_hlt_nonimm) <- functional_nonimm_markers

rownames(hm_pah_nonimm) <- rownames(hm_hlt_nonimm)
colnames(hm_pah_nonimm) <- colnames(hm_hlt_nonimm)

rownames(hm_nonimm) <- rownames(hm_hlt_nonimm)
colnames(hm_nonimm) <- colnames(hm_hlt_nonimm)

# get color order for imm plots

colorkey_imm<-droplevels(colorkey[colorkey$cell_lineage %in% imm, ])
colorkey_imm$cell_lineage<-factor(colorkey_imm$cell_lineage, levels = imm)
colorkey_imm<-colorkey_imm[order(colorkey_imm$cell_lineage),]
color<-as.vector(colorkey_imm$codes)

heatmap.2(hm_imm,
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
          colsep=0:ncol(hm_hlt_imm),
          rowsep=0:nrow(hm_hlt_imm),
          cexRow=1,
          cexCol=1,
          RowSideColors = color)

# non-immune cells

colorkey_nonimm<-droplevels(colorkey[colorkey$cell_lineage %in% nonimm, ])
colorkey_nonimm$cell_lineage<-factor(colorkey_nonimm$cell_lineage, levels = nonimm)
colorkey_nonimm<-colorkey_nonimm[order(colorkey_nonimm$cell_lineage),]
color<-as.vector(colorkey_nonimm$codes)

heatmap.2(hm_nonimm,
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
          colsep=0:ncol(hm_pah_nonimm),
          rowsep=0:nrow(hm_pah_nonimm),
          cexRow=1,
          cexCol=1,
          RowSideColors = color)

##..Make separate healthy and PAH and separate immune and non-immune matrices for functional markers..##

data_neighborhood_markers<-merge(data_neighborhood, data, by=c("Point_num","label",'PID'))

functional_markers<-c('CD45RO','GranzB','CD163','HLA.DR','IDO1','IFNg','iNOS','KI67',
                      'PDL1b','TIM.3','CD141')

clusters<-c(1,2,3,4,5,6,7,8,9,10)

hm_neighborhood <- matrix(, nrow = length(clusters), ncol = length(functional_markers))

for(i in 1:length(clusters)) {
  temp_mat <- data_neighborhood_markers[data_neighborhood_markers$cluster == clusters[i], functional_markers]
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


