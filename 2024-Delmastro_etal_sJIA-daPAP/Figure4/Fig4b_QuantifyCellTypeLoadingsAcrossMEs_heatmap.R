# Fig4b_QuantifyCellTypeLoadingsAcrossMEs.R
# Author: Erin McCaffrey
# Date created: 191211
# Date updated: 240720
# Overview: This script reads in the csv for the ME loadings for each cell and visualizes the mean
# probability for each cell type across MEs as a heatmap

library(dplyr)
library(gplots)
require(pals)
library(colorspace)

##...Load in data..##

setwd("/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points")
topic_wts<-read.csv("Neighborhood_Analysis/Spatial_LDA/model_phenotype_75px/topic_weights.csv")
data <- read.csv("celltable_05032024.csv")

##..Merge topic loadings with cell table..##

ME_data <- merge(topic_wts, data[,colnames(data) %in% c('point','label','name')], by = c("point", "label"))

##...Define MEs and cell_types...##

phenos<-as.vector(unique(ME_data$name))
MEs<-as.vector(colnames(ME_data[,3:11]))

##...For each topic and cell type, produce heatmap with the mean loading...##

hm_alltopics <- matrix(, nrow = length(phenos) , ncol=length(MEs))
for(i in 1:length(phenos)) {
  temp_mat <- ME_data[ME_data$name==phenos[i], MEs]
  hm_alltopics[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

##..Redefine cell order..##

cell_order <- c('Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg","CD11c+_mDC","CD14+_Mono",
                            "CD209+_Mac","iNOS+_Mac","M2_Mac","Giant_cell", "Eosinophil","Mast_cell","Neutrophil","CD16+_ImmuneOther",
                            "CD57+_ImmuneOther","Endothelial", "Fibroblast","iNOS+_Pneumocyte","Lung_Epithelium","Mesenchymal")

rownames(hm_alltopics) <- phenos
colnames(hm_alltopics) <- MEs
hm_alltopics <- hm_alltopics[rev(cell_order),]

##..Visualize..##

pdf(file = "Sub-panels/Figure4/Fig4b_ME_heatmap.pdf")
heatmap.2(hm_alltopics, 
          Colv = F, Rowv = F,
          #hclustfun = hclust,
          scale = "row",
          dendrogram = "none",
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          sepcolor="#999999",
          colsep=0:ncol(hm_alltopics),
          rowsep=0:nrow(hm_alltopics),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '', 
          cexRow = 1, cexCol = 1, margins = c(8,14), 
          breaks=seq(-2, 2, length.out=101))
dev.off()

