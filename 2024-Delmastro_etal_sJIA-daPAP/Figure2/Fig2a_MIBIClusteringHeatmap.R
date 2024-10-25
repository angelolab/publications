# Fig2a_MIBIClusteringHeatmap.R
# Author: Erin McCaffrey 
# Date updated: 07/15/2024

require(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(viridis)
library(colorspace)
library(pals)
library(pheatmap)

##..Import data..##
setwd("/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points")
data<-read.csv("celltable_05032024.csv")

##..Define markers and lineages..##

heatmap_channels<-c('Tryptase','SMA','CD31','VIM','CD8','CD57','CD20','Foxp3','CD4','CD3','CD45','HLA.DR','PanCK','EPOX','Calprotectin',
           'iNOS','CD11c','CD68','CD14','CD209','CD206','CD163','CD16')
row_order<-c('Bcell','CD4+_Tcell','CD8+_Tcell','CD57+_CD8+_Tcell','Treg','CD11c+_mDC','CD14+_Mono','CD209+_Mac','iNOS+_Mac','M2_Mac',
             'Giant_cell','Eosinophil', 'Mast_cell', 'Neutrophil', 'CD16+_ImmuneOther','CD57+_ImmuneOther','Endothelial','Fibroblast',
             'iNOS+_Pneumocyte','Lung_Epithelium', 'Mesenchymal')

##..Make heatmap matrix..##

hm_allclusters <- matrix(, nrow = length(unique(data$name)), ncol = length(heatmap_channels))
for(i in 1:length(row_order)) {
  temp_mat <- data[data$name==row_order[i], heatmap_channels]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

hm_allclusters
rownames(hm_allclusters) <- paste(row_order)
colnames(hm_allclusters) <- heatmap_channels

##..Get vector of color codes for annotation column..##

colorkey<-read.csv("cellpheno_num_colorkey.csv")
colorkey_cells<-droplevels(colorkey[colorkey$Pheno %in% row_order,])
colorkey_cells$Pheno<-factor(colorkey_cells$Pheno, levels = row_order)
colorkey_cells<-colorkey_cells[order(colorkey_cells$Pheno),]
color<-as.vector(colorkey_cells$Hex)


##..Make heatmap..##
# with breaks between major lineages

meta_data = as.data.frame(row_order)
colnames(meta_data) <-'cell_type' 
meta_data$cell_type<-factor(meta_data$cell_type, levels = row_order)
rownames(meta_data) <- rownames(hm_allclusters)
names(color) <- row_order
mycolors <- list(cell_type = color)

pheatmap(hm_allclusters, 
         cluster_cols = T,
         cluster_rows = F,
         scale = 'none',
         annotation_row = meta_data, 
         annotation_colors = mycolors,
         gaps_row = c(5,11,14,16),
         border_color = 'grey',
         annotation_legend = FALSE,
         color = viridis(256),
         breaks = seq(0,1,length.out=256), 
         filename = "Sub-panels/Figure2/Fig2a_CellPhenotype_Heatmap.pdf")
