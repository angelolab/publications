# PAP_MIBIFunctionalheatmapPerSubset.R
# Author: Erin McCaffrey 
# Date created: 220110

require(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(viridis)
library(colorspace)
library(pals)
library(pheatmap)

##..Load in Data..##

setwd("/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points")
data <- read.csv("celltable_05032024.csv")
fov_labels = read.csv("fov_labels.csv")

##..Define Disease Category of Interest..##

cats <- c("sJIA-daPAP", "non-sJIA-PAP")

##..Merge topic loadings with cell table..##

data <- merge(data, fov_labels, by = c("point"))

##..Define Cell Phenotypes and Functional Markers..##

cell_order = c('Bcell', "CD4+_Tcell", "CD8+_Tcell", "CD57+_CD8+_Tcell", "Treg", "CD11c+_mDC", "CD14+_Mono", "CD209+_Mac", "iNOS+_Mac", "M2_Mac", "Eosinophil", "Mast_cell", "Neutrophil", "CD16+_ImmuneOther", "CD57+_ImmuneOther", "Endothelial", "Fibroblast", "iNOS+_Pneumocyte", "Lung_Epithelium", "Mesenchymal")
markers = c('IFNg','HLA.DR', 'CD45RO','TIM3','pS6','CD57','MMP9','Ki67','GrzB','HO.1','IDO','iNOS','H3K9Ac','H3K27me3')
lst = c("name", markers)

##..Iterate through disease categories and Concatenate tables..##
cell_all <- data.frame(matrix(NA, nrow = length(markers), ncol = 0))

for (cat in cats) {
  ##..Lipoproteinosis cells..##
  
  cell_Lipo <- data %>% dplyr::filter(lipoprotein == 1 & FOV_Category == cat) %>% 
    dplyr::select(all_of(lst)) %>% 
    dplyr::mutate(across(-name, ~ . + 0.01)) %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarize_all(list(mean = mean))
  colnames(cell_Lipo) <- append("name",markers)
  
  ##..Non-Lipoproteinosis cells..##
  
  cell_noLipo <- data %>% dplyr::filter(lipoprotein == 0 & FOV_Category == cat) %>% 
    dplyr::select(all_of(lst)) %>% 
    dplyr::mutate(across(-name, ~ . + 0.01)) %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarize_all(list(mean = mean)) %>%
    dplyr::filter(name != "Giant_cell")
  colnames(cell_noLipo) <- append("name",markers)
  
  orig_rowNames <- cell_Lipo$name
  
  ##..Calculate Log2 Fold Change of Average Expression..##
  
  cell_FC = cell_Lipo[,markers] / cell_noLipo[,markers]
  cell_log2FC = log2(cell_FC)
  
  ##..Transpose and Reorder Columns by Phenotype Groupings..##
  
  cell_log2FC = t(cell_log2FC)
  colnames(cell_log2FC) <- orig_rowNames
  cell_log2FC = cell_log2FC[,cell_order]
  
  ##..Concatenate..##
  
  cell_all <- cbind(cell_all, cell_log2FC)
}

##..Plot Heatmap..##

z_max <- 2
paletteLength <- 100
myColor <- diverging_hcl(paletteLength, palette = 'Blue-Red 3')
myBreaks <- c(seq(-z_max, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0 + (z_max/(paletteLength/2)), z_max, length.out=floor(paletteLength/2)))
pheatmap(cell_all, 
         color=myColor, 
         breaks=myBreaks, 
         cluster_cols = FALSE, 
         cellheight=10, 
         cellwidth = 10, 
         gaps_col = c(20),
         filename=paste("Sub-panels/SupplementalFigures/SuppFig5d_Lipoprotein_FC_heatmap.pdf", sep = ""))

