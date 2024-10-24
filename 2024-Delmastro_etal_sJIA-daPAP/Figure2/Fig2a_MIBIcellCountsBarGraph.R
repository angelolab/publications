# Fig2a_MIBIcellCountsBarGraph.R
# Author: Erin McCaffrey 
# Date updated: 07/15/2024
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it plots a bar graph of the total counts for each cell type in the same order as the heatmap for visualization
# purposes.

library(ggplot2)

##..Import data..##
setwd("/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points")
data<-read.csv("celltable_05032024.csv")

##..Specify order of populations..##

order<-c('Bcell','CD4+_Tcell','CD8+_Tcell','CD57+_CD8+_Tcell','Treg','CD11c+_mDC','CD14+_Mono','CD209+_Mac','iNOS+_Mac','M2_Mac',
          'Giant_cell','Eosinophil', 'Mast_cell','Neutrophil', 'CD16+_ImmuneOther','CD57+_ImmuneOther','Endothelial','Fibroblast',
         'iNOS+_Pneumocyte','Lung_Epithelium', 'Mesenchymal')


##..Create table of cell type and total number..##

cell_counts<-as.data.frame(table(data$name))
colnames(cell_counts)<-c('name','count')

##..Reorder to match specified order..##

cell_counts$name<-factor(cell_counts$name, levels = order)
cell_counts<-cell_counts[order(cell_counts$name),]

##..Add group for faceting..##

cell_counts$group = 5
cell_counts[1:5,]$group = 1
cell_counts[6:11,]$group = 2
cell_counts[12:14,]$group = 3
cell_counts[15:16,]$group = 4

##..Get vector of color codes for annotation column..##

colorkey<-read.csv("cellpheno_num_colorkey.csv")
colorkey_cells<-droplevels(colorkey[colorkey$Pheno %in% order,])
colorkey_cells$Pheno<-factor(colorkey_cells$Pheno, levels = order)
colorkey_cells<-colorkey_cells[order(colorkey_cells$Pheno),]
color<-as.vector(colorkey_cells$Hex)

##..Produce barplot with counts in decreasing order..##

ggplot(data=cell_counts, aes(x=name, y=count, fill = name)) +
  geom_bar(position="dodge",stat="identity") + 
  scale_fill_manual(values = color) +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = 'None')

ggsave('Sub-panels/Figure2/Fig2a_CellPhenotype_BarChart.pdf')
