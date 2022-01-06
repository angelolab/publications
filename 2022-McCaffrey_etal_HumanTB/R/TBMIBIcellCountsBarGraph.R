# TBMIBIcellCountsBarGraph.R
# Author: Erin McCaffrey 
# Date created: 190830
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it plots a bar graph of the total counts for each cell type in the same order as the heatmap for visualization
# purposes.

library(ggplot2)

##..Import data..##
data_gran<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Specify order of populations..##

order<-c("CD4_T","B_cell","CD11b/c_CD206_Mac/Mono","CD8_T","CD14_Mono","CD11c_DC/Mono","imm_other","CD16_CD14_Mono","endothelial","fibroblast",
         "CD68_Mac","CD163_Mac","CD206_Mac","neutrophil","epithelial","CD209_DC","Treg","mast","gdT_cell","giant_cell")
         

##..Create table of cell type and total number..##

cell_counts<-as.data.frame(table(data_gran$cell_type))

##..Produce barplot with counts in decreasing order..##

ggplot(data=cell_counts, aes(x=reorder(Var1, Freq),y=Freq)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  theme_bw() + 
  labs(x='Count') + 
  labs(y='Cell Type')
