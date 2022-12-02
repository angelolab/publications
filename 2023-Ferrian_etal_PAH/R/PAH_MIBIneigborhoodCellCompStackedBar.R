# PAH_MIBIneigborhoodCellCompStackedBar.R
# Author: Erin McCaffrey 

require(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(viridis)
library(colorspace)
library(pals)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data_neighborhood<-read.csv("PAH_AllCells_Neighborhood_K=10_PatientAnnotated.csv")

##..Create a table of cell type counts per cluster..##

neighborhood_freqs<-as.data.frame(table(data_neighborhood$cell_type, data_neighborhood$cluster))
colnames(neighborhood_freqs)<-c('cell_type','NC','Count')
totals_cell<-as.numeric(table(data_neighborhood$cell_type))
neighborhood_freqs$celltype_total<-rep(totals_cell,10)
neighborhood_freqs$freq_of_celltype<-as.numeric(neighborhood_freqs$Count/neighborhood_freqs$celltype_total)

##..Reorder cell types..##

cell_types<-c('Th','Tc','Bcell','NK','Macro','Mono','Neutro','DC','epithelial','endothelial','mesenchymal','fibroblast')
neighborhood_freqs$cell_type<-factor(neighborhood_freqs$cell_type, levels = cell_types)
neighborhood_freqs<-neighborhood_freqs[order(neighborhood_freqs$cell_type),]

##..Stacked bar all clusters..##

ggplot(neighborhood_freqs, aes(x=cell_type, y=freq_of_celltype, fill=NC)) + 
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") 

##..Create a table of cluster frequency per cell type..##

cell_freqs<-as.data.frame(table(data_neighborhood$cluster, data_neighborhood$cell_type))
colnames(cell_freqs)<-c('NC','cell_type','Count')
totals_NC<-as.numeric(table(data_neighborhood$cluster))
cell_freqs$NC_total<-rep(totals_NC,12)
cell_freqs$freq_of_NC<-as.numeric(cell_freqs$Count/cell_freqs$NC_total)

##..Reorder cell types..##

cell_types<-c('Th','Tc','Bcell','NK','Macro','Mono','Neutro','DC','epithelial','endothelial','mesenchymal','fibroblast')
cell_freqs$cell_type<-factor(cell_freqs$cell_type, levels = cell_types)
cell_freqs<-cell_freqs[order(cell_freqs$cell_type),]

##..Read and order color key..##

colorkey<-read.csv('colorkey.csv')
colorkey$cell_lineage<-factor(colorkey$cell_lineage, levels = cell_types)
colorkey<-colorkey[order(colorkey$cell_lineage),]
color<-as.vector(colorkey$codes)

##..Stacked bar all clusters..##

ggplot(cell_freqs, aes(x=NC, y=freq_of_NC, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values=color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") 
