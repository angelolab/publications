# PAH_immuneCellFreqsPerPatient.R
# Author: Erin McCaffrey 
# Date created: 200307
# Overview: This script reads in the csv for the PAH dataset, determines the percent of 
# total for all patients and produces a stacked bar chart


library(dplyr)
library(viridis)
library(ggplot2)
library(forcats)
library(reshape2)
library(tidyr)
library(ggpubr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")

##..Subset just the immune cells..##

data_immune<-droplevels(data[!data$cell_lineage %in% c('fibroblast','endothelial','mesenchymal','epithelial'),])


##..Create a table of cell frequencies per patient..##

imm_freqs<-as.data.frame(table(data_immune$PID,data_immune$cell_lineage))
colnames(imm_freqs)<-c('Patient','Cell_Type','Count')
totals_patient<-as.numeric(table(data_immune$PID))
imm_freqs$total<-rep(totals_patient,8)
imm_freqs$freq<-as.numeric(imm_freqs$Count/imm_freqs$total)

##..Order and produce color key..##

order<-c('Tc','Th','Mono','DC','Neutro','Macro','NK','Bcell')

imm_freqs$Cell_Type <- factor(imm_freqs$Cell_Type, levels=rev(order))
imm_freqs<-imm_freqs[order(imm_freqs$Cell_Type),]

cell_color_key<-read.csv('colorkey.csv')
cell_color_key<-droplevels(cell_color_key[cell_color_key$cell_lineage %in% order, ])
cell_color_key$cell_lineage<-factor(cell_color_key$cell_lineage, levels = order)
cell_color_key<-cell_color_key[order(cell_color_key$cell_lineage),]
color<-as.vector(cell_color_key$codes)

##..Stacked bars all patients..##

ggplot(imm_freqs, aes(x=Patient, y=freq, fill=Cell_Type)) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
