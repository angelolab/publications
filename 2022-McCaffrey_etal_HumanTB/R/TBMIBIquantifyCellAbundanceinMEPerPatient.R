# TBMIBIquantifyCellAbundanceinMEPerPatient.R
# Author: Erin McCaffrey
# Date created: 191205
# Overview: This script reads in the csv for the ME loadings for each cell. It assigns each cell to an ME 
# (max loading) and then plots the abundance of each cell type in each topic across each patient 
# (ie. one subplot per ME with 1 stacked bar per patient)

library(ggplot2)
library(dplyr)

##...Load in data..##
ME_data<-read.csv('data/allTB-ME_annotated.csv')

##...Quantify cell type abundances in each ME...##

ME_freqs<-as.data.frame(table(ME_data$cell_type, ME_data$maxME, ME_data$SampleID))
colnames(ME_freqs)<-c('cell_type','ME','SampleID','count')

##...Plot broken down by cell type..##

imm_order <- unique(ME_freqs$cell_type)
colorkey<-read.csv('data/colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% imm_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = imm_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color<-as.vector(colorkey_imm$code)

point_order_gran<-c(21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,26,27,40,61,47,48,54,55,92,93)
ME_freqs$SampleID <- factor(ME_freqs$SampleID, levels=point_order_gran)

ggplot(data = ME_freqs, aes(x = SampleID, y = count, fill = cell_type)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x = 'SampleID') + 
  labs(y = 'Cell Count') +
  scale_fill_manual(name = 'ME', values = color) +
  facet_wrap(.~ME, scale='free_y', ncol = 2)

