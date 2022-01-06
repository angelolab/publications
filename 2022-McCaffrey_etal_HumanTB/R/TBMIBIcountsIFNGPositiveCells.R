# TBMIBIcountsIFNGPositiveCells.R
# Author: Erin McCaffrey
# Date created: 201117
# Overview: This script reads in the csv for the single cell data and counts the number of IFNG+ cells per patient 
# and plots it as a bar plot. 

library(ggplot2)
library(forcats)
library(dplyr)

##...Load in data..##

data<-read.csv('data/allTB-sarcoid-scdata.csv')

##..Counts of all IFNG positive cells across cell types..##

IFNG_pos<-as.data.frame(table(data[data$IFNg>0.43, ]$SampleID, data[data$IFNg>0.43, ]$cell_type))
colnames(IFNG_pos)<-c('SampleID','cell_type','Count')

colorkey<-read.csv('data/colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% unique(IFNG_pos$cell_type),])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = unique(IFNG_pos$cell_type))
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color<-as.vector(colorkey_imm$code)

ggplot(data = IFNG_pos, aes(x = as.factor(SampleID), y = Count, fill = cell_type)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = color) +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = 'Sample') + 
  labs(y = 'Count IFNG+ Cells')


