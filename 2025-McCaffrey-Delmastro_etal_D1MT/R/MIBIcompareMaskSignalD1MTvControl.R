# MIBIcompareMaskSignalD1MTvControl.R
# Author: Erin McCaffrey 
# Date created: 201129


require(dplyr)
library(ggplot2)
library(tibble)
library(reshape2)
library(ggpubr)

##..Import data_mask..##

setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell/")
data_mask<-read.csv("mask_intensity.csv")


##..Add group and subgroup annotation..##

gran_D1MT<-c(8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38) 

data_mask$Group<-'D1MT'
data_mask[data_mask$SampleID %in% gran_con,]$Group<-'control'

##..create melted version..##
data_mask.m<-melt(data_mask, id.vars = c("SampleID",'Group'))


##..Plot and compare..##
ggplot(data = data_mask.m, aes(x = Group, y = value, fill = Group)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Total Counts') +
  facet_wrap(~variable,scales=c("free"))

