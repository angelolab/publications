# PAH_MIBIPercentPositiveAcrossGroups.R
# Author: Erin McCaffrey 
# Date created: 200410
# Overview: This script reads in the csv for all PAH cell data with the percent positive for all
# markers across all cell subsets per patient. It next summarizes hlt v pah and ipah v hpah for all
# subsets and markers

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)


##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("allpatients_freq-of-positivecells-persubset.csv")

##..Reorder cell types and define functional markers..##

order<-c('Th','Tc','Bcell','NK','Macro','Mono','Neutro','DC','endothelial','epithelial',
         'mesenchymal','fibroblast')
data$cell_lineage<-factor(data$cell_lineage, levels=order)

functional_markers<-colnames(data[,c(4,7,9,12,19,23,24,25,26,27,28,29,30,34,35,37)])


##..Load color key..##

colorkey<-read.csv('colorkey.csv')
colorkey$cell_lineage<-factor(colorkey$cell_lineage, levels = levels(data$cell_lineage))
colorkey<-colorkey[order(colorkey$cell_lineage),]
color<-as.vector(colorkey$codes)

##..Plot hlt v pah..##

ggplot(data = data, aes(x = Group, y = TIM.3)) + 
  geom_boxplot(aes(fill = cell_lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_lineage, scales='free_y', ncol=4) + 
  theme(legend.position = 'none')


##..Plot ipah v hpah..##
plot_data<-droplevels(data[!data$Group=='hlt',])
ggplot(data = plot_data, aes(x = Subgroup, y = TIM.3)) + 
  geom_boxplot(aes(fill = cell_lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Condition') + 
  labs(y = 'Percent Positive') +
  facet_wrap(.~cell_lineage, scales='free_y', ncol=4) + 
  theme(legend.position = 'none')

