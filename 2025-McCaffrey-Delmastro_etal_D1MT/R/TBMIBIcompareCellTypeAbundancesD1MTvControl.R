# TBMIBIcompareCellTypeAbundancesD1MTvControl.R
# Author: Erin McCaffrey 
# Date created: 201124

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

##...Read in the data...##

setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
freq_data<-read.csv("lineage_freqs.csv")

##..Drop the excluded samples..##
drop_samples <- c(17,22,26)
freq_data <- freq_data[!freq_data$SampleID %in% drop_samples,]

##...Plot paired boxplot of all cell-types between sarcoid and TB...##

# plot paired boxplots with stats 
ggplot(data = freq_data, aes(x = fct_reorder(as.factor(lineage),frequency,.fun=median,.desc=TRUE), frequency, fill=Tissue)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Tissue), method= "wilcox.test", label = "p.format") +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

ggplot(data = freq_data, aes(x = fct_reorder(as.factor(lineage),count,.fun=median,.desc=TRUE), count, fill=Tissue)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Tissue), method= "wilcox.test", label = "p.format") +
  theme_bw() +
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

