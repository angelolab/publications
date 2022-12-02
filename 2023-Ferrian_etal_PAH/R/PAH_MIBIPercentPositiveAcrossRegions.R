# PAH_MIBIPercentPositiveAcrossRegions.R
# Author: Erin McCaffrey 
# Date created: 200410
# Overview: This script reads in the csv for all PAH cell data with the percent positive for all
# markers across all cell subsets per patient. It next summarizes hlt v pah and ipah v hpah for all
# subsets and markers

library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggpubr)


##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("allpatients_freq-of-positivecells-perregion.csv")

##..Define functional markers and colors..##

functional_markers<-colnames(data[,c(4,7,9,12,19,23,24,25,26,27,28,29,30,34,35,37)])

color<-c('#651E38','#C57050','#C8C8C7')

##..Change the format of the dataframe..##

data.m<-melt(data, id.vars = c('region_modified','PID','Subgroup'), measure.vars = functional_markers)

##..Plot pah all markers all zones..##

my_comparisons<-list(c("vessel_associated","non_vascular"),
                     c("vessel_associated","vessel_proximal"),
                     c("non_vascular","vessel_proximal"))

plot_data<-droplevels(data.m[data.m$Subgroup %in% c('ipah','hpah'),])

ggplot(data = plot_data, aes(x = region_modified, y = value)) + 
  geom_boxplot(aes(fill = region_modified), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~variable, scales='free_y', ncol=4) + 
  theme(legend.position = 'none')

##..Plot pah all markers all zones..##

plot_data<-droplevels(data.m[data.m$Subgroup %in% c('hlt'),])

ggplot(data = plot_data, aes(x = region_modified, y = value)) + 
  geom_boxplot(aes(fill = region_modified), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~variable, scales='free_y', ncol=4) + 
  theme(legend.position = 'none')

##..Plot ipah v hpah all markers all zones..##

plot_data<-droplevels(data.m[data.m$Subgroup %in% c('hpah','ipah') & data.m$region_modified == 'non_vascular',])

ggplot(data = plot_data, aes(x = Subgroup, y = value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~variable, scales='free_y', ncol=4) + 
  theme(legend.position = 'none')
