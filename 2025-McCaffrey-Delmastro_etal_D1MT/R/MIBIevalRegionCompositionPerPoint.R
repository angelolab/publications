# MIBIevalRegionCompositionPerPoint.R
# Author: Erin McCaffrey 
# Date created: 200310
# Overview: Script reads in a csv of all cell data with vessel assignment

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

##..Import data..##
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv("celldata_region_annotated.csv")

##..Drop excluded samples..##
drop_samples <- c(17,22,26)
data <- data[!data$SampleID %in% drop_samples,]

##..Determine the frequency of a given cell type in that region out of all cells in that region..##

data_summ <- as.data.frame(table(data$cell_type, data$region, data$SampleID))
colnames(data_summ) <- c("cell_type","region","SampleID","count")

#get cell totals for each region across samples
region_totals <- as.data.frame(table(data$SampleID,data$region))
colnames(region_totals) <- c("SampleID","region","total")

#sort and merge dataframes
region_totals <- region_totals %>% slice(rep(1:n(), each = 13))
region_totals <- region_totals[order(region_totals$SampleID), ] 
region_totals <- region_totals[order(region_totals$region), ] 
data_summ <- data_summ[order(data_summ$SampleID), ]
data_summ <- data_summ[order(data_summ$region), ] 
data_summ$region_total <- region_totals$total

#determine frequency per sample per region
data_summ$freq <- data_summ$count / data_summ$region_total

#add annotation for sample type
gran_D1MT<-c(8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,23,24,28,32,33,34,35,36,37,38)
data_summ$group<-'gran_D1MT'
data_summ[data_summ$SampleID %in% gran_con, ]$group <- 'gran_con'

##..Plot..##

ggplot(data = data_summ, aes(x = region, y = freq)) + 
  geom_boxplot() +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y', ncol = 7)

ggplot(data = data_summ, aes(x = region, y = freq)) + 
  geom_boxplot(aes(fill = group)) +
  stat_compare_means(aes(group = group),label = "p.signif", method= "wilcox.test") +
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y', ncol = 7)

plot_data<-data_summ[data_summ$cell_type=='CD8_T',]
ggplot(data = plot_data, aes(x = region, y = freq)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = group)) +
  geom_point(size = 2, aes(fill = group), position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = group),label = "p.format", method= "wilcox.test") +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') 

##..For CD8 freq get log2 FC for each point from median and mean control..##

median_con_core<-median(plot_data[plot_data$group=="gran_con" & plot_data$region=='core',]$freq, na.rm = T)
mean_con_core<-mean(plot_data[plot_data$group=="gran_con" & plot_data$region=='core',]$freq, na.rm = T)
median_con_peri<-median(plot_data[plot_data$group=="gran_con" & plot_data$region=='peripheral',]$freq, na.rm = T)
mean_con_peri<-mean(plot_data[plot_data$group=="gran_con" & plot_data$region=='peripheral',]$freq, na.rm = T)

CD8_D1MT<-plot_data[plot_data$group=="gran_D1MT" ,]
CD8_D1MT$FC_median<-NA
CD8_D1MT$FC_mean<-NA
CD8_D1MT[CD8_D1MT$region=='core',]$FC_median<-log2(CD8_D1MT[CD8_D1MT$region=='core',]$freq/median_con_core)
CD8_D1MT[CD8_D1MT$region=='core',]$FC_mean<-log2(CD8_D1MT[CD8_D1MT$region=='core',]$freq/mean_con_core)
CD8_D1MT[CD8_D1MT$region=='peripheral',]$FC_median<-log2(CD8_D1MT[CD8_D1MT$region=='peripheral',]$freq/median_con_peri)
CD8_D1MT[CD8_D1MT$region=='peripheral',]$FC_mean<-log2(CD8_D1MT[CD8_D1MT$region=='peripheral',]$freq/mean_con_peri)

ggplot(data = CD8_D1MT, aes(x = region, y = FC_median)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, aes(fill = group), position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency')


