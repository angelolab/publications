# MIBIevalCellDistributionAcrossRegionsPerPoint.R
# Author: Erin McCaffrey 
# Date created: 200310
# Overview: Script reads in a csv of all cell data with vessel assignment

library(ggplot2)
library(tidyverse)
library(dplyr)

##..Import data..##
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv("celldata_region_annotated.csv")

##..Drop excluded samples..##
drop_samples <- c(17,22,26)
data <- data[!data$SampleID %in% drop_samples,]

##..Mutate matrix to have all cell types and the region per point (not patient)..##

data_summ <- as.data.frame(table(data$SampleID,data$cell_type, 
                                     data$region))
colnames(data_summ) <- c("SampleID","cell_type","region","count")

##..Determine the frequency of a given cell type in that region out of all cells of that type in sample..##

#get cell totals for each type across samples
cell_totals <- as.data.frame(table(data$SampleID,data$cell_type))
colnames(cell_totals) <- c("SampleID","cell_type","total")

#sort and merge dataframes
cell_totals <- cell_totals %>% slice(rep(1:n(), each = 2))
cell_totals <- cell_totals[order(cell_totals$SampleID), ] 
cell_totals <- cell_totals[order(cell_totals$cell_type), ] 
data_summ <- data_summ[order(data_summ$SampleID), ]
data_summ <- data_summ[order(data_summ$cell_type), ] 
data_summ$cell_total <- cell_totals$total

#determine frequency per sample per region
data_summ$freq <- data_summ$count / data_summ$cell_total
data_summ$region <- factor(data_summ$region, 
                               levels = c("core","peripheral"))

#add annotation for sample type
gran_D1MT<-c(8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,23,24,28,32,33,34,35,36,37,38)
data_summ$group<-'gran_D1MT'
data_summ[data_summ$SampleID %in% gran_con, ]$group <- 'gran_con'

##..Plot summary stats..##

ggplot(data = data_summ, aes(x = region, y = freq)) + 
  geom_boxplot(aes(fill = cell_type)) +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y')

ggplot(data = data_summ, aes(x = region, y = freq)) + 
  geom_boxplot(aes(fill = group)) +
  stat_compare_means(aes(group = group),label = "p.signif", method= "wilcox.test") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  scale_y_continuous(limits=c(0,1.5)) +
  facet_wrap(.~cell_type, scales='free_y')


plot_data<-data_summ[data_summ$cell_type=='CD8_T',]
ggplot(data = plot_data, aes(x = region, y = freq)) + 
  geom_boxplot(aes(fill = group)) +
  stat_compare_means(aes(group = group),label = "p.signif", method= "wilcox.test") +
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') 




