# MIBI_neutrophil_radial_analysis.R
# Date created: 07/28/2024
# This script takes the output of the neutrophil radial classification
# and compares distributions based on CFU status. 

library(reshape2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

##..Read in the data..##

setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
neut_data <- read.csv('./spatial_analysis/neutrophil/Neutrophil_ratio_25um.csv')
meta_data <- read.csv('./cohort_metadata/study_cohort_metadata.csv')

##..Append metadata to data..##

samples <- unique(neut_data$sample)
meta_data_sub <- meta_data[meta_data$sample %in% samples,]
neut_data_anno <- left_join(neut_data, meta_data_sub, by = c('sample'))

##..Plot and compare..##

plot_data <- neut_data_anno[is.finite(neut_data_anno$log2FC),]
  
ggplot(plot_data, aes(burden, log2FC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", width = 0.1) +
  theme_bw() +
  stat_compare_means(method= "t.test") +
  theme(legend.position = 'none')


