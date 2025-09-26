# TBMIBIcompareMixingScore
# Author: Erin McCaffrey 
# Date created: 201124
# Overview: This script reads in the mixing score for the NHP D1-MT and control, and human cohort. 
# It merges the datasets and then between each group compares the mixing score and number of 
# myeloid-lymphocyte interactions.

require(gtools)
library(ggplot2)
library(dplyr)
library(ggpubr)

##..Read in NHP and human data..##

# setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
# human_mixing<-read.csv("mixingscore_10um.csv")

setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
nhp_mixing<-read.csv("mixingscore_10um_CD4.csv")

##..Drop excluded samples..##

drop_samples <- c(17,22,26)
nhp_mixing <- nhp_mixing[!nhp_mixing$SampleID %in% drop_samples,]


##..Add species data..##

# human_mixing$species<-'human'
nhp_mixing$species<-'nhp'

##..Assign tissue type to sample IDs..##

gran_D1MT<-c(7,8,9,13,14,15,16,17,18,25,26,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38) 
# AHRI<-c(64,65,21,84,42,88,28,89,85,13,35,36)
# Stanford<-c(14,15,98,99,6,7,33,34,26,27,40,61,47,48,54,55,92,93)
# Autopsy<-c(90,91,94,95,96,97)
# 
# human_mixing$tissue<-NA
# human_mixing[human_mixing$SampleID %in% AHRI,]$tissue<-'Resection'
# human_mixing[human_mixing$SampleID %in% Stanford,]$tissue<-'Biopsy'
# human_mixing[human_mixing$SampleID %in% Autopsy,]$tissue<-'Autopsy'

nhp_mixing$tissue<-NA
nhp_mixing[nhp_mixing$SampleID %in% gran_D1MT,]$tissue<-'IDO_inhibitor'
nhp_mixing[nhp_mixing$SampleID %in% gran_con,]$tissue<-'control'

##..Merge into single dataset..##

# combo_mixing<-rbind(nhp_mixing, human_mixing)

##..Plot NHP only..##

# Mixing score
data_summary <- nhp_mixing %>%
  group_by(tissue) %>%
  summarize(combo_median = median(Mixing_Score))

ggplot(data=nhp_mixing, aes(x=tissue, y=Mixing_Score)) +
  geom_point(data=nhp_mixing, aes(x=tissue, y=Mixing_Score), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Mixing Score")

ggplot(data=nhp_mixing, aes(x=tissue, y=Mixing_Score, fill= tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Mixing Score")

# Myeloid-T cell interactions
data_summary <- nhp_mixing %>%
  group_by(tissue) %>%
  summarize(combo_median = median(n_Myeloid.Lymph))

ggplot(data=nhp_mixing, aes(x=tissue, y=n_Myeloid.Lymph)) +
  geom_point(data=nhp_mixing, aes(x=tissue, y=n_Myeloid.Lymph), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="n_Myeloid.Lymph")

##..Plot All..##

my_comparisons = list(c('Resection','Biopsy'),
                      c('Resection','Autopsy'),
                      c('Autopsy','Biopsy'),
                      c('IDO_inhibitor','control'))

data_summary <- combo_mixing %>%
  group_by(tissue) %>%
  summarize(combo_median = median(Mixing_Score))

ggplot(data=combo_mixing, aes(x=tissue, y=Mixing_Score)) +
  geom_point(data=combo_mixing, aes(x=tissue, y=Mixing_Score), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Mixing Score")







