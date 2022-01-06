# TBMIBIpulmvextrapulm.R
# Author: Erin McCaffrey 
# Date created: 210707
# Overview: This script reads in the frequency data for granuloma immune cell population. Itthen plots the proportion 
# of total cells of all cell types in extra-pulmonary v pulmonary samples to determine organ-site specific cell traits

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(ggpmisc)

##..Data importation for immune frequencies..## 

setwd("/Volumes/GoogleDrive/My Drive/AngeloLab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("immune_cell_freqs_noHIV.csv")

##..Drop sarcoid and add organ site annotation..##

data_TB<-data[data$Tissue == 'mycobacteria',]
pulm<-c(21,84,42,88,28,89,14,15,98,99,90,91,94,95,96,97)
expulm<-c(6,7,33,34,26,27,40,61,47,48,54,55,92,93)

data_TB$site<-'pulmonary'
data_TB[data_TB$SampleID %in% expulm,]$site <- 'extra-pulm'

##..Plot across groups..##

ggplot(data = data_TB, aes(site, frequency, fill=site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitterdodge()) +
  stat_compare_means(aes(group = site), method= "wilcox.test", label = "p.format") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) +
  facet_wrap(~cell_type,scales=c("free")) 

##..Make single plot..##

ggplot(data = data_TB, aes(x = fct_reorder(as.factor(cell_type),frequency,.fun=median,.desc=TRUE), frequency, fill = site)) +
  geom_boxplot() +
  stat_compare_means(aes(group = site), method= "wilcox.test", label = "p.signif") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

##..Data importation for lineage frequencies..## 

data<-read.csv("lineage_freqs_noHIV.csv")

##..Drop sarcoid and add organ site annotation..##

data_TB<-data[data$Tissue == 'mycobacteria',]
pulm<-c(21,84,42,88,28,89,14,15,98,99,90,91,94,95,96,97)
expulm<-c(6,7,33,34,26,27,40,61,47,48,54,55,92,93)

data_TB$site<-'pulmonary'
data_TB[data_TB$SampleID %in% expulm,]$site <- 'extra-pulm'

##..Plot across groups..##

ggplot(data = data_TB, aes(site, frequency, fill=site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitterdodge()) +
  stat_compare_means(aes(group = site), method= "wilcox.test", label = "p.format") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Cell Lineage') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) +
  facet_wrap(~lineage,scales=c("free")) 

##..Make single plot..##

ggplot(data = data_TB, aes(x = fct_reorder(as.factor(lineage),frequency,.fun=median,.desc=TRUE), frequency, fill = site)) +
  geom_boxplot() +
  stat_compare_means(aes(group = site), method= "wilcox.test", label = "p.signif") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Cell Lineage') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 


