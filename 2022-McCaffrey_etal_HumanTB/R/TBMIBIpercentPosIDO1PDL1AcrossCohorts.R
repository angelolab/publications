# TBMIBIpercentPosIDO1PDL1AcrossCohorts.R
# Author: Erin McCaffrey 
# Date created: 200317
# Overview: This script reads in the normalized data for IDO1 and PDL1. Next it plots the frequency of
# positive cells for IDO1 and PDL1 across the cohort broken down by sample type. 

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(ggpubr)

##..Import data, get just myeloid cells in TB and sarcoid..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")
sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
data<-droplevels(data[!data$SampleID %in% c(sarcoid), ])
data<-droplevels(data[data$cell_lin == "myeloid", ])

##..Get the percent IDO and PDL1 pos cells in each point..#

IDO_thresh = 0.26
PDL1_thresh = 0.25

data_IDO<-droplevels(data[data$IDO>=IDO_thresh, ])
data_PDL1<-droplevels(data[data$PD.L1>=PDL1_thresh, ])

##..Across samples..##

freq_sample<-as.data.frame(table(data$SampleID))

#IDO
freq_IDO_sample<-as.data.frame(table(data_IDO$SampleID))
freq_sample$IDOpos<-freq_IDO_sample$Freq
names(freq_sample)<-c("SampleID","Total","Total_IDOpos")
freq_sample$percentIDO<-as.numeric(format((freq_sample$Total_IDOpos / freq_sample$Total)*100),digits=3)

#PDL1
freq_PDL1_sample<-as.data.frame(table(data_PDL1$SampleID))
freq_sample$Total_PDL1pos<-0
freq_sample[freq_sample$SampleID %in% freq_PDL1_sample$Var1,]$Total_PDL1pos<-freq_PDL1_sample$Freq
freq_sample$percentPDL1<-as.numeric(format((freq_sample$Total_PDL1pos / freq_sample$Total)*100),digits=3)

##..Add sample type annotation..##

# Define categories

resection<-c(21,84,42,88,28,89)
diagnostic_pulm<-c(14,15,98,99)
diagnostic_expulm<-c(6,7,33,34,26,27,40,61,47,48,54,55,92,93)
autopsy<-c(90,91,94,95,96,97)

freq_sample$Tissue<-'lung_resection'
freq_sample[freq_sample$SampleID %in% c(diagnostic_pulm,diagnostic_expulm),]$Tissue<-'biopsy'
freq_sample[freq_sample$SampleID %in% autopsy,]$Tissue<-'autopsy'

freq_sample$Organ<-'lung'
freq_sample[freq_sample$SampleID %in% diagnostic_expulm,]$Organ<-'extrapulm'

##..Plot..##

my_comparisons <- list(c("biopsy","lung_resection"),
                       c("biopsy","autopsy"))

#IDO1

ggplot(data = freq_sample, aes(x = Tissue, y = percentIDO, fill = Tissue)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3, position = position_jitterdodge(), aes(shape=Organ)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency')


#PDL1

ggplot(data = freq_sample, aes(x = Tissue, y = percentPDL1, fill = Tissue)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3, position = position_jitterdodge(), aes(shape=Organ)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency')

