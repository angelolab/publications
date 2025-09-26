# MIBIcompareBulkSignalD1MTvControl.R
# Author: Erin McCaffrey 
# Date created: 201129


require(dplyr)
library(ggplot2)
library(tibble)
library(ggpubr)

##..Import data..##

setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv("bulk_intensity.csv")

##..Drop excluded samples..##

drop_samples <- c(17,22,26)
data <- data[!data$SampleID %in% drop_samples,]

##..Add group and subgroup annotation..##

gran_D1MT<-c(8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38) 

data$Group<-'D1MT'

data[data$SampleID %in% gran_con,]$Group<-'control'

##..create melted version..##
data.m<-melt(data, id.vars = c("SampleID",'Group'))


##..Plot and compare..##
ggplot(data = data.m, aes(x = Group, y = value, fill = Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Total Counts') +
  facet_wrap(~variable,scales=c("free"))

plot_data<-data.m[data.m$variable=='IFNg',]
ggplot(data = plot_data, aes(x = Group, y = value, fill = Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method= "t.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Total Counts') 

