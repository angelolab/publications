# TBMIBIcompareCellTypeAbundancesSarcvTB.R
# Author: Erin McCaffrey 
# Date created: 191226
# Overview: This script reads in a dataframe with the counts for all cell types across all
# sarcoid and TB samples. It converts this to frequencies and then plots a paired boxplots for all
# cell types to compare their abundances between thet two sample types. 

library(ggplot2)
library(viridis)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggpubr)

##...Read in the data...##

count_data<-read.csv("data/sarcvTB_counts.csv")


##...Create matrix of frequencies...##

freq_data <- count_data[,-c(23,24)] # remove the diversity index data
freq_data[,2:21] <- freq_data[,2:21]/rowSums(freq_data[,2:21])


##...Create a melted version...##

freq_data_melt <- melt(freq_data, id.vars = c('SampleID', 'Tissue'))
colnames(freq_data_melt)<-c('SampleID', 'Tissue', 'cell_type','Freq')

##...Plot paired boxplot of all cell-types between sarcoid and TB...##

color<-c("#00A59C","#9CD9D3")

signif<-c('CD8_T', 'endothelial','neutrophil','CD16_CD14_Mono', 'giant_cell','epithelial')
plot_data<-droplevels(freq_data_melt[freq_data_melt$cell_type %in% signif,])

# plot paired boxplots with stats 

ggplot(data = freq_data_melt, aes(x = fct_reorder(as.factor(cell_type),Freq,.fun=median,.desc=TRUE), Freq, fill=Tissue)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Tissue), method= "wilcox.test", label = "p.signif") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = color)  +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 


# plot individual boxplots for all subsets

ggplot(data = freq_data_melt, aes(Tissue, Freq, fill=Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitterdodge()) +
  stat_compare_means(aes(group = Tissue), method= "wilcox.test", label = "p.signif") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = color)  +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) +
  facet_wrap(~cell_type,scales=c("free")) 

#...Determine FC of mean freq between sarcoid and TB...##

data_summary <- freq_data %>%
  group_by(Tissue) %>%
  summarize_all(list(mean=mean))

data_summary<-as.data.frame(t(data_summary))
names(data_summary)<-c('TB','Sarcoid')
data_summary<-droplevels(data_summary[-c(1,2),])
data_summary$cell_type <- rownames(data_summary)
data_summary$TB<-as.numeric(as.character(data_summary$TB))
data_summary$Sarcoid<-as.numeric(as.character(data_summary$Sarcoid))
data_summary$FC<-log2(data_summary$TB/data_summary$Sarcoid)


ggplot(data=data_summary[is.finite(data_summary$FC),], aes(x=reorder(cell_type, -FC), y=FC, fill=FC)) +
  geom_bar(stat="Identity") +
  scale_fill_viridis() +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  labs(y="log2(mean freq TB / mean freq sarc)") 


