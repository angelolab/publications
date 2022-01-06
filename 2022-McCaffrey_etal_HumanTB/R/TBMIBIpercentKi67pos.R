# TBMIBIpercentKi67pos.R
# Author: Erin McCaffrey 
# Date created: 190501
# Overview: This script reads in the csv for annotated, cell-size normalized, scaled 
# expression data. It asinh transforms the data, then looks at the percent of subsets that
# are Ki67 positive. 

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(exactRankTests)

##..Import data, subset TB only..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")
data_myco<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])

##..Get just the lymphocytes..##

lymphs<-c("B_cell","CD4_T","CD8_T","Treg")

##..Get the percent Ki67 pos cells in each subset..#

data_lymph<-droplevels(data_myco[data_myco$cell_type %in% lymphs, ])
data_lymph_Ki67<-droplevels(data_lymph[data_lymph$Ki67>0.27578920, ])

#bulk
freq<-as.data.frame(table(data_lymph$cell_type))
freq_Ki67<-as.data.frame(table(data_lymph_Ki67$cell_type))
freq$Ki67pos<-freq_Ki67$Freq
names(freq)<-c("Cell_Type","Total","Total_Ki67pos")
freq$percent<-as.numeric(format((freq$Total_Ki67pos/freq$Total*100), digits=3))


bulk<-ggplot(data=freq, aes(x=reorder(Cell_Type, -percent), y=percent, color=Cell_Type)) +
  geom_text(aes(label=percent), vjust=-0.3, size=3.5)+
  geom_bar(stat="identity", fill="white", lwd=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="Cell Type") +
  labs(y="% Ki67+ of Total") +
  guides(fill=FALSE) 
bulk

#Across samples
freq_sample<-as.data.frame(table(data_lymph$cell_type, data_lymph$SampleID))
freq_Ki67_sample<-as.data.frame(table(data_lymph_Ki67$cell_type,data_lymph_Ki67$SampleID))
freq_sample$Ki67pos<-0
freq_sample[freq_sample$Var2 %in% freq_Ki67_sample$Var2 & freq_sample$Var1 %in% freq_Ki67_sample$Var1,]$Ki67pos<-freq_Ki67_sample$Freq
names(freq_sample)<-c("Cell_Type","SampleID","Total","Total_Ki67pos")
freq_sample$percent<-as.numeric(format((freq_sample$Total_Ki67pos / freq_sample$Total)*100),digits=3)

##..Plot..##

my_comparisons <- list(c("Treg","CD4_T"),c("Treg","CD8_T"),c("Treg","B_cell"))

per_sample<-ggplot(na.omit(freq_sample), aes(x=fct_reorder(Cell_Type,percent,.fun=median,.desc=TRUE), y=percent, fill=Cell_Type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  stat_compare_means(data=na.omit(freq_sample),mapping=aes(x=fct_reorder(Cell_Type,percent,.fun=median,.desc=TRUE),
  y=percent) ,label = "p.signif", method= "wilcox.test",
  comparisons=my_comparisons) +
  theme(legend.position = "none") +
  labs(x="Cell Type") +
  labs(y="% Ki67+ of Total") +
  guides(fill=guide_legend(title="Cell Type"))
per_sample

compare_means(percent ~ Cell_Type, freq_sample, method= "wilcox.test")

wilcox.exact(percent ~ Cell_Type, freq_sample, Cell_Type %in% c('Treg','B_cell'))

