# TBMIBIpercentIDOPDL1pos.R
# Author: Erin McCaffrey 
# Date created: 190626
# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# IDO and PDL1 pos cells for each sample, split by sarcoid and mycobacteria and exports csvs for
# later feature comparison.

library(ggplot2)
library(dplyr)
library(ggpubr)

##..Import data..##

data_norm<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Get percent IDO and PDL1 positive.##

IDO_thresh = 0.26
PDL1_thresh = 0.25

#IDO
data_IDO<-data.frame(table(data_norm$SampleID))
IDOpos<-data_norm[data_norm$IDO >= IDO_thresh,]
data_IDOpos<-data.frame(table(IDOpos$SampleID))

data_IDO$IDOpos<-0
data_IDO[data_IDO$Var1 %in% data_IDOpos$Var1,]$IDOpos<-data_IDOpos$Freq
names(data_IDO)<-c("SampleID","Total","Total_IDOpos")
data_IDO$percentIDO<-as.numeric(format((data_IDO$Total_IDOpos / data_IDO$Total)*100),digits=3)

#PDL1
data_PD.L1<-data.frame(table(data_norm$SampleID))
PD.L1pos<-data_norm[data_norm$PD.L1 >= PDL1_thresh,]
data_PD.L1pos<-data.frame(table(PD.L1pos$SampleID))

data_PD.L1$PD.L1pos<-0
data_PD.L1[data_PD.L1$Var1 %in% data_PD.L1pos$Var1,]$PD.L1pos<-data_PD.L1pos$Freq
names(data_PD.L1)<-c("SampleID","Total","Total_PD.L1pos")
data_PD.L1$percentPDL1<-as.numeric(format((data_PD.L1$Total_PD.L1pos / data_PD.L1$Total)*100),digits=3)

# Add tissue annotation

sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
data_PD.L1$Tissue<-'mycobacteria'
data_PD.L1[data_PD.L1$SampleID %in% sarcoid, ]$Tissue<-'sarcoid'
data_IDO$Tissue<-'mycobacteria'
data_IDO[data_IDO$SampleID %in% sarcoid, ]$Tissue<-'sarcoid'


# Plot sarcoid v TB PD-L1 and IDO1

ggplot(data = data_PD.L1, aes(x = Tissue, y = percentPDL1, fill = Tissue)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = 3, position = position_jitterdodge()) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency')

ggplot(data = data_IDO, aes(x = Tissue, y = percentIDO, fill = Tissue)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = 3, position = position_jitterdodge()) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency')
