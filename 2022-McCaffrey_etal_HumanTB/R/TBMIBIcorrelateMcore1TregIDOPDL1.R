# TBMIBIcorrelateMcore1TregIDOCD8.R
# Author: Erin McCaffrey
# Date created: 210803
# Script reads in the cellular data for MEs. For ME Mcore1, it enumerates the frequency of
# IDO1+ cells in Mcore1, the frequency of Tregs in Mcore1, and IDO1+ cells in Mcore1. Next
# it assesses the linear relationship between these features.

library(ggplot2)
library(dplyr)

##...Load in data..##

data<-read.csv("data/allTB-ME_annotated.csv")
mcore1_data<-data[data$maxME==1,]

##..Create results dataframe to store all data..##

freq_sample<-as.data.frame(table(mcore1_data$SampleID))
names(freq_sample)<-c("SampleID","Total_Mcore1")

##..Get # IDO1+ cells all MEs per sample..##

IDO_thresh <- 0.26
data_marker_IDO<-data[data$IDO>IDO_thresh, ]
freq_IDO_sample<-as.data.frame(table(data_marker_IDO$SampleID))
freq_sample$Total_IDO1pos<-0
freq_sample[freq_sample$SampleID %in% freq_IDO_sample$Var1,]$Total_IDO1pos <- freq_IDO_sample[freq_IDO_sample$Var1 %in% freq_sample$SampleID, ]$Freq

##..Get # IDO1+ cells in Mcore1 per sample..##

data_marker_IDO_mcore1<-mcore1_data[mcore1_data$IDO>IDO_thresh, ]
freq_IDO_sample_mcore1<-as.data.frame(table(data_marker_IDO_mcore1$SampleID))
freq_sample$Total_IDO1pos_Mcore1<-0
freq_sample[freq_sample$SampleID %in% freq_IDO_sample_mcore1$Var1,]$Total_IDO1pos_Mcore1 <- freq_IDO_sample_mcore1$Freq
freq_sample$percentIDO1pos_Mcore1<-as.numeric(format((freq_sample$Total_IDO1pos_Mcore1 / freq_sample$Total_Mcore1)*100),digits=3)


##..Get # PDL1+ cells all MEs per sample..##

PDL1_thresh <- 0.25
data_marker_PDL1<-data[data$PD.L1>PDL1_thresh, ]
freq_PDL1_sample<-as.data.frame(table(data_marker_PDL1$SampleID))
freq_sample$Total_PDL1pos<-0
freq_sample[freq_sample$SampleID %in% freq_PDL1_sample$Var1,]$Total_PDL1pos <- freq_PDL1_sample[freq_PDL1_sample$Var1 %in% freq_sample$SampleID, ]$Freq

##..Get # PDL1+ cells in Mcore1 per sample..##

data_marker_PDL1_mcore1<-mcore1_data[mcore1_data$PD.L1>PDL1_thresh, ]
freq_PDL1_sample_mcore1<-as.data.frame(table(data_marker_PDL1_mcore1$SampleID))
freq_sample$Total_PDL1pos_Mcore1<-0
freq_sample[freq_sample$SampleID %in% freq_PDL1_sample_mcore1$Var1,]$Total_PDL1pos_Mcore1 <- freq_PDL1_sample_mcore1$Freq
freq_sample$percentPDL1pos_Mcore1<-as.numeric(format((freq_sample$Total_PDL1pos_Mcore1 / freq_sample$Total_Mcore1)*100),digits=3)

##..Get frequency and count of Tregs in Mcore1 per sample..##

freq_treg<-as.data.frame(table(data[data$cell_type=='Treg',]$SampleID))
freq_sample$n_Treg<-0
freq_sample[freq_sample$SampleID %in% freq_treg$Var1,]$n_Treg <- freq_treg[freq_treg$Var1 %in% freq_sample$SampleID,]$Freq

freq_treg_mcore1<-as.data.frame(table(mcore1_data[mcore1_data$cell_type=='Treg',]$SampleID))
freq_sample$n_Treg_Mcore1<-0
freq_sample[freq_sample$SampleID %in% freq_treg_mcore1$Var1,]$n_Treg_Mcore1 <- freq_treg_mcore1[freq_treg_mcore1$Var1 %in% freq_sample$SampleID,]$Freq

freq_sample$percentTreg_of_Mcore1<-as.numeric(format((freq_sample$n_Treg_Mcore1 / freq_sample$Total_Mcore1)*100),digits=3)
freq_sample$percentTreg_in_Mcore1<-as.numeric(format((freq_sample$n_Treg_Mcore1 / freq_sample$n_Treg)*100),digits=3)

##..Plot individual relationships..##

# Total_IDO1pos_Mcore1 with n_Treg_Mcore1
summary(lm(n_Treg_Mcore1~Total_IDO1pos_Mcore1, freq_sample))
corr<-cor.test(freq_sample$n_Treg_Mcore1,freq_sample$Total_IDO1pos_Mcore1, method="pearson")

ggplot(freq_sample,aes(x=	Total_IDO1pos_Mcore1,y=n_Treg_Mcore1)) +
  geom_smooth(method='lm', formula= y~x) +
  geom_point(aes(group=1, size = 1)) + 
  labs(x="# IDO1+ Cells in Mcore1") + 
  labs(y="# Tregs in Mcore1") + theme_bw() + 
  theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

# Total_PDL1pos_Mcore1 with	n_Treg_Mcore1
summary(lm(n_Treg_Mcore1~Total_PDL1pos_Mcore1, freq_sample))
corr<-cor.test(freq_sample$Total_PDL1pos_Mcore1,freq_sample$n_Treg_Mcore1, method="pearson")

ggplot(freq_sample,aes(x=	Total_PDL1pos_Mcore1,y=n_Treg_Mcore1)) +
  geom_smooth(method='lm', formula= y~x) +
  geom_point(aes(group=1, size = 1)) + 
  labs(x="# PDL1+ Cells in Mcore1") + 
  labs(y="# Tregs in Mcore1") + theme_bw() + 
  theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


