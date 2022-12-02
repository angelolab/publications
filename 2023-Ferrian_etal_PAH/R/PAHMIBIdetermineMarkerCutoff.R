# PAHMIBIdetermineMarkerCutoff.R
# Author: Erin McCaffrey
# Date created: 200327
# Overview: Script reads in the normalized expression matrix for all cohort samples (PAH and Healthy)
# It plots a histogram for every single marker in order to determine a cutoff value for positivity on a 
# per marker basis for downstream analysis. Also implements silhouette-scanning approach from Hu et al 2018 
# Cell Reports to suggest a cutoff for markers from a 1D distribution

##..Libraries needed..##

library(Hmisc)
library("MetaCyto")
library(tidyr)
library(plyr)
library(ggplot2)

##..Import data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")

##..Get silhouette scanning threshold of all markers..##

channelCols <- c(4:39) #channels to assess
channelNames <- colnames(data[,channelCols]) #channel names

#create dataframe to store thresholds for each marker
cutoffs <- as.data.frame(colnames(data[,channelCols])) 
cutoffs$t <- 0
colnames(cutoffs) <- c("marker","Threshold")

for (i in 1:length(channelNames)) {
  channel <- channelNames[i]
  channelData <- data[,channel]
  cutoff <- findCutoff(channelData,useBL=FALSE) #find cutoff function from metacyto
  cutoffs[cutoffs$marker==channelNames[i],]$Threshold <- cutoff
}

##..Manually adjust several based on histogram inspection..##

cutoffs[cutoffs$marker == "GranzB",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "HLA.Class.I",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "HLA.DR",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "IDO1",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "IFNg",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "iNOS",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "iNOS",]$Threshold <- 0.15
cutoffs[cutoffs$marker == "NaK.ATPase",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "PanCK",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "PDL1b",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "SMA",]$Threshold <- 0.25
cutoffs[cutoffs$marker == "TIM.3",]$Threshold <- 0.25


##..Plot histograms of all markers..##
data_melted <- data[,channelCols] %>% gather(colnames(data[,channelCols]),key="marker",value="value")
ggplot(data_melted[which(data_melted$value>0),],aes(x=value)) + 
  geom_histogram(aes(x = value, y = ..density..), bins=100) + 
  geom_density(aes(x = value), color = "blue") +
  geom_vline(data = cutoffs, aes(xintercept=Threshold, color ="red")) +
  facet_wrap(~marker,scales=c("free")) 

##..Save the thresholds as a csv for additional analysis..##

write.csv(cutoffs, file="markerThresholds.csv",row.names = FALSE)
