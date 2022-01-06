# TBMIBIdetermineMarkerCutoff.R
# Author: Erin McCaffrey
# Date created: 190728
# Overview: Script reads in the normalized expression matrix for all cohort samples (myco and sarcoid). 
# It plots a histogram for every single marker in order to determine a cutoff value for positivity on a 
# per marker basis for downstream analysis. Also implements silhouette-scanning approach from Hu et al 2018 
# Cell Reports to suggest a cutoff for markers from a 1D distribution

##..Libraries needed..##

library(Hmisc)
library("MetaCyto")
library(tidyr)
library(plyr)
library(ggplot2)

##..Set directory and read in data..##

setwd("~/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

data<-read.csv("granA_cellpheno_CS-asinh-norm.csv") #cell size normalized, scaled by 100, asinh trans, 99th percentile scaled

##..Get silhouette scanning threshold of all markers..##

channelCols <- c(7:9,11:44) #channels to assess
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

##..Plot histograms of all markers..##
data_melted <- data[,channelCols] %>% gather(colnames(data[,channelCols]),key="marker",value="value")
ggplot(data_melted[which(data_melted$value>0),],aes(x=value)) + 
  geom_histogram(aes(x = value, y = ..density..), bins=100) + 
  geom_density(aes(x = value), color = "blue") +
  geom_vline(data = cutoffs, aes(xintercept=Threshold, color ="red")) +
  facet_wrap(~marker,scales=c("free")) 

##..Save the thresholds as a csv for additional analysis..##

write.csv(cutoffs, file="markerThresholds.csv",row.names = FALSE)

