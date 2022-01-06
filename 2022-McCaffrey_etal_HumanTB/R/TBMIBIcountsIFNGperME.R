# TBMIBIcountsIFNGiNOSPositiveCellsperME.R
# Author: Erin McCaffrey
# Date created: 201117
# Overview: This script reads in the csv for the single cell data and counts the number of IFNG+ cells FOV and
# per ME and plots it as a barplot.

library(ggplot2)
library(forcats)
library(dplyr)

##...Load in data..##

data<-read.csv('data/allTB-ME_annotated.csv')

##..Counts of all IFNG positive cells across MEs..##

IFNG_pos<-as.data.frame(table(data[data$IFNg>0.43, ]$SampleID, data[data$IFNg>0.43, ]$maxME))
colnames(IFNG_pos)<-c('SampleID','ME','Count')

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

ggplot(data = IFNG_pos, aes(x = as.factor(SampleID), y = Count, fill = ME)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = colors) +
  theme_bw() + 
  labs(x = 'ME') + 
  labs(y = 'Count IFNG+ Cells')

##..Counts of all IFNG positive cells across just ME1, include all ROIs..##

ME1_data<-data[data$maxME == 1, ]
IFNG_pos_ME1<-as.data.frame(table(ME1_data[ME1_data$IFNg>0.43, ]$SampleID, ME1_data[ME1_data$IFNg>0.43, ]$maxME))
colnames(IFNG_pos_ME1)<-c('SampleID','ME','Count')

# add zeros for ROIs with none
SampleID<-setdiff(unique(data$SampleID), IFNG_pos_ME1$SampleID)
ME<-rep(1,length(SampleID))
Count<-rep(0,length(SampleID))

neg_ROIs<-data.frame(SampleID, ME, Count)
neg_ROIs$SampleID<-as.factor(neg_ROIs$SampleID)

IFNG_pos_ME1_all<-rbind(IFNG_pos_ME1, neg_ROIs)

IFNG_pos_ME1_all$SampleID<-as.numeric(as.character(IFNG_pos_ME1_all$SampleID))
IFNG_pos_ME1_all<-arrange(IFNG_pos_ME1_all,SampleID)

ggplot(data = IFNG_pos_ME1_all, aes(x = reorder(SampleID, -Count), y = Count, fill = ME)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#CA6627')) +
  theme_bw() + 
  labs(x = 'Patient-ROI') + 
  labs(y = 'Count IFNG+ Cells')

