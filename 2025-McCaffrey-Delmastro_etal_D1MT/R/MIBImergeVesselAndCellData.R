# MIBImergeVesselAndCellData.R
# Author: Erin McCaffrey 
# Date created: 200318
# Overview: This script reads in the csv for all NHP cell data and the core annotation data. 
# It produces a merged dataframe, assigns each cell to a mask category and exports as a csv. 

library(dplyr)


##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv('NHP_cohort_data_norm_annotated.csv')
data_mask<-read.csv("core-mask-data.csv")

##..Merge the dataframes..##
merged<-data %>% left_join(data_mask, by = c('SampleID','cellLabelInImage'))

##..Add mock 0s to point 22..##

merged[merged$SampleID==22, ]$In_Mask<-0
merged$region<-NA
merged[merged$In_Mask==0,]$region<-'peripheral'
merged[merged$In_Mask==1,]$region<-'core'

write.csv(merged, file="celldata_region_annotated.csv",row.names = FALSE)
