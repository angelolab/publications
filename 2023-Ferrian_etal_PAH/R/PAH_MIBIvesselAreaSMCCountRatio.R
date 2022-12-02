# PAH_MIBIvesselAreaSMCCountRatio.R
# Author: Erin McCaffrey 
# Date created: 200908

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data<-read.csv("vessel_area_data.csv")

##..Sum per point..##
data_summary <- data %>%
  group_by(Point_num) %>%
  summarize(total_vessel_area_px = sum(area))

##..Append smooth muscle cell count..##
annotated_data<-read.csv("celldata_region_annotated.csv")
SMC_data<-droplevels(annotated_data[annotated_data$cell_lineage=='fibroblast',])
SMC_counts<-as.data.frame(table(SMC_data$Point_num,SMC_data$cell_lineage))

#drop points withut vessel data
SMC_counts_vessel<-droplevels(SMC_counts[SMC_counts$Var1 %in% data_summary$Point_num,])

#add SMC counts to summary of vessel area
data_summary$SMC_counts<-SMC_counts_vessel$Freq

##..Calculate vessel area:SMC count ratio..##
data_summary$vessel_area_SMC_count_ratio<-data_summary$total_vessel_area_px/data_summary$SMC_counts

##..Get log2 of ratio..##
data_summary$log2_vessel_area_SMC_count_ratio<-log2(data_summary$vessel_area_SMC_count_ratio)

##..Append condition info..##
ipah<-unique(annotated_data[annotated_data$Subgroup=='IPAH',]$Point_num)
hpah<-unique(annotated_data[annotated_data$Subgroup=='HPAH',]$Point_num)
hlt<-unique(annotated_data[annotated_data$Subgroup=='Healthy Controls',]$Point_num)
data_summary$group<-'hlt'
data_summary[data_summary$Point_num %in% hpah, ]$group <- 'hpah'
data_summary[data_summary$Point_num %in% ipah, ]$group <- 'ipah'

##..Export..##
write.csv(data_summary, 'vesselarea_SMCcount_ratio.csv',row.names = F)




