# PAH_MIBIMPOinNeutroPerNC.R
# Author: Erin McCaffrey
# Date created: 200909

library(ggpubr)
library(ggplot2)
library(colorspace)
library(forcats)
library(dplyr)
library(tibble)
library(reshape2)
library(gplots)
library(RColorBrewer)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")
data_neighborhood<-read.csv("PAH_AllCells_Neighborhood_K=10_PatientAnnotated.csv")

##..Drop healthy..##

data_pah<-droplevels(data[data$Subgroup %in% c('HPAH','IPAH'),])

##..Merge the neighborhood data with the expression data..##

NC_data<-merge(data_neighborhood, data_pah, by=c("Point_num","label",'PID'))

##...Filter to only include Neutro cells...##

data_Neutro<-droplevels(NC_data[NC_data$cell_lineage =='Neutro',])

##...Get percent of cells positive for TIM3 and ratio of pos:neg...##

marker_thresh <- 0.46
data_marker<-droplevels(data_Neutro[data_Neutro$MPO>marker_thresh, ])

#Across PID
freq_sample<-as.data.frame(table(data_Neutro$PID, data_Neutro$cluster))
freq_marker_sample<-as.data.frame(table(data_marker$PID, data_marker$cluster))

freq_sample$MPOpos<-0
freq_sample[freq_sample$Var1 %in% freq_marker_sample$Var1,]$MPOpos<-freq_marker_sample$Freq
names(freq_sample)<-c("PID","cluster","Total_Neutro","Total_MPOpos")
freq_sample$percentMPOpos<-as.numeric(format((freq_sample$Total_MPOpos / freq_sample$Total_Neutro)*100),digits=3)

write.csv(freq_sample, "MPO_Neutro_NC_PerPatient.csv", row.names = F)

#Across Point_num
freq_sample<-as.data.frame(table(data_Neutro$Point_num, data_Neutro$cluster))
freq_marker_sample<-as.data.frame(table(data_marker$Point_num, data_marker$cluster))

freq_sample$MPOpos<-0
freq_sample[freq_sample$Var1 %in% freq_marker_sample$Var1,]$MPOpos<-freq_marker_sample$Freq
names(freq_sample)<-c("Point_num","cluster","Total_Neutro","Total_MPOpos")
freq_sample$percentMPOpos<-as.numeric(format((freq_sample$Total_MPOpos / freq_sample$Total_Neutro)*100),digits=3)

write.csv(freq_sample, "MPO_Neutro_NC_PerPoint.csv", row.names = F)


# Mean per PID per NC

data_summary <- data_Neutro %>%
  group_by(PID, cluster) %>%
  summarize(mean = mean(MPO), 
            median = median(MPO))

write.csv(data_summary, 'MPO_mean-median_perNC_PerPatient.csv', row.names = F)

# Mean per Point_num per NC

data_summary <- data_Neutro %>%
  group_by(Point_num, cluster) %>%
  summarize(mean = mean(MPO), 
            median = median(MPO))

write.csv(data_summary, 'MPO_mean-median_perNC_PerPoint.csv', row.names = F)



