# PAH_MIBIDCFunctionalMarkersPerDCSubtype.R
# Author: Erin McCaffrey 

library(dplyr)
library(ggplot2)
library(ggpubr)

##..Import data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data_DC<-read.csv("DCMetaclust_K=20.csv")

##..Define thresholds for TIM3 and IDO1 and SAMHD1..##

IDO_thresh<-0.25
TIM3_thresh<-0.25
SAMHD1_thresh<-0.39

##..Filter to include only TIM3 or IDO1+ DCs or TIM3+ IDO1+ SAMHD1+ DCs..##

data_DC_IDO1_TIM3<-droplevels(data_DC[data_DC$IDO1>IDO_thresh & data_DC$TIM.3>TIM3_thresh, ])
data_DC_IDO1_TIM3_SAMHD1<-droplevels(data_DC[data_DC$IDO1>IDO_thresh & data_DC$TIM.3>TIM3_thresh
                                             & data_DC$SAMHD1>SAMHD1_thresh, ])

##..Create a table with total DC counts and DC counts for each functional marker category and subtype per patient..##

DC_counts<-as.data.frame(table(data_DC$PID, data_DC$cell_lineage2))
colnames(DC_counts)<-c("PID","DC_Subtype","DC_count")

DC_IDO_TIM3_counts<-as.data.frame(table(data_DC_IDO1_TIM3$PID, data_DC_IDO1_TIM3$cell_lineage2))
colnames(DC_IDO_TIM3_counts)<-c("PID","DC_Subtype","DC_count")

DC_IDO_TIM3_SAMHD1_counts<-as.data.frame(table(data_DC_IDO1_TIM3_SAMHD1$PID, data_DC_IDO1_TIM3_SAMHD1$cell_lineage2))
colnames(DC_IDO_TIM3_SAMHD1_counts)<-c("PID","DC_Subtype","DC_count")

##..Determine frequency for each DC type of total..##

DC_freqs<-DC_counts

#TIM3 and IDO1
DC_freqs$count_TIM3_IDO1<-0
DC_freqs[DC_freqs$PID %in% DC_IDO_TIM3_counts$PID,]$count_TIM3_IDO1<-DC_IDO_TIM3_counts$DC_count
DC_freqs$freq_TIM3_IDO1<-as.numeric(format((DC_freqs$count_TIM3_IDO1 / DC_freqs$DC_count)),digits=3)

#TIM3 and IDO1 and SAMHD1
DC_freqs$count_IDO_TIM3_SAMHD1<-0
DC_freqs[DC_freqs$PID %in% DC_IDO_TIM3_SAMHD1_counts$PID,]$count_IDO_TIM3_SAMHD1<-DC_IDO_TIM3_SAMHD1_counts$DC_count
DC_freqs$freq_IDO_TIM3_SAMHD1<-as.numeric(format((DC_freqs$count_IDO_TIM3_SAMHD1 / DC_freqs$DC_count)),digits=3)

##..Export..##
write.csv(DC_freqs, "DC_Subtype_FunctionalMarkers_PerPID.csv", row.names = F)


##..Create a table with total DC counts and DC counts for each functional marker category and subtype per point..##

DC_counts<-as.data.frame(table(data_DC$Point_num, data_DC$cell_lineage2))
colnames(DC_counts)<-c("Point_num","DC_Subtype","DC_count")

DC_IDO_TIM3_counts<-as.data.frame(table(data_DC_IDO1_TIM3$Point_num, data_DC_IDO1_TIM3$cell_lineage2))
colnames(DC_IDO_TIM3_counts)<-c("Point_num","DC_Subtype","DC_count")

DC_IDO_TIM3_SAMHD1_counts<-as.data.frame(table(data_DC_IDO1_TIM3_SAMHD1$Point_num, data_DC_IDO1_TIM3_SAMHD1$cell_lineage2))
colnames(DC_IDO_TIM3_SAMHD1_counts)<-c("Point_num","DC_Subtype","DC_count")

##..Determine frequency for each DC type of total..##

DC_freqs<-DC_counts

#TIM3 and IDO1
DC_freqs$count_TIM3_IDO1<-0
DC_freqs[DC_freqs$Point_num %in% DC_IDO_TIM3_counts$Point_num,]$count_TIM3_IDO1<-DC_IDO_TIM3_counts$DC_count
DC_freqs$freq_TIM3_IDO1<-as.numeric(format((DC_freqs$count_TIM3_IDO1 / DC_freqs$DC_count)),digits=3)

#TIM3 and IDO1 and SAMHD1
DC_freqs$count_IDO_TIM3_SAMHD1<-0
DC_freqs[DC_freqs$Point_num %in% DC_IDO_TIM3_SAMHD1_counts$Point_num,]$count_IDO_TIM3_SAMHD1<-DC_IDO_TIM3_SAMHD1_counts$DC_count
DC_freqs$freq_IDO_TIM3_SAMHD1<-as.numeric(format((DC_freqs$count_IDO_TIM3_SAMHD1 / DC_freqs$DC_count)),digits=3)

##..Export..##
write.csv(DC_freqs, "DC_Subtype_FunctionalMarkers_PerPoint.csv", row.names = F)

