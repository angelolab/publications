# PAH_MIBIDCFunctionalMarkersPerCellSubtype_Neg_SinglePos.R
# Author: Erin McCaffrey 

library(dplyr)
library(ggplot2)
library(ggpubr)

##..Import data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data_DC<-read.csv("celldata_region_annotated.csv")

##..Define thresholds for TIM3 and IDO1 and SAMHD1..##

IDO_thresh<-0.25
TIM3_thresh<-0.25
SAMHD1_thresh<-0.39

##..Filter to include only TIM3 or IDO1+ DCs or TIM3+ IDO1+ SAMHD1+ DCs..##

data_DC_IDO1<-droplevels(data_DC[data_DC$IDO1>IDO_thresh & data_DC$TIM.3<=TIM3_thresh & data_DC$SAMHD1<=SAMHD1_thresh, ])
data_DC_TIM3<-droplevels(data_DC[data_DC$IDO1<=IDO_thresh & data_DC$TIM.3>TIM3_thresh & data_DC$SAMHD1<=SAMHD1_thresh, ])
data_DC_SAMHD1<-droplevels(data_DC[data_DC$IDO1<=IDO_thresh & data_DC$TIM.3<=TIM3_thresh & data_DC$SAMHD1>SAMHD1_thresh, ])

data_DC_neg<-droplevels(data_DC[data_DC$IDO1<=IDO_thresh & data_DC$TIM.3<=TIM3_thresh & data_DC$SAMHD1<=SAMHD1_thresh, ])
data_DC_pos<-droplevels(data_DC[data_DC$IDO1>IDO_thresh & data_DC$TIM.3>TIM3_thresh & data_DC$SAMHD1>SAMHD1_thresh, ])

data_DC_IDO1_TIM3<-droplevels(data_DC[data_DC$IDO1>IDO_thresh & data_DC$TIM.3>TIM3_thresh & data_DC$SAMHD1<=SAMHD1_thresh, ])
data_DC_IDO1_SAMHD1<-droplevels(data_DC[data_DC$IDO1>IDO_thresh & data_DC$TIM.3<=TIM3_thresh & data_DC$SAMHD1>SAMHD1_thresh, ])
data_DC_TIM3_SAMHD1<-droplevels(data_DC[data_DC$IDO1<=IDO_thresh & data_DC$TIM.3>TIM3_thresh & data_DC$SAMHD1>SAMHD1_thresh, ])


##..Create a table with total DC counts and DC counts for each functional marker category and subtype per patient..##

cell_counts<-as.data.frame(table(data_DC$PID, data_DC$cell_lineage))
colnames(cell_counts)<-c("PID","cell_lineage","cell_count")

#single pos
DC_IDO_counts<-as.data.frame(table(data_DC_IDO1$PID, data_DC_IDO1$cell_lineage))
colnames(DC_IDO_counts)<-c("PID","cell_lineage","cell_count")
DC_TIM3_counts<-as.data.frame(table(data_DC_TIM3$PID, data_DC_TIM3$cell_lineage))
colnames(DC_TIM3_counts)<-c("PID","cell_lineage","cell_count")
DC_SAMHD1_counts<-as.data.frame(table(data_DC_SAMHD1$PID, data_DC_SAMHD1$cell_lineage))
colnames(DC_SAMHD1_counts)<-c("PID","cell_lineage","cell_count")

#triple neg/pos
DC_neg_counts<-as.data.frame(table(data_DC_neg$PID, data_DC_neg$cell_lineage))
colnames(DC_neg_counts)<-c("PID","cell_lineage","cell_count")
DC_pos_counts<-as.data.frame(table(data_DC_pos$PID, data_DC_pos$cell_lineage))
colnames(DC_pos_counts)<-c("PID","cell_lineage","cell_count")

#double pos
DC_IDO1_TIM3_counts<-as.data.frame(table(data_DC_IDO1_TIM3$PID, data_DC_IDO1_TIM3$cell_lineage))
colnames(DC_IDO1_TIM3_counts)<-c("PID","cell_lineage","cell_count")
DC_IDO1_SAMHD1_counts<-as.data.frame(table(data_DC_IDO1_SAMHD1$PID, data_DC_IDO1_SAMHD1$cell_lineage))
colnames(DC_IDO1_SAMHD1_counts)<-c("PID","cell_lineage","cell_count")
DC_TIM3_SAMHD1_counts<-as.data.frame(table(data_DC_TIM3_SAMHD1$PID, data_DC_TIM3_SAMHD1$cell_lineage))
colnames(DC_TIM3_SAMHD1_counts)<-c("PID","cell_lineage","cell_count")


##..Determine frequency for each DC type of total..##

DC_freqs<-cell_counts

#IDO1
DC_freqs$count_IDO1_singlepos<-0
DC_freqs[DC_freqs$PID %in% DC_IDO_counts$PID &
           DC_freqs$cell_lineage %in% DC_IDO_counts$cell_lineage,]$count_IDO1_singlepos<-DC_IDO_counts$cell_count
DC_freqs$freq_IDO1_singlepos<-as.numeric(format((DC_freqs$count_IDO1_singlepos / DC_freqs$cell_count)),digits=3)
#TIM3
DC_freqs$count_TIM3_singlepos<-0
DC_freqs[DC_freqs$PID %in% DC_TIM3_counts$PID & 
           DC_freqs$cell_lineage %in% DC_TIM3_counts$cell_lineage,]$count_TIM3_singlepos<-DC_TIM3_counts$cell_count
DC_freqs$freq_TIM3_singlepos<-as.numeric(format((DC_freqs$count_TIM3_singlepos / DC_freqs$cell_count)),digits=3)
#SAMHD1
DC_freqs$count_SAMHD1_singlepos<-0
DC_freqs[DC_freqs$PID %in% DC_SAMHD1_counts$PID &
           DC_freqs$cell_lineage %in% DC_SAMHD1_counts$cell_lineage,]$count_SAMHD1_singlepos<-DC_SAMHD1_counts$cell_count
DC_freqs$freq_SAMHD1_singlepos<-as.numeric(format((DC_freqs$count_SAMHD1_singlepos / DC_freqs$cell_count)),digits=3)
#neg
DC_freqs$triple_neg_count<-0
DC_freqs[DC_freqs$PID %in% DC_neg_counts$PID & 
           DC_freqs$cell_lineage %in% DC_neg_counts$cell_lineage,]$triple_neg_count<-DC_neg_counts$cell_count
DC_freqs$freq_triple_neg<-as.numeric(format((DC_freqs$triple_neg_count / DC_freqs$cell_count)),digits=3)
#pos
DC_freqs$triple_pos_count<-0
DC_freqs[DC_freqs$PID %in% DC_pos_counts$PID & 
           DC_freqs$cell_lineage %in% DC_pos_counts$cell_lineage,]$triple_pos_count<-DC_pos_counts$cell_count
DC_freqs$freq_triple_pos<-as.numeric(format((DC_freqs$triple_pos_count / DC_freqs$cell_count)),digits=3)
#IDO1 TIM3
DC_freqs$IDO1_TIM3_count<-0
DC_freqs[DC_freqs$PID %in% DC_IDO1_TIM3_counts$PID & 
           DC_freqs$cell_lineage %in% DC_IDO1_TIM3_counts$cell_lineage,]$IDO1_TIM3_count<-DC_IDO1_TIM3_counts$cell_count
DC_freqs$freq_IDO1_TIM3<-as.numeric(format((DC_freqs$IDO1_TIM3_count / DC_freqs$cell_count)),digits=3)
#IDO1 SAMHD1
DC_freqs$IDO1_SAMHD1_count<-0
DC_freqs[DC_freqs$PID %in% DC_IDO1_SAMHD1_counts$PID & 
           DC_freqs$cell_lineage %in% DC_IDO1_SAMHD1_counts$cell_lineage,]$IDO1_SAMHD1_count<-DC_IDO1_SAMHD1_counts$cell_count
DC_freqs$freq_IDO1_SAMHD1<-as.numeric(format((DC_freqs$IDO1_SAMHD1_count / DC_freqs$cell_count)),digits=3)
#TIM3 SAMHD1
DC_freqs$TIM3_SAMHD1_count<-0
DC_freqs[DC_freqs$PID %in% DC_TIM3_SAMHD1_counts$PID & 
           DC_freqs$cell_lineage %in% DC_TIM3_SAMHD1_counts$cell_lineage,]$TIM3_SAMHD1_count<-DC_TIM3_SAMHD1_counts$cell_count
DC_freqs$freq_TIM3_SAMHD1<-as.numeric(format((DC_freqs$TIM3_SAMHD1_count / DC_freqs$cell_count)),digits=3)


##..Export..##
write.csv(DC_freqs, "cell_lineage_FunctionalMarkers_PerPID_UniqueCombos.csv", row.names = F)


