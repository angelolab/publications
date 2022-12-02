# PAH_MIBIPercentPositiveAllMarkersAllSubsets.R
# Author: Erin McCaffrey 
# Date created: 200310


require(dplyr)


##..Import data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")

##..Create a dataframe with all points and all cells counts..##

cell_counts<-as.data.frame(table(data$PID, data$cell_lineage))
colnames(cell_counts)<-c("PID","cell_lineage","cell_count")
cell_freqs<-cell_counts

##..Read in markers thresholds..##

marker_thresholds<-read.csv("markerThresholds.csv")
markers<-as.vector(marker_thresholds$marker)
thresholds<-marker_thresholds$Threshold

##..Append frequency of positive cells for all patients and subsets..##

for(i in 1:length(markers)) {
  #define marker and threshold
  marker<-markers[i]
  threshold<-thresholds[i]
  #subset positive cells
  data_marker<-data  %>% filter(get(marker)>=threshold)
  #get table of positive cell counts per patient per cell type
  freq_marker_patient<-as.data.frame(table(data_marker$PID,data_marker$cell_lineage))
  #append positive counts and percent positive to summary dataframes
  cell_counts[[marker]]<-0
  cell_counts[cell_counts$PID %in% freq_marker_patient$Var1 & cell_counts$cell_lineage
              %in% freq_marker_patient$Var2,][[marker]]<-freq_marker_patient$Freq
  cell_freqs[[marker]]<-as.numeric(format((cell_counts[[marker]] / cell_counts$cell_count)),digits=3)
}

##..Add group and subgroup annotation..##

PAH<-unique(data[data$Group=="PAH",]$PID)
Hlt<-unique(data[data$Group=="Healthy Controls",]$PID)

cell_counts$Group<-'hlt'
cell_counts[cell_counts$PID %in% PAH,]$Group<-'pah'
cell_freqs$Group<-cell_counts$Group


IPAH<-unique(data[data$Subgroup=="IPAH",]$PID)
HPAH<-unique(data[data$Subgroup=="HPAH",]$PID)

cell_counts$Subgroup<-'hlt'
cell_counts[cell_counts$PID %in% IPAH,]$Subgroup<-'ipah'
cell_counts[cell_counts$PID %in% HPAH,]$Subgroup<-'hpah'
cell_freqs$Subgroup<-cell_counts$Subgroup

##..Export csv..##

write.csv(cell_counts,"allpatients_count-of-positivecells-persubset.csv",row.names = FALSE)
write.csv(cell_freqs,"allpatients_freq-of-positivecells-persubset.csv",row.names = FALSE)
