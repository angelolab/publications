# TBMIBIgiantcellIDOPDL1.R
# Author: Erin McCaffrey 
# Date created: 190925
# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# IDO and PDL1 pos cells in giant cells. It also plots the count of GCs per fov.

library(ggplot2)

##..Import data..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Subset just TB samples..##

data_norm<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])

##..Get percent IDO and PDL1 positive..##

# modified threshold due to large cell size
IDO_thresh = 0
PDL1_thresh = 0

##..Get total number of IDO1 and PDL1+ giant cells..##

GC<-as.data.frame(as.numeric(sum(data_norm$cell_type=="giant_cell")))
colnames(GC)<-"GC_total"
GC$IDOpos<-as.numeric(sum(data_norm[data_norm$cell_type=="giant_cell",]$IDO>IDO_thresh))
GC$PDL1pos<-as.numeric(sum(data_norm[data_norm$cell_type=="giant_cell",]$PD.L1>PDL1_thresh))

##..Get percent IDO1 and PDL1 positive..##

GC$IDOpos_percent<-GC$IDOpos/GC$GC_total
GC$PDL1pos_percent<-GC$PDL1pos/GC$GC_total

markers<-c("IDO","PDL1")
percent<-c(GC$IDOpos_percent,GC$PDL1pos_percent)
GC_markers<-data.frame(markers,percent)

##..Plot..##

bulk<-ggplot(data=GC_markers, aes(x=markers,y=percent)) +
  geom_bar(stat="identity", fill="grey", lwd=0.5, width = 0.75) +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x="Marker") +
  labs(y="% Giant Cells Positive") +
  guides(fill=FALSE) 
bulk

##..Plot number of GCs per Sample..##

freq_data<-read.csv("immune_cell_freqs.csv")
data_GC<-freq_data[freq_data$cell_type=="giant_cell",]
data_GC_TB<-data_GC[data_GC$Tissue=='mycobacteria',]

# order by sample
h_order<-c(21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,26,27,40,61,47,48,54,55,92,93)
data_GC_TB$SampleID <- factor(data_GC_TB$SampleID, levels=h_order)
data_GC_TB<-data_GC_TB[order(data_GC_TB$SampleID),]

# add annotation of sample type

AHRI<-c(21,84,42,88,28,89)
Stanford<-c(14,15,98,99,6,7,33,34,26,27,40,61,47,48,54,55,92,93)
Autopsy<-c(90,91,94,95,96,97)

data_GC_TB<-data_GC_TB %>% mutate(Origin=case_when(data_GC_TB$SampleID %in% AHRI ~ "AHRI",
                                                   data_GC_TB$SampleID %in% Stanford ~ "Stanford",
                                                   data_GC_TB$SampleID %in% Autopsy ~ "Autopsy"))

# plot
color_specimen<-c('#2E3192','#006838','#9E1F63') #resection, biopsy, autopsy
total_GC<-ggplot(data=data_GC_TB, aes(x=as.factor(SampleID),y=count, fill = Origin)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=color_specimen) +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x="Sample") +
  labs(y="$ Giant Cells") +
  guides(fill=FALSE) 
total_GC
