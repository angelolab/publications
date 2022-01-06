# TBMIBI_SuppCells_IDO1PDL1.R
# Author: Erin McCaffrey 
# Date created: 200306
# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# IDO and PDL1 pos cells in neutrophils and epithelium

library(ggplot2)

##..Import data..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Subset just TB samples..##

data_norm<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])

##..Get percent IDO and PDL1 positive.##

IDO_thresh = 0.26
PDL1_thresh = 0.25

##..Get total number of IDO1 and PDL1+ neutrophils..##

neut<-as.data.frame(as.numeric(sum(data_norm$cell_type=="neutrophil")))
colnames(neut)<-"neut_total"
neut$IDOpos<-as.numeric(sum(data_norm[data_norm$cell_type=="neutrophil",]$IDO>IDO_thresh))
neut$PDL1pos<-as.numeric(sum(data_norm[data_norm$cell_type=="neutrophil",]$PD.L1>PDL1_thresh))

##..Get percent IDO1 and PDL1 positive..##

neut$IDOpos_percent<-neut$IDOpos/neut$neut_total
neut$PDL1pos_percent<-neut$PDL1pos/neut$neut_total

markers<-c("IDO","PDL1")
percent<-c(neut$IDOpos_percent,neut$PDL1pos_percent)
neut_markers<-data.frame(markers,percent)

##..Get total number of IDO1 and PDL1+ epithelial..##

epi<-as.data.frame(as.numeric(sum(data_norm$cell_type=="epithelial")))
colnames(epi)<-"epi_total"
epi$IDOpos<-as.numeric(sum(data_norm[data_norm$cell_type=="epithelial",]$IDO>IDO_thresh))
epi$PDL1pos<-as.numeric(sum(data_norm[data_norm$cell_type=="epithelial",]$PD.L1>PDL1_thresh))

##..Get percent IDO1 and PDL1 positive..##

epi$IDOpos_percent<-epi$IDOpos/epi$epi_total
epi$PDL1pos_percent<-epi$PDL1pos/epi$epi_total

markers<-c("IDO","PDL1")
percent<-c(epi$IDOpos_percent,epi$PDL1pos_percent)
epi_markers<-data.frame(markers,percent)

##..Plot..##

neut<-ggplot(data=neut_markers, aes(x=markers,y=percent)) +
  geom_bar(stat="identity", fill="grey", lwd=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="Marker") +
  labs(y="% neutrophils") +
  guides(fill=FALSE) 
neut

epi<-ggplot(data=epi_markers, aes(x=markers,y=percent)) +
  geom_bar(stat="identity", fill="grey", lwd=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="Marker") +
  labs(y="% epithelials") +
  guides(fill=FALSE) 
epi

