# PAH_neighborhoodFreqsPerPatient.R
# Author: Erin McCaffrey 
# Date created: 200403
# Overview: This script reads in the csv for the PAH dataset, determines the percent of 
# total for all patients of the neighborhoods
library(dplyr)
library(ggplot2)
library(RColorBrewer)


##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("PAH_AllCells_Neighborhood_K=10.csv")
cell_data<-read.csv("celldata_region_annotated.csv")


##..Add patient data for point data..##
pah_data<-droplevels(cell_data[cell_data$Group=='PAH',])
data$PID <- pah_data$PID

# write.csv(data,"PAH_AllCells_Neighborhood_K=10_PatientAnnotated.csv", row.names = FALSE)

##..Create a table of cell frequencies per point..##

neighborhood_freqs<-as.data.frame(table(data$Point_num,data$cluster))
colnames(neighborhood_freqs)<-c('Point','Cluster','Count')
totals_patient<-as.numeric(table(data$Point_num))
neighborhood_freqs$total<-rep(totals_patient,10)
neighborhood_freqs$freq<-as.numeric(neighborhood_freqs$Count/neighborhood_freqs$total)

##..Stacked bars all points..##

ggplot(neighborhood_freqs, aes(x=Point, y=freq, fill=Cluster)) + 
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") 

##..Create a table of cell frequencies per patient..##

neighborhood_freqs_patient<-as.data.frame(table(data$PID,data$cluster))
colnames(neighborhood_freqs_patient)<-c('PID','Cluster','Count')
totals_patient<-as.numeric(table(data$PID))
neighborhood_freqs_patient$total<-rep(totals_patient,10)
neighborhood_freqs_patient$freq<-as.numeric(neighborhood_freqs_patient$Count/neighborhood_freqs_patient$total)

##..Stacked bars all points..##

ggplot(neighborhood_freqs_patient, aes(x=PID, y=freq, fill=Cluster)) + 
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") 

##..Assess frequency between HPAH and IPAH..##

# get the IPAH and HPAH labels

ipah<-unique(pah_data[pah_data$Subgroup=='IPAH',]$PID)
hpah<-unique(pah_data[pah_data$Subgroup=='HPAH',]$PID)
neighborhood_freqs_patient$group<-'hpah'
neighborhood_freqs_patient[neighborhood_freqs_patient$PID %in% ipah, ]$group <- 'ipah'

my_comparisons<-list(c('hpah','ipah'))

ggplot(neighborhood_freqs_patient,aes(x=group, y=freq)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Cluster), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'PAH') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Cluster, scales='free_y', ncol=5) +
  theme(legend.position = 'none')

# ##..Export neighborhood colorkey..##
# 
# 
# colorkey<-as.data.frame(c(1,2,3,4,5,6,7,8,9,10))
# colorkey$code<-brewer.pal(n=10, "Paired")
# colnames(colorkey)<-c('neighborhood','code')
# write.csv(colorkey, "neighborhood_colorkey.csv", row.names = FALSE)
