# PAH_MIBIevalCellDistributionPerPatient_ImmuneBaseline.R
# Author: Erin McCaffrey 
# Date created: 200331

library(ggplot2)
library(tidyverse)
library(forcats)
library(dplyr)
library(ggpubr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")

##..Subset just the immune cells..##

data_immune<-droplevels(data[!data$cell_lineage %in% c('fibroblast','endothelial','mesenchymal','epithelial'),])

##..Determine the frequency of a given cell type in that region out of all cells of that type in the patient.##

data_immune_summ <- as.data.frame(table(data_immune$cell_lineage, data_immune$region_modified, data_immune$PID))
colnames(data_immune_summ) <- c("cell_type","region_modified","PID","count")

#get cell totals for each type across samples
cell_totals <- as.data.frame(table(data_immune$PID,data_immune$cell_lineage))
colnames(cell_totals) <- c("PID","cell_type","total")

#sort and merge dataframes
cell_totals <- cell_totals %>% slice(rep(1:n(), each = 3))
cell_totals <- cell_totals[order(cell_totals$PID), ] 
cell_totals <- cell_totals[order(cell_totals$cell_type), ] 
data_immune_summ <- data_immune_summ[order(data_immune_summ$PID), ]
data_immune_summ <- data_immune_summ[order(data_immune_summ$cell_type), ] 
data_immune_summ$cell_total <- cell_totals$total

#determine frequency per sample per region
data_immune_summ$freq <- data_immune_summ$count / data_immune_summ$cell_total


##..Determine the baseline immune cell distribution..##

data_region_summ <- as.data.frame(table(data_immune$region_modified, data_immune$PID))
colnames(data_region_summ) <- c("region_modified","PID","count")

#get cell totals for each region across samples
totals <- as.data.frame(table(data_immune$PID))
colnames(totals) <- c("PID","total")

#sort and merge dataframes
totals <- totals %>% slice(rep(1:n(), each = 3))
totals <- totals[order(totals$PID), ] 
data_region_summ<- data_region_summ[order(data_region_summ$PID), ]
data_region_summ$patient_total <- totals$total

#determine frequency per sample per region
data_region_summ$freq <- data_region_summ$count / data_region_summ$patient_total

##..Append the baseline data to the cell-type specific dataframe..##

data_immune_summ$imm_baseline<-rep(data_region_summ$freq, times=8)

##..Get log2 fold-change for cell type relative to total immune baseline..##

data_immune_summ$FC<-log2(data_immune_summ$freq / data_immune_summ$imm_baseline)

##..Append the patient group..##

#add annotation for sample type
ipah<-unique(data_immune[data_immune$Subgroup=='IPAH',]$PID)
hpah<-unique(data_immune[data_immune$Subgroup=='HPAH',]$PID)
hlt<-unique(data_immune[data_immune$Subgroup=='Healthy Controls',]$PID)
data_immune_summ$group<-'hlt'
data_immune_summ[data_immune_summ$PID %in% hpah, ]$group <- 'hpah'
data_immune_summ[data_immune_summ$PID %in% ipah, ]$group <- 'ipah'

##..Drop patient 21..##

data_immune_summ<-droplevels(data_immune_summ[!data_immune_summ$PID == 21, ])

##..Import color key..##

colorkey<-read.csv('colorkey.csv')
colorkey<-droplevels(colorkey[colorkey$cell_lineage %in% levels(data_immune_summ$cell_type), ])
colorkey$cell_lineage<-factor(colorkey$cell_lineage, levels = levels(data_immune_summ$cell_type))
colorkey<-colorkey[order(colorkey$cell_lineage),]
color<-as.vector(colorkey$codes)


##..Plot..##

#pah only
plot_data<- droplevels(data_immune_summ[data_immune_summ$group %in% c('hpah','ipah'),])
plot_data<-plot_data[is.finite(plot_data$FC),]

vessel_data<- droplevels(plot_data[plot_data$region == 'vessel_associated',])
vessel<-ggplot(vessel_data, 
               aes(fct_reorder(cell_type,FC,.fun=median, .desc=TRUE), y=FC, fill=cell_type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="log2(Cell Type Freq / Immune Cell Freq)") +
  theme(legend.position = 'none') +
  ggtitle("Vessel-Associated") +
  guides(fill=guide_legend(title="Cell Type"))
vessel

vessel_prox_data<- droplevels(plot_data[plot_data$region == 'vessel_proximal',])
vessel_prox<-ggplot(vessel_prox_data, 
               aes(fct_reorder(cell_type,FC,.fun=median, .desc=TRUE), y=FC, fill=cell_type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="log2(Cell Type Freq / Immune Cell Freq)") +
  theme(legend.position = 'none') +
  ggtitle("Vessel-Proximal") +
  guides(fill=guide_legend(title="Cell Type"))
vessel_prox

non_vasc_data<- droplevels(plot_data[plot_data$region == 'non_vascular',])
non_vessel<-ggplot(non_vasc_data, 
                    aes(fct_reorder(cell_type,FC,.fun=median, .desc=TRUE), y=FC, fill=cell_type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="log2(Cell Type Freq / Immune Cell Freq)") +
  theme(legend.position = 'none') +
  ggtitle("Non-Vascular") +
  guides(fill=guide_legend(title="Cell Type"))
non_vessel



