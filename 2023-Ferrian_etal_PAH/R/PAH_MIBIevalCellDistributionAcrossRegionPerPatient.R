# PAH_MIBIevalCellDistributionAcrossRegionPerPatient.R
# Author: Erin McCaffrey 
# Date created: 200322

library(ggplot2)
library(tidyverse)
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

#add annotation for sample type
ipah<-unique(data_immune[data_immune$Subgroup=='IPAH',]$PID)
hpah<-unique(data_immune[data_immune$Subgroup=='HPAH',]$PID)
hlt<-unique(data_immune[data_immune$Subgroup=='Healthy Controls',]$PID)
data_immune_summ$group<-'hlt'
data_immune_summ[data_immune_summ$PID %in% hpah, ]$group <- 'hpah'
data_immune_summ[data_immune_summ$PID %in% ipah, ]$group <- 'ipah'

##..Drop patient 21..##

data_immune_summ<-droplevels(data_immune_summ[!data_immune_summ$PID == 21, ])

##..Reorder regions..##

data_immune_summ$region_modified<-factor(data_immune_summ$region_modified, levels=c('vessel_associated',
                                                                                    'vessel_proximal',
                                                                                    'non_vascular'))

##..Reorder cell types..##

order<-c('Th','Tc','Bcell','NK','Macro','Mono','Neutro','DC')
data_immune_summ$cell_type<-factor(data_immune_summ$cell_type, levels=order)


##..Plot..##

colorkey<-read.csv('colorkey.csv')
colorkey<-droplevels(colorkey[colorkey$cell_lineage %in% levels(data_immune_summ$cell_type), ])
colorkey$cell_lineage<-factor(colorkey$cell_lineage, levels = levels(data_immune_summ$cell_type))
colorkey<-colorkey[order(colorkey$cell_lineage),]
color<-as.vector(colorkey$codes)

my_comparisons<-list(c("vessel_associated","non_vascular"),
                     c("vessel_associated","vessel_proximal"),
                     c("non_vascular","vessel_proximal"))

all<-ggplot(data = data_immune_summ, aes(x = region_modified, y = freq)) + 
  geom_boxplot(aes(fill = cell_type)) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method= "wilcox.test") +
  scale_fill_manual(values=color) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y')
all


pah<-ggplot(data = data_immune_summ[data_immune_summ$group %in% c('hpah','ipah'),], 
            aes(x = region_modified, y = freq)) + 
  geom_boxplot(aes(fill = cell_type), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method= "wilcox.test") +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Region') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y', ncol=4) + 
  theme(legend.position = 'none')
pah


# IPAH v HPAH

my_comparisons<-list(c("hpah","ipah"))

plot_data<-droplevels(data_immune_summ[data_immune_summ$group %in% c('ipah','hpah'),])

vessel<-ggplot(plot_data[plot_data$region_modified == 'vessel_associated',],aes(x=group, y=freq)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = cell_type), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'PAH') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y', ncol=4) +
  theme(legend.position = 'none')
vessel


vessel_proximal<-ggplot(plot_data[plot_data$region_modified == 'vessel_proximal',],aes(x=group, y=freq)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = cell_type), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'PAH') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y', ncol=4) +
  theme(legend.position = 'none')
vessel_proximal

non_vascular<-ggplot(plot_data[plot_data$region_modified == 'non_vascular',],aes(x=group, y=freq)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = cell_type), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'PAH') + 
  labs(y = 'Frequency') +
  facet_wrap(.~cell_type, scales='free_y', ncol=4) +
  theme(legend.position = 'none')
non_vascular
