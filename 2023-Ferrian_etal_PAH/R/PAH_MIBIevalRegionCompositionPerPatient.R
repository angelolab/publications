# PAH_MIBIevalRegionCompositionPerPatient.R
# Author: Erin McCaffrey 
# Date created: 200310

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")

##..Subset just the immune cells..##

data_immune<-droplevels(data[!data$cell_lineage %in% c('fibroblast','endothelial','mesenchymal','epithelial'),])

##..Determine the frequency of immune cells across all zones..##

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

#add annotation for sample type
ipah<-unique(data_immune[data_immune$Subgroup=='IPAH',]$PID)
hpah<-unique(data_immune[data_immune$Subgroup=='HPAH',]$PID)
hlt<-unique(data_immune[data_immune$Subgroup=='Healthy Controls',]$PID)
data_region_summ$group<-'hlt'
data_region_summ[data_region_summ$PID %in% hpah, ]$group <- 'hpah'
data_region_summ[data_region_summ$PID %in% ipah, ]$group <- 'ipah'

##..Drop patient 21..##

data_region_summ<-droplevels(data_region_summ[!data_region_summ$PID == 21, ])

##..Reorder regions..##

data_region_summ$region_modified<-factor(data_region_summ$region_modified, levels=c('vessel_associated',
                                                                                    'vessel_proximal',
                                                                                    'non_vascular'))
##..Reorder cell types..##

order<-c('Th','Tc','Bcell','NK','Macro','Mono','Neutro','DC')
data_region_summ$cell_type<-factor(data_region_summ$cell_type, levels=order)

##..Plot mean and points for distribution..##

my_comparisons<-list(c("vessel_associated","non_vascular"),
                     c("vessel_associated","vessel_proximal"),
                     c("non_vascular","vessel_proximal"))

data_summary <- droplevels(data_region_summ[data_region_summ$group %in% c('hpah','ipah'),]) %>%
  group_by(region_modified) %>%
  summarize(combo_median = median(freq),
            combo_se = sqrt(var(freq)/length(freq)))

ggplot(data_region_summ[data_region_summ$group %in% c('hpah','ipah'),],
       aes(x =region_modified, y = freq)) + 
  geom_point(aes(x = region_modified, y = freq), 
             position = position_jitter(width = 0.1, height = 0.0),
             size = 4) + 
  geom_point(aes(y = combo_median), color = "black", size = 2, data =data_summary) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method= "wilcox.test")

##..Determine the frequency of a given cell type in that region out of all cells in that region..##

data_immune_summ <- as.data.frame(table(data_immune$cell_lineage, data_immune$region_modified, data_immune$PID))
colnames(data_immune_summ) <- c("cell_type","region_modified","PID","count")

#get cell totals for each region across samples
region_totals <- as.data.frame(table(data_immune$PID,data_immune$region_modified))
colnames(region_totals) <- c("PID","region","total")

#sort and merge dataframes
region_totals <- region_totals %>% slice(rep(1:n(), each = 8))
region_totals <- region_totals[order(region_totals$PID), ] 
region_totals <- region_totals[order(region_totals$region), ] 
data_immune_summ <- data_immune_summ[order(data_immune_summ$PID), ]
data_immune_summ <- data_immune_summ[order(data_immune_summ$region_modified), ] 
data_immune_summ$region_total <- region_totals$total


#determine frequency per sample per region
data_immune_summ$freq <- data_immune_summ$count / data_immune_summ$region_total

#add annotation for sample type
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
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method= "t.test") +
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


# Region composition in PAH (pooled)

vessel_box<-ggplot(plot_data[plot_data$region == 'vessel_associated',], 
                aes(x=fct_reorder(cell_type,freq,.fun=median,.desc=TRUE), y=freq, fill=cell_type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Immune Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
vessel_box

vessel_proximal_box<-ggplot(plot_data[plot_data$region == 'vessel_proximal',], 
                   aes(x=fct_reorder(cell_type,freq,.fun=median,.desc=TRUE), y=freq, fill=cell_type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Immune Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
vessel_proximal_box

perivascular_box<-ggplot(plot_data[plot_data$region == 'perivascular',], 
                            aes(x=fct_reorder(cell_type,freq,.fun=median,.desc=TRUE), y=freq, fill=cell_type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Immune Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
perivascular_box

nonvascular_box<-ggplot(plot_data[plot_data$region == 'non_vascular',], 
                         aes(x=fct_reorder(cell_type,freq,.fun=median,.desc=TRUE), y=freq, fill=cell_type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Immune Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
nonvascular_box

# Stacked bar
order<-c('Tc','Th','Mono','DC','Neutro','Macro','NK','Bcell')

plot_data$cell_type <- factor(plot_data$cell_type, levels=rev(order))
plot_data<-plot_data[order(plot_data$cell_type),]

colorkey<-droplevels(colorkey[colorkey$cell_lineage %in% levels(plot_data$cell_type), ])
colorkey$cell_lineage<-factor(colorkey$cell_lineage, levels = levels(plot_data$cell_type))
colorkey<-colorkey[order(colorkey$cell_lineage),]
color<-as.vector(colorkey$codes)

vesselbar<-ggplot(plot_data[plot_data$region_modified == 'vessel_associated',], 
                  aes(x=PID, y=freq, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = 'none') +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency of Total Vessel Cells") +
  guides(fill=guide_legend(title="Cell Type"))
vesselbar


vessel_proximalbar<-ggplot(plot_data[plot_data$region_modified == 'vessel_proximal',], 
                           aes(x=PID, y=freq, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = 'none') +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency of Total Vessel Proximal Cells") +
  guides(fill=guide_legend(title="Cell Type"))
vessel_proximalbar

non_vascularbar<-ggplot(plot_data[plot_data$region_modified == 'non_vascular',], 
                        aes(x=PID, y=freq, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = 'none') +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency of Total Non-Vascular Cells") +
  guides(fill=guide_legend(title="Cell Type"))
non_vascularbar