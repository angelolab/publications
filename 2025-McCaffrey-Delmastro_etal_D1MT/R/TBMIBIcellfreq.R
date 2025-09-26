# TBMIBIcellfreq.R
# Author: Erin McCaffrey 
# Date created: 190318
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it plots the frequency of cell types across all granulomas with boxplots and stacked barplots. Done for both major
# lineage and then within the immune cell compartment.

library(dplyr)
library(viridis)
library(ggplot2)
library(forcats)
library(reshape2)
library(tidyr)

##..Import data..##
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv('NHP_cohort_data_norm_annotated.csv')

data_lineage<-data %>% select(SampleID,lineage)

data_immlineage<-data %>% select(SampleID,cell_lin)
data_immlineage<-droplevels(data_immlineage[!data_immlineage$cell_lin=="nonimmune",])

data_immune<-data %>% select(SampleID,cell_type)
data_immune<-droplevels(data_immune[!data_immune$cell_type  %in% c('epithelial','endothelial','fibroblast'),])


##..Create data frame with SampleID, cell_type, and count of that cell type in the given ROI..##

#lineage
Freqs_lineage <- table(data_lineage$SampleID, data_lineage$lineage)
Freqs_lineage <- as.data.frame(Freqs_lineage)
names(Freqs_lineage) <- c("SampleID","lineage","count")

#immune
Freqs_immune <- table(data_immune$SampleID, data_immune$cell_type)
Freqs_immune <- as.data.frame(Freqs_immune)
names(Freqs_immune) <- c("SampleID","cell_type","count")

#all cell types (total frequency)
Freqs_total <- table(data$SampleID, data$cell_type)

Freqs_total <- as.data.frame(Freqs_total)
names(Freqs_total) <- c("SampleID","cell_type","count")

#immune lineage
Freqs_immunelin <- table(data_immlineage$SampleID, data_immlineage$cell_lin)
Freqs_immunelin <- as.data.frame(Freqs_immunelin)
names(Freqs_immunelin) <- c("SampleID","cell_lin","count")

##..Get overall totals on a per sample basis..##
totals_lineage<-aggregate(Freqs_lineage$count, by=list(Category=Freqs_lineage$SampleID), FUN=sum)
totals_immune<-aggregate(Freqs_immune$count, by=list(Category=Freqs_immune$SampleID), FUN=sum)
totals_all<-aggregate(Freqs_total$count, by=list(Category=Freqs_total$SampleID), FUN=sum)
total_immune_lineage<-aggregate(Freqs_immunelin$count, by=list(Category=Freqs_immunelin$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

#lineage
for(i in unique(Freqs_lineage$SampleID)) {
  frequencies <- Freqs_lineage[Freqs_lineage$SampleID==i,"count"] / totals_lineage[totals_lineage$Category==i,2]
  Freqs_lineage[Freqs_lineage$SampleID==i,"frequency"] <- frequencies
}

#immune
for(i in unique(Freqs_immune$SampleID)) {
  frequencies_imm <- Freqs_immune[Freqs_immune$SampleID==i,"count"] / totals_immune[totals_immune$Category==i,2]
  Freqs_immune[Freqs_immune$SampleID==i,"frequency"] <- frequencies_imm
}

#all cell types
for(i in unique(Freqs_total$SampleID)) {
  frequencies_total <- Freqs_total[Freqs_total$SampleID==i,"count"] / totals_all[totals_all$Category==i,2]
  Freqs_total[Freqs_total$SampleID==i,"frequency"] <- frequencies_total
}

#immune lineage
for(i in unique(Freqs_immunelin$SampleID)) {
  frequencies_immlin <- Freqs_immunelin[Freqs_immunelin$SampleID==i,"count"] / total_immune_lineage[total_immune_lineage$Category==i,2]
  Freqs_immunelin[Freqs_immunelin$SampleID==i,"frequency"] <- frequencies_immlin
}

##..Create boxplot across cell types..##

imm_order<-c("Mac_Mono","CD4_T","neutrophil","imm_other","B_cell","CD8_T","T_other","mast","Treg","giant_cell")
lin_order<-c("immune","endothelial","fibroblast","epithelial")
immlin_order<-c("myeloid","lymphocyte","other","granulocyte")
point_order_all<-c(1,2,10,11,12,19,20,21,23,24,28,32,33,34,35,36,37,38,
                   8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)

##..Create a barplot of cell totals across FOVs..##

plot_data<-droplevels(totals_all[totals_all$Category %in% point_order_all,])

cell_bar<-ggplot(plot_data, aes(x=reorder(as.factor(Category),-x), y=x)) +
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_manual("grey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank()) 
cell_bar

#reorder by cell type and point order
Freqs_lineage$SampleID <- factor(Freqs_lineage$SampleID, levels=point_order_all)
Freqs_immune$SampleID <- factor(Freqs_immune$SampleID, levels=point_order_all)
Freqs_total$SampleID <- factor(Freqs_total$SampleID, levels=point_order_all)
Freqs_immunelin$SampleID <- factor(Freqs_immunelin$SampleID, levels=point_order_all)

Freqs_lineage$lineage <- factor(Freqs_lineage$lineage, levels=lin_order)
Freqs_immune$cell_type <- factor(Freqs_immune$cell_type, levels=rev(imm_order))
Freqs_immunelin$cell_lin <- factor(Freqs_immunelin$cell_lin, levels=immlin_order)

##..Subset the data to drop excluded samples..##
# point_order_all <- c(1,2,10,11,12,19,20,21,23,24,28,32,33,34,35,36,37,38,
#                    8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)

point_order_all <- c(23,24,33,34,35,10,11,12,19,20,21,1,2,36,37,38,28,32,
                     29,30,31,39,40,41,16,18,25,27,8,9,13,14,15)

Freqs_lineage <- Freqs_lineage[Freqs_lineage$SampleID %in% point_order_all,]
Freqs_immune <- Freqs_immune[Freqs_immune$SampleID %in% point_order_all,]
Freqs_immunelin <- Freqs_immunelin[Freqs_immunelin$SampleID %in% point_order_all,]

##..Get the color key..##

color_key <- read.csv("cell_color_key.csv")
plot_colors <- droplevels(color_key[color_key$Phenotype %in% imm_order,])
plot_colors$Phenotype <- factor(plot_colors$Phenotype, levels = imm_order)
plot_colors<-plot_colors[order(plot_colors$Phenotype),]
color<-as.vector(plot_colors$Hex)

##..Create box plots of pooled ROIs..##

lin_box<-ggplot(Freqs_lineage, aes(x=lineage, y=frequency, fill=lineage)) + geom_boxplot() +
  theme_bw() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
lin_box

imm_box<-ggplot(Freqs_immune, aes(x=fct_reorder(cell_type,frequency,.fun=median,.desc=TRUE), y=frequency, fill=cell_type)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = rev(color)) +
  theme(axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Immune Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
imm_box

##..Create bar plots lineage and immune for myco of pooled ROIs..##

lin_bar<-ggplot(Freqs_lineage, aes(x=lineage, y=frequency)) +
  theme_bw() +
  scale_fill_manual("white") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank()) +
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3)
lin_bar

immlin_bar<-ggplot(Freqs_immunelin, aes(x=cell_lin, y=frequency)) +
  theme_bw() +
  scale_fill_manual("white") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank()) +
  stat_summary(geom = "bar", fun.y = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3)
immlin_bar

##..Create stacked bar plot across ROIs..##

Freqs_lineage$lineage <- factor(Freqs_lineage$lineage, levels=rev(lin_order))
Freqs_lineage$SampleID <- factor(Freqs_lineage$SampleID, levels=point_order_all)
Freqs_lineage_<-Freqs_lineage[order(Freqs_lineage$SampleID),]


lin_bar<-ggplot(Freqs_lineage, aes(x=SampleID, y=frequency, fill=lineage)) + 
  theme_bw() +
  scale_fill_manual(values =c('#FF9966','#99FF99','#CC3333','#3366FF')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
lin_bar

#re-order to have the samples in order specified and immune phenos in order of decreasing median
Freqs_immune$SampleID <- factor(Freqs_immune$SampleID, levels=point_order_all)
Freqs_immune<-Freqs_immune[order(Freqs_immune$SampleID),]

Freqs_immune$cell_type <- factor(Freqs_immune$cell_type, levels=rev(imm_order))
Freqs_immune<-Freqs_immune[order(Freqs_immune$cell_type),]

plot_data <- Freqs_immune[Freqs_immune$SampleID %in% point_order_all,]

imm_bar<-ggplot(plot_data, aes(x=SampleID, y=frequency, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values=rev(color))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") 
imm_bar

##..Make a pie chart..##
imm_bar_pie<-ggplot(Freqs_immune, aes(x="", y=frequency, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values=rev(color))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") 
imm_bar_pie

##..Append sample type info and save frequency data..##

gran_D1MT<-c(7,8,9,13,14,15,16,17,18,25,26,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38) 

Freqs_immune<-Freqs_immune %>% mutate(Tissue=case_when(Freqs_immune$SampleID %in% gran_D1MT~ "IDO_Inhibitor",
                                                       Freqs_immune$SampleID %in% gran_con ~ "Control"))


Freqs_lineage<-Freqs_lineage %>% mutate(Tissue=case_when(Freqs_lineage$SampleID %in% gran_D1MT~ "IDO_Inhibitor",
                                                         Freqs_lineage$SampleID %in% gran_con ~ "Control"))

Freqs_total<-Freqs_total %>% mutate(Tissue=case_when(Freqs_total$SampleID %in% gran_D1MT~ "IDO_Inhibitor",
                                                     Freqs_total$SampleID %in% gran_con ~ "Control"))

Freqs_immunelin<-Freqs_immunelin %>% mutate(Tissue=case_when(Freqs_immunelin$SampleID %in% gran_D1MT~ "IDO_Inhibitor",
                                                             Freqs_immunelin$SampleID %in% gran_con ~ "Control"))


write.csv(Freqs_immune, file="immune_cell_freqs.csv",row.names = FALSE)
write.csv(Freqs_lineage, file="lineage_freqs.csv",row.names = FALSE)
write.csv(Freqs_total, file="total_freqs.csv",row.names = FALSE)
write.csv(Freqs_immunelin, file="immlineage_freqs.csv",row.names = FALSE)

write.csv(totals_lineage, file="lineage_cell_totals.csv",row.names = FALSE)
write.csv(totals_immune, file="immune_cell_totals.csv",row.names = FALSE)

