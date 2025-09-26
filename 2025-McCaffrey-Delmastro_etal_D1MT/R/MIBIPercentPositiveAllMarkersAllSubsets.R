# MIBIPercentPositiveAllMarkersAllSubsets.R
# Author: Erin McCaffrey 
# Date created: 201129

require(dplyr)
library(ggplot2)
library(tibble)
library(ggpubr)

##..Import data..##

setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv("NHP_cohort_data_norm_annotated.csv")

##..Create a dataframe with all points and all cells counts..##

cell_counts<-as.data.frame(table(data$SampleID, data$cell_type))
colnames(cell_counts)<-c("SampleID","cell_type","cell_count")
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
  freq_marker_patient<-as.data.frame(table(data_marker$SampleID,data_marker$cell_type))
  #append positive counts and percent positive to summary dataframes
  cell_counts[[marker]]<-0
  cell_counts[cell_counts$SampleID %in% freq_marker_patient$Var1 & cell_counts$cell_type
              %in% freq_marker_patient$Var2,][[marker]]<-freq_marker_patient$Freq
  cell_freqs[[marker]]<-as.numeric(format((cell_counts[[marker]] / cell_counts$cell_count)),digits=3)
}

##..Get stats for total sublineage..##

cell_counts_sublin<-as.data.frame(table(data$SampleID, data$cell_lin))
colnames(cell_counts_sublin)<-c("SampleID","cell_type","cell_count")
cell_freqs_sublin<-cell_counts_sublin

for(i in 1:length(markers)) {
  #define marker and threshold
  marker<-markers[i]
  threshold<-thresholds[i]
  #subset positive cells
  data_marker<-data  %>% filter(get(marker)>=threshold)
  #get table of positive cell counts per patient per cell type
  freq_marker_patient<-as.data.frame(table(data_marker$SampleID,data_marker$cell_lin))
  #append positive counts and percent positive to summary dataframes
  cell_counts_sublin[[marker]]<-0
  cell_counts_sublin[cell_counts_sublin$SampleID %in% freq_marker_patient$Var1 & cell_counts_sublin$cell_type
              %in% freq_marker_patient$Var2,][[marker]]<-freq_marker_patient$Freq
  cell_freqs_sublin[[marker]]<-as.numeric(format((cell_counts_sublin[[marker]] / cell_counts_sublin$cell_count)),digits=3)
}

##..Get stats for total immune..##

data_immune<-droplevels(data[data$lineage == "immune",])
cell_counts_imm<-as.data.frame(table(data_immune$SampleID))
colnames(cell_counts_imm)<-c("SampleID","cell_count")
cell_freqs_imm<-cell_counts_imm

for(i in 1:length(markers)) {
  #define marker and threshold
  marker<-markers[i]
  threshold<-thresholds[i]
  #subset positive cells
  data_marker<-data_immune  %>% filter(get(marker)>=threshold)
  #get table of positive cell counts per patient per cell type
  freq_marker_patient<-as.data.frame(table(data_marker$SampleID))
  #append positive counts and percent positive to summary dataframes
  cell_counts_imm[[marker]]<-0
  cell_counts_imm[cell_counts_imm$SampleID %in% freq_marker_patient$Var1,][[marker]]<-freq_marker_patient$Freq
  cell_freqs_imm[[marker]]<-as.numeric(format((cell_counts_imm[[marker]] / cell_counts_imm$cell_count)),digits=3)
}

cell_counts_imm<-add_column(cell_counts_imm, cell_type = "immune", .after = 1)
cell_freqs_imm<-add_column(cell_freqs_imm, cell_type = "immune", .after = 1)

##..Merge into single data set..##

cell_counts_combo<-do.call("rbind", list(cell_counts, cell_counts_sublin, cell_counts_imm))
cell_freqs_combo<-do.call("rbind", list(cell_freqs, cell_freqs_sublin, cell_freqs_imm))

##..Add group and subgroup annotation..##

D1MT<-unique(data[data$Tissue=="gran_D1MT",]$SampleID)
con<-unique(data[data$Tissue=="gran_control",]$SampleID)

cell_counts_combo$Group<-'D1MT'
cell_counts_combo[cell_counts_combo$SampleID %in% con,]$Group<-'control'
cell_freqs_combo$Group<-'D1MT'
cell_freqs_combo[cell_freqs_combo$SampleID %in% con,]$Group<-'control'

##..Export csv..##

write.csv(cell_counts_combo,"count-of-positivecells-persubset.csv",row.names = FALSE)
write.csv(cell_freqs_combo,"freq-of-positivecells-persubset.csv",row.names = FALSE)

##..Plot..##
freq_data <- cell_freqs_combo

##..Drop excluded samples..##

drop_samples <- c(17,22,26)
freq_data_cohort <- freq_data[!freq_data$SampleID %in% drop_samples,]
data_cohort <- data[!data$SampleID %in% drop_samples,]

##..Choose marker and plot D1MT v control for all subsets..##

plot_data<-droplevels(cell_freqs_combo[!cell_freqs_combo$cell_type=="giant_cell",])
ggplot(data = plot_data, aes(x = Group, y = Ki67, fill = Group)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency') +
  facet_wrap(~cell_type,scales=c("free")) 

plot_data<-droplevels(freq_data_cohort[freq_data_cohort$cell_type=="CD4_T",])
ggplot(data = plot_data, aes(x = Group, y = Ki67, fill = Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency') 

data_noGC<-droplevels(data[!data$cell_type=="giant_cell",])
ggplot(data = data_noGC, aes(x = Tissue, y =pS6, fill = Tissue)) + 
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="black") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency') +
  facet_wrap(~cell_type,scales=c("free")) 

plot_data<-droplevels(data_cohort[data_cohort$cell_type=="Mac_Mono",])
plot_data<-droplevels(data_cohort[data_cohort$lineage=="immune",])
# plot_data <-data
ggplot(data = plot_data, aes(x = Tissue, y = IDO, fill = Tissue)) + 
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  stat_summary(fun=median, geom="point", size=2, color="black") +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Tissue') + 
  labs(y = 'IDO1 Expression') 



