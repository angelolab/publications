# TBMIBICD4CD8ratio.R
# Author: Erin McCaffrey 
# Date created: 190501
# Overview: This script reads in the frequency data for granuloma immune cell population. It then determines
# the ratio of CD4 T cells to CD8 T cells (log2(CD4/CD8)) and plots them in descending order (ie. more CD4s on 
# the left side and less CD4s on the right side due to higher CD8). Runs downstream stats based on ratio.

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(ggpmisc)

##..Data importation, clean-up and reshaping..## 

# import
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/D1MT_NHP_TB/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
data<-read.csv("immune_cell_freqs.csv")

# drop excluded samples
drop_samples <- c(17,22,26)
data <- data[!data$SampleID %in% drop_samples,]

# reshape matrix so that each row is a point, each column is a frequency for that cell type
freq_data <- reshape2::dcast(data, SampleID ~ cell_type, value.var = "frequency") #frequency data
count_data <- reshape2::dcast(data, SampleID ~ cell_type, value.var = "count") #count data

# create T cell count only matrix 
tcell_data<-count_data[,c(1,3,4)]

# annotate origin of sample based on sample number
gran_D1MT<-c(8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,23,24,28,32,33,34,35,36,37,38) 

tcell_data<-tcell_data %>% mutate(Origin=case_when(tcell_data$SampleID %in% gran_D1MT ~ "gran_D1MT",
                                                   tcell_data$SampleID %in% gran_con ~ "gran_con"))
freq_data<-freq_data %>% mutate(Origin=case_when(freq_data$SampleID %in% gran_D1MT ~ "gran_D1MT",
                                                 freq_data$SampleID %in% gran_con ~ "gran_con"))
count_data<-count_data %>% mutate(Origin=case_when(count_data$SampleID %in% gran_D1MT ~ "gran_D1MT",
                                                   count_data$SampleID %in% gran_con ~ "gran_con"))

# append the NHP ID
data_gran <- read.csv('NHP_cohort_data_norm_annotated.csv')
sample_key <- unique(data_gran[,c('SampleID','PatientID')])
sample_key[sample_key$SampleID == 15, ]$PatientID <- '16A105'
tcell_data <- left_join(tcell_data, sample_key, by = c('SampleID'))
freq_data <- left_join(freq_data, sample_key, by = c('SampleID'))
count_data <- left_join(count_data, sample_key, by = c('SampleID'))

##..Calculate log2(CD4/CD8)..##

log2fc<-log2(tcell_data$CD4_T/tcell_data$CD8_T)
tcell_data$FC<-log2fc

##..Plot ratio and color by sample of origin or sample site..##

# T cell ratio in all TB FOVs in descending order
tcell_ratio<-ggplot(data=tcell_data, aes(x=reorder(SampleID, -FC), y=FC)) +
  geom_bar(stat="Identity", aes(fill=as.factor(Origin))) +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  # theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="CD4:CD8 Ratio (log2 [CD4 T cells / CD8 T cells])")
tcell_ratio

##..Plot and compare the absolute CD4 and CD8 counts as well as the frequencies..##

# Total CD4 counts per sample
CD4_count_persample<-ggplot(data=tcell_data, aes(x=reorder(SampleID, -CD4_T), y=CD4_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total Number of CD4 T cells") 
CD4_count_persample


# Total CD8 counts per sample
CD8_count_persample<-ggplot(data=tcell_data, aes(x=reorder(SampleID, -CD8_T), y=CD8_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total Number of CD8 T cells")
CD8_count_persample


# Origin comparisons: frequency

# CD4 T cells
data_summary <- freq_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD4_T))

CD4_freq_bulk<-ggplot(data=freq_data, aes(x=as.factor(Origin), y=CD4_T, color=as.factor(Origin))) +
  geom_point(data=freq_data, aes(x=as.factor(Origin), 
                                 y=CD4_T, 
                                 fill=as.factor(Origin), 
                                 color=as.factor(Origin),
                                 shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total immune cells)")
CD4_freq_bulk


# CD8 T cells
data_summary <- freq_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD8_T))


CD8_freq_bulk <- ggplot(data=freq_data, aes(x=as.factor(Origin), y=CD8_T, color=as.factor(Origin))) +
  geom_point(data=freq_data, aes(x=as.factor(Origin), 
                                 y=CD8_T, 
                                 fill=as.factor(Origin), 
                                 color=as.factor(Origin),
                                 shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total immune cells)")
CD8_freq_bulk

# fold-change
data_summary <- tcell_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(FC))


FC_group<-ggplot(data=tcell_data[!tcell_data$FC == Inf,], aes(x=as.factor(Origin), y=FC, color=as.factor(Origin))) +
  geom_point(data=tcell_data[!tcell_data$FC == Inf,],
             aes(x=as.factor(Origin), 
                 y=FC, 
                 fill=as.factor(Origin), 
                 color=as.factor(Origin),
                 shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y="CD4:CD8 Ratio (log2 [CD4 T cells / CD8 T cells])")
FC_group


# Origin comparisons: count

# CD4 T cells
data_summary <- tcell_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD4_T))


CD4_count<-ggplot(data=tcell_data, aes(x=as.factor(Origin), y=CD4_T, color=as.factor(Origin))) +
  geom_point(data=tcell_data, aes(x=as.factor(Origin), 
                                  y=CD4_T, 
                                  fill=as.factor(Origin), 
                                  color=as.factor(Origin),
                                  shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y="Count CD4 T") 
CD4_count


# CD8 T cells
data_summary <- tcell_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD8_T))


CD8_count<-ggplot(data=tcell_data, aes(x=as.factor(Origin), y=CD8_T, color=as.factor(Origin))) +
  geom_point(data=tcell_data, 
             aes(x=as.factor(Origin), 
                 y=CD8_T, 
                 fill=as.factor(Origin), 
                 color=as.factor(Origin),
                 shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y="Count CD8 T")
CD8_count

# Evaluate T cell frequencies of total T cells

tcell_count_data <- count_data[,c(1,3,4,10:13)]
tcell_count_data$total_T <- rowSums(tcell_count_data[,c(2:5)])
tcell_freq_data <- tcell_count_data
tcell_freq_data[,c(2:5,8)] <- tcell_freq_data[,c(2:5,8)]/tcell_count_data$total_T

# CD4 T cells
data_summary <- tcell_freq_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD4_T))

CD4_freq_tcell<-ggplot(data=tcell_freq_data, aes(x=as.factor(Origin), y=CD4_T, color=as.factor(Origin))) +
  geom_point(data=tcell_freq_data, aes(x=as.factor(Origin), 
                                 y=CD4_T, 
                                 fill=as.factor(Origin), 
                                 color=as.factor(Origin),
                                 shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total T cells)")
CD4_freq_tcell

# CD8 T cells
data_summary <- tcell_freq_data %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD8_T))

CD8_freq_tcell<-ggplot(data=tcell_freq_data, aes(x=as.factor(Origin), y=CD8_T, color=as.factor(Origin))) +
  geom_point(data=tcell_freq_data, aes(x=as.factor(Origin), 
                                       y=CD8_T, 
                                       fill=as.factor(Origin), 
                                       color=as.factor(Origin),
                                       shape = as.factor(PatientID)), 
             position = position_jitter(width = 0.1, height = 0.0), 
             size = 3) + 
  scale_fill_manual(values = c('#343434','#323CA1')) +
  scale_color_manual(values = c('#343434','#323CA1')) +
  scale_shape_manual(values = c(8,21,22,4,23,24,25)) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total T cells)")
CD8_freq_tcell
