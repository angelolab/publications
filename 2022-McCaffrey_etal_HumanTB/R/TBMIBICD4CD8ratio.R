# TBMIBICD4CD8ratio.R
# Author: Erin McCaffrey 
# Date created: 190501
# Overview: This script reads in the frequency data for granuloma immune cell population. It then determines
# the ratio of CD4 T cells to CD8 T cells (log2(CD4/CD8)) and plots them in descending order (ie. more CD4s on 
# the left side and less CD4s on the right side due to higher CD8). Runs downstream stats based on ratio.

library(ggplot2)
library(dplyr) # make sure plyr is not loaded or dplyr is loaded first to prevent issue ith summarize()
library(forcats)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(ggpmisc)

##..Data importation, clean-up and reshaping..## 
data<-read.csv("immune_cell_freqs.csv")

# reshape matrix so that each row is a point, each column is a frequency for that cell type
freq_data <- reshape2::dcast(data, SampleID ~ cell_type, value.var = "frequency") #frequency data
count_data <- reshape2::dcast(data, SampleID ~ cell_type, value.var = "count") #count data

# create T cell count only matrix 
keep_cols<-c('SampleID','CD4_T','CD8_T')
tcell_data<-count_data[,keep_cols]

# add sample code information to all matrices
Samp_num<-c(30,30,17,17,2,20,20,4,18,18,31,3,29,29,12,12,31,67,68,69,70,71,72,73,74,75,76,2,3,4,
            40,40,41,41,42,42,43,43,16,16)
tcell_data$samp_num<-Samp_num
freq_data$samp_num<-Samp_num
count_data$samp_num<-Samp_num

#annotate origin of sample based on sample number
AHRI<-c(2,3,4)
Stanford<-c(17,30,18,20,31,29,12,41,16)
Sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
Autopsy<-c(40,42,43)

tcell_data<-tcell_data %>% mutate(Origin=case_when(tcell_data$samp_num %in% AHRI ~ "AHRI",
                                                   tcell_data$samp_num %in% Stanford ~ "Stanford",
                                                   tcell_data$samp_num %in% Sarcoid ~ "Sarcoid",
                                                   tcell_data$samp_num %in% Autopsy ~ "Autopsy"))
freq_data<-freq_data %>% mutate(Origin=case_when(freq_data$samp_num %in% AHRI ~ "AHRI",
                                                  freq_data$samp_num %in% Stanford ~ "Stanford",
                                                 freq_data$samp_num %in% Sarcoid ~ "Sarcoid",
                                                 tcell_data$samp_num %in% Autopsy ~ "Autopsy"))
count_data<-count_data %>% mutate(Origin=case_when(count_data$samp_num %in% AHRI ~ "AHRI",
                                                   count_data$samp_num %in% Stanford ~ "Stanford",
                                                   count_data$samp_num %in% Sarcoid ~ "Sarcoid",
                                                   tcell_data$samp_num %in% Autopsy ~ "Autopsy"))
##..Add the organ site..##

eptb<-c(30,18,20,31,29,12,41) #extrrapulm

tcell_data$organ<-'lung'
tcell_data[tcell_data$samp_num %in% eptb,]$organ<-'expulm'
tcell_data[tcell_data$samp_num %in% Sarcoid,]$organ<-'sarcoid'

freq_data$organ<-'lung'
freq_data[freq_data$samp_num %in% eptb,]$organ<-'expulm'
freq_data[freq_data$samp_num %in% Sarcoid,]$organ<-'sarcoid'

count_data$organ<-'lung'
count_data[count_data$samp_num %in% eptb,]$organ<-'expulm'
count_data[count_data$samp_num %in% Sarcoid,]$organ<-'sarcoid'

##..Calculate log2(CD4/CD8)..##

log2fc<-log2(tcell_data$CD4_T/tcell_data$CD8_T)
tcell_data$FC<-log2fc

##..Plot ratio and color by sample of origin or sample site..##

# TB only
tcell_data_TB<-droplevels(tcell_data[!tcell_data$Origin=="Sarcoid",])

# TB and Sarcoid
tcell_data_all<-tcell_data
tcell_data_all<-tcell_data_all %>% mutate(Tissue=case_when(tcell_data_all$samp_num %in% AHRI ~ "TB",
                                                           tcell_data_all$samp_num %in% Stanford ~ "TB",
                                                           tcell_data_all$samp_num %in% Autopsy ~ "TB",
                                                           tcell_data_all$samp_num %in% Sarcoid ~ "Sarcoid"))

# Assign color palettes

color_disease<-c("#00A59C","#9CD9D3") #TB and sarcoid
color_specimen<-c('#2E3192','#006838','#9E1F63') #resection, biopsy, autopsy

# T cell ratio in all TB FOVs in descending order
tcell_ratio<-ggplot(data=tcell_data_TB, aes(x=reorder(SampleID, -FC), y=FC)) +
  geom_bar(stat="Identity", aes(fill=Origin)) +
  scale_fill_manual(values=color_specimen) +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="CD4:CD8 Ratio (log2 [CD4 T cells / CD8 T cells])")
tcell_ratio

# CD4:CD8 Ratio in Sarcoid v TB
ratio_sarcvTB<-ggplot(data=tcell_data_all[!tcell_data_all$FC == Inf,], aes(x=Tissue, y=FC, fill=Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() + 
  scale_fill_manual(values=color_disease) +
  stat_compare_means(data=tcell_data_all,label = "p.signif", method= "wilcox.test") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(legend.position = 'none') +
  labs(y="FC Sarcoid v TB")
ratio_sarcvTB

##..Plot and compare the absolute CD4 and CD8 counts..##

# Total CD4 counts per sample
CD4_count_persample<-ggplot(data=tcell_data_TB, aes(x=reorder(SampleID, -CD4_T), y=CD4_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_specimen) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total Number of CD4 T cells") 
CD4_count_persample


# Total CD8 counts per sample
CD8_count_persample<-ggplot(data=tcell_data_TB, aes(x=reorder(SampleID, -CD8_T), y=CD8_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_specimen) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total Number of CD8 T cells")
CD8_count_persample


##..Plot and compare theCD4 and CD8 freq of total immune..##

# order by descending CD4:CD8 ratio
FC_order<-c(93,96,97,92,91,95,90,94,15,27,55,6,33,34,26,61,40,7,98,47,88,42,84,99,21,28,54,14,89,48)
freq_data_TB<-droplevels(freq_data[freq_data$samp_num %in% c(Stanford,AHRI,Autopsy),])
freq_data_TB$Origin<-tcell_data_TB$Origin
freq_data_TB$SampleID<-factor(freq_data_TB$SampleID, levels = FC_order)
freq_data_TB<-freq_data_TB[order(freq_data_TB$SampleID),]

# CD4 frequency per sample
CD4_freq_persample<-ggplot(data=freq_data_TB, aes(x=as.factor(SampleID), y=CD4_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_specimen) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total immune cells)")
CD4_freq_persample

# CD8 frequency per sample
CD8_freq_persample<-ggplot(data=freq_data_TB, aes(x=as.factor(SampleID), y=CD8_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_specimen) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total immune cells)")
CD8_freq_persample


##..Compare cell frequencies (T cells and others that appear related) by origin..##

my_comparisons = list(c('AHRI','Stanford'),
                   c('AHRI','Autopsy'),
                   c('Autopsy','Stanford'))

# CD4 T cell freq
data_summary <- freq_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD4_T))

CD4_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(Origin), y=CD4_T, color=Origin)) +
  geom_point(data=freq_data_TB, aes(x=as.factor(Origin), y=CD4_T, color=as.factor(Origin), shape = organ), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_shape_manual(values=c(17, 16)) +
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total immune cells)")
CD4_freq_bulk

# CD8 T cell freq
data_summary <- freq_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD8_T))


CD8_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(Origin), y=CD8_T, color=Origin)) +
  geom_point(data=freq_data_TB, aes(x=as.factor(Origin), y=CD8_T, color=as.factor(Origin), shape= organ), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_shape_manual(values=c(17, 16)) +
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total immune cells)")
CD8_freq_bulk


# CD14 monos freq
data_summary <- freq_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD14_Mono))


CD14_Mono_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(Origin), y=CD14_Mono, color=Origin)) +
  geom_point(data=freq_data_TB, aes(x=as.factor(Origin), y=CD14_Mono, color=as.factor(Origin), shape = organ), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_shape_manual(values=c(17, 16)) +
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD14+ Monos (of total immune cells)")
CD14_Mono_freq_bulk

# MDSC-like macrophage freq
data_summary <- freq_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(`CD11b/c_CD206_Mac/Mono`))

MDSC_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(Origin), y=`CD11b/c_CD206_Mac/Mono`, color=Origin)) +
  geom_point(data=freq_data_TB, aes(x=as.factor(Origin), y=`CD11b/c_CD206_Mac/Mono`, color=as.factor(Origin), shape = organ), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_shape_manual(values=c(17, 16)) +
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of MDSC (of total immune cells)")
MDSC_freq_bulk

##..Format the cell count data..##
count_data_TB<-droplevels(count_data[count_data$samp_num %in% c(Stanford,AHRI,Autopsy),])
count_data_TB$Origin<-tcell_data_TB$Origin
count_data_TB$SampleID<-factor(count_data_TB$SampleID, levels = FC_order)
count_data_TB<-count_data_TB[order(count_data_TB$SampleID),]

##..Compare counts across specimen origin for same cell types as above..##

# CD4 T cell count
data_summary <- count_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD4_T))

CD4_count_bulk<-ggplot(data=count_data_TB, aes(x=as.factor(Origin), y=CD4_T, color=Origin)) +
  geom_point(data=count_data_TB, aes(x=as.factor(Origin), y=CD4_T, color=as.factor(Origin)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total immune cells)")
CD4_count_bulk

# CD8 T cell count
data_summary <- count_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD8_T))

CD8_count_bulk<-ggplot(data=count_data_TB, aes(x=as.factor(Origin), y=CD8_T, color=Origin)) +
  geom_point(data=count_data_TB, aes(x=as.factor(Origin), y=CD8_T, color=as.factor(Origin)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total immune cells)")
CD8_count_bulk


# CD14 mono count
data_summary <- count_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(CD14_Mono))


CD14_Mono_count_bulk<-ggplot(data=count_data_TB, aes(x=as.factor(Origin), y=CD14_Mono, color=Origin)) +
  geom_point(data=count_data_TB, aes(x=as.factor(Origin), y=CD14_Mono, color=as.factor(Origin)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD14+ Monos (of total immune cells)")
CD14_Mono_count_bulk

# MDSC-like macrophage count
data_summary <- count_data_TB %>%
  group_by(Origin) %>%
  summarize(combo_median = median(`CD11b/c_CD206_Mac/Mono`))

MDSC_count_bulk<-ggplot(data=count_data_TB, aes(x=as.factor(Origin), y=`CD11b/c_CD206_Mac/Mono`, color=Origin)) +
  geom_point(data=count_data_TB, aes(x=as.factor(Origin), y=`CD11b/c_CD206_Mac/Mono`, color=as.factor(Origin)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_specimen) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  stat_compare_means(label = "p.signif", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of MDSC (of total immune cells)")
MDSC_count_bulk

##..Correlate MDSC:Mono ratio with CD4:CD8 ratio..##

# count
keep_cols<-c('SampleID','CD11b/c_CD206_Mac/Mono','CD14_Mono','CD4_T','CD8_T','samp_num','Origin','organ')
count_data_corr<-count_data_TB[,c(1,3,5,10,12,19,20,21)]
count_data_corr$tcell_ratio<-log2(count_data_corr$CD4_T/count_data_corr$CD8_T)
count_data_corr$mye_ratio<-log2(count_data_corr$`CD11b/c_CD206_Mac/Mono`/count_data_corr$CD14_Mono)

data_noInf<-do.call(data.frame,lapply(count_data_corr, function(x) replace(x, is.infinite(x), NA)))
summary(lm(mye_ratio~tcell_ratio, data_noInf))

cor.test(data_noInf$mye_ratio,data_noInf$tcell_ratio,method="pearson")

ggplot(data_noInf,aes(x=tcell_ratio,y=mye_ratio, color = Origin)) +
  geom_smooth(method='lm', formula= y~x, aes(group=1)) +
  geom_point(aes(group=1, size = 1, shape = organ)) + 
  scale_shape_manual(values=c(17, 16)) +
  scale_color_manual(values=color_specimen) +
  labs(x="T cell ratio (log2 [CD4 T cells / CD8 T cells])") + 
  labs(y="Myeloid Ratio (log2 [11b/c Mac / CD14 Mono])") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

##..Save dataframe SampleID, Samp_num, CD4 total, CD8 total, CD4/CD8 ratio and save..##

write.csv(tcell_data, file="CD4CD8ratio.csv",row.names = FALSE)

