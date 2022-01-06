#TBMIBIexamineMETissueType.R
#Date created: 03/05/20
#Author: Erin McCaffrey
#This script reads in the ME data and assesses freqeucny of MEs across three groups:
# 1. Lung resection 2. lung biopsy 3. extrapulmonary biopsy

library(ggpubr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(forcats)

##...Load in data..##

# ME data

ME_data<-read.csv('data/allTB-ME_annotated.csv')

# Define categories

# Plot individual examples and add origin data

resection<-c(21,84,42,88,28,89)
diagnostic_pulm<-c(14,15,98,99)
diagnostic_expulm<-c(6,7,33,34,26,27,40,61,47,48,54,55,92,93)
autopsy<-c(90,91,94,95,96,97)

# Get ME frequencies

freq_MEs <- as.data.frame(table(ME_data$SampleID, ME_data$maxME))
colnames(freq_MEs) <- c('SampleID', 'ME', 'Count')
sample_totals<-as.data.frame(table(ME_data$SampleID))
freq_MEs$Total <- rep(sample_totals$Freq, 8)
freq_MEs$ME_freq <- freq_MEs$Count / freq_MEs$Total
freq_MEs$Tissue<-'lung_resection'
freq_MEs[freq_MEs$SampleID %in% diagnostic_pulm,]$Tissue<-'lung_biopsy'
freq_MEs[freq_MEs$SampleID %in% diagnostic_expulm,]$Tissue<-'extrapulm_biopsy'
freq_MEs[freq_MEs$SampleID %in% autopsy,]$Tissue<-'lung_postmortem'

my_comparisons <- list(c("extrapulm_biopsy","lung_resection"),
                       c("extrapulm_biopsy","lung_biopsy"),
                       c("lung_biopsy","lung_resection"),
                       c("lung_biopsy","lung_postmortem"),
                       c("lung_postmortem","lung_resection"),
                       c("extrapulm_biopsy","lung_postmortem"))


##..Produce plots with median line and dots..## 

data_summary <- freq_MEs %>%
  group_by(Tissue,ME) %>%
  summarize(combo_median = median(ME_freq),
            combo_se = sqrt(var(ME_freq)/length(ME_freq)))


ggplot(data = freq_MEs, aes(x = Tissue, y = ME_freq)) + 
  geom_point(aes(x = Tissue, y = ME_freq), 
             position = position_jitter(width = 0.1, height = 0.0),size = 1) + 
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "red", width = 0.5, data = data_summary) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  facet_wrap(~ME, scale='free_y', ncol=4)
