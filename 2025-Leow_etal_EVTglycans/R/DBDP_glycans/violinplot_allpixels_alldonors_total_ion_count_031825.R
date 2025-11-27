# Get ratio of each glycan
# Author: Ke Leow
# Date: 03/18/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggpubr) 

#--------------------------------
# Load MALDI data
#--------------------------------
#data before scaling/batch correction
data = read_csv("data/DBDP_transition/MALDI/library_matched/DBDP_transition_all_pixels.csv")
#metadata
meta <- read_csv("data/DBDP_transition/MALDI/library_matched/fiji/Overlay Elements of H5N4F1_annotated.csv") %>% 
  separate(Name, into = c("tissue","week","day", "id"), sep = "_", remove = FALSE) %>% 
  unite("patient", tissue, week, day, sep = "_")

data <- inner_join(data[2:ncol(data)], meta %>% select(Name, Region, patient, id), by = c('mask' = 'Name')) 

#--------------------------------
# Calculate total ion count and plot
#--------------------------------
# Sum glycans for each sample
data$glycan_sum <- rowSums(data[, !(names(data) %in% c("mask", "Region","patient","id"))])

#get unique patients
patients = unique(meta$patient)

# Perform Wilcoxon test - per patient
wilcox.test(glycan_sum ~ Region, data = data %>% filter(patient == patients[2])) #"DBDP_8_2"
wilcox.test(glycan_sum ~ Region, data = data %>% filter(patient == patients[3])) #"DBDP_10_0"

# Create scatter plot - for each donor
patient_sub = patients[1]
data_sub <- data%>% filter(patient == patient_sub)

pdf(paste("R_plots/DBDP_transition/vioplot_totalioncount_", patient_sub, ".pdf", sep = ""), width = 3, height = 3)
ggplot(data_sub, aes(x = Region, y = glycan_sum)) +
  geom_violin(fill = "darkgrey") +
  # geom_boxplot(width=0.02, outlier.shape = NA, fill="white", col = "black")+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x=1.5) +  # Adds significance
  labs(x = NULL,
       y = "Total ion count",
       title = paste("Donor:", patient_sub, sep = " ")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
dev.off()

# Create scatter plot - for each donor
patient_sub = patients[2]
data_sub <- data%>% filter(patient == patient_sub)

pdf(paste("R_plots/DBDP_transition/vioplot_totalioncount_", patient_sub, ".pdf", sep = ""), width = 3, height = 3)
ggplot(data_sub, aes(x = Region, y = glycan_sum)) +
  geom_violin(fill = "darkgrey") +
  # geom_boxplot(width=0.02, outlier.shape = NA, fill="white", col = "black")+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x=1.5) +  # Adds significance
  labs(x = NULL,
       y = "Total ion count",
       title = paste("Donor:", patient_sub, sep = " ")
       ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
dev.off()

# Create scatter plot - for each donor
patient_sub = patients[3]
data_sub <- data%>% filter(patient == patient_sub)

pdf(paste("R_plots/DBDP_transition/vioplot_totalioncount_", patient_sub, ".pdf", sep = ""), width = 3, height = 3)
ggplot(data_sub, aes(x = Region, y = glycan_sum)) +
  geom_violin(fill = "darkgrey") +
  # geom_boxplot(width=0.02, outlier.shape = NA, fill="white", col = "black")+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x=1.5) +  # Adds significance
  labs(x = NULL,
       y = "Total ion count",
       title = paste("Donor:", patient_sub, sep = " ")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
dev.off()