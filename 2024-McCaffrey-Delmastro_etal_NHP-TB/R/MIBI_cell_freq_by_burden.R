# MIBI_cell_freq_analysis.R
# Date created: 12/21/2023
# Author: Erin McCaffrey
#  
# Based on the correlation between CFU and the abudnace of CD11c+ Macrophages, 
# here we perform some additional analysis of the macrophage populations found
# in the TB granulomas

library(forcats)
library(viridis)
library(dplyr)
library(stringr)
library(ggpubr)
devtools::install_github("psyteachr/introdataviz")
library(introdataviz)

## Read in data ##
setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
data<-read.csv('cell_stats_all_samples_meta_data.csv')
data<-droplevels(data[data$category == 'pheno_of_total',])
data<-tibble::rowid_to_column(data, "ID")

## Create quartiles for CFU ##
data <- data %>% mutate(CFU_q = cut(log_CFU, quantile(log_CFU), include.lowest=TRUE, labels=FALSE))

## Separate count and density data ##
freq_data<-reshape2::dcast(data, sample + log_CFU + burden + CFU_q ~ variable, value.var = "freq_of_total", fun.aggregate = sum)
density_data<-reshape2::dcast(data, sample + log_CFU + burden + CFU_q ~ variable, value.var = "cell_density", fun.aggregate = sum)

## Melt for easy plotting ##
freq_data.m <- reshape2::melt(freq_data, id.vars = c('sample', 'burden', 'log_CFU','CFU_q'))
density_data.m <- reshape2::melt(density_data, id.vars = c('sample', 'burden', 'log_CFU','CFU_q'))


## Plot ##
plot_data <- freq_data.m
ggplot(data = plot_data, aes(x = fct_reorder(as.factor(variable), value,.fun=median,.desc=TRUE), 
                             value, fill=burden)) +
  geom_boxplot(width = 0.75, alpha = .6, fatten = NULL, show.legend = TRUE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.75)) +
  stat_compare_means(aes(group = burden), method= "wilcox.test", label = "p.format") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total Macrophages') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

## Plot between quartiles
plot_data <- freq_data.m
plot_data <- plot_data[plot_data$CFU_q %in% c(1,4),]
ggplot(data = plot_data, aes(x = fct_reorder(as.factor(variable), value,.fun=median,.desc=TRUE), 
                             value, fill=as.factor(CFU_q))) +
  geom_boxplot(width = 0.75, alpha = .6, fatten = NULL, show.legend = TRUE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.75)) +
  stat_compare_means(aes(group = burden), method= "wilcox.test", label = "p.format") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

