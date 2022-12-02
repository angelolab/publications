# PAH_MIBIevalVesselArea.R
# Author: Erin McCaffrey 
# Date created: 200327
# Script reads in vessel size data for all PAH and healthy FOVs. Visualizes data
# as histograms and violin plots.


library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data<-read.csv("vessel_area_data.csv")


##..Append condition info..##
annotated_data<-read.csv("celldata_region_annotated.csv")

ipah<-unique(annotated_data[annotated_data$Subgroup=='IPAH',]$Point_num)
hpah<-unique(annotated_data[annotated_data$Subgroup=='HPAH',]$Point_num)
hlt<-unique(annotated_data[annotated_data$Subgroup=='Healthy Controls',]$Point_num)
data$group<-'hlt'
data[data$Point_num %in% hpah, ]$group <- 'hpah'
data[data$Point_num %in% ipah, ]$group <- 'ipah'

##..Plot Histogram..##

ggplot(data, aes(area, fill = group)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity', 
                 center = 5000, bins = 50)

##..Plot Violin..##

my_comparisons<-list(c("hlt","hpah"),
                     c("hlt","ipah"),
                     c("ipah","hpah"))

ggplot(data, aes(x = group, y = area, colour = group)) + 
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method= "wilcox.test") +
  geom_violin() + 
  geom_point(position = "jitter", aes(colour = group)) + 
  theme_bw() +
  scale_color_manual(values=c('#03CD03','#0700F5','#E950F5'))

##..Sum within each point..##

data_summary <- data %>%
  group_by(Point_num) %>%
  summarize(sum = sum(area),
            mean = mean(area),
            median = median(area),
            max = max(area),
            min = min(area))

data_summary$group<-'hlt'
data_summary[data_summary$Point_num %in% hpah, ]$group <- 'hpah'
data_summary[data_summary$Point_num %in% ipah, ]$group <- 'ipah'


##..Append PID..##

key<-annotated_data[,c('PID','Point_num')]
key_unique<-key[!duplicated(t(apply(key, 1, sort))),]
PID_data<-merge(data, key_unique, by.x="Point_num", by.y="Point_num")

##..Summarize within each patient..##

data_summary_PID <- PID_data %>%
  group_by(PID) %>%
  summarize(sum = sum(area),
            mean = mean(area),
            median = median(area),
            max = max(area),
            min = min(area))

ipah_PID<-unique(annotated_data[annotated_data$Subgroup=='IPAH',]$Point_num)
hpah_PID<-unique(annotated_data[annotated_data$Subgroup=='HPAH',]$Point_num)
hlt_PID<-unique(annotated_data[annotated_data$Subgroup=='Healthy Controls',]$Point_num)

data_summary_PID$group<-'hlt'
data_summary_PID[data_summary_PID$PID %in% hpah_PID, ]$group <- 'hpah'
data_summary_PID[data_summary_PID$PID %in% ipah_PID, ]$group <- 'ipah'


##..Plot histogram..##

ggplot(data_summary, aes(sum, fill = group)) + 
  scale_fill_manual(values=c('#03CD03','#0700F5','#E950F5')) +
  theme(legend.position = "none") + 
  geom_density(alpha = 0.5)

##..Plot violin with stats..##
  
ggplot(data_summary, aes(x = group, y = sum, colour = group)) + 
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method= "wilcox.test") +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.1, height = 0.0),size = 3, aes(colour = group)) + 
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('#03CD03','#0700F5','#E950F5'))



