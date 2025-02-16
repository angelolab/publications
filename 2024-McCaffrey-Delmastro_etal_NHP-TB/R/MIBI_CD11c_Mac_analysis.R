# MIBI_macrophage_composition_analysis.R
# Author: Erin McCaffrey
#  
# Overview: Based on the correlation between CFU and the abudnace of CD11c+ Macrophages, 
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
# data <- data %>% mutate(CFU_q = cut(log_CFU, quantile(log_CFU, probs = c(0, 1/3, 2/3, 1)), include.lowest=TRUE, labels=FALSE))
data <- data %>% mutate(CFU_q = cut(log_CFU, quantile(log_CFU), include.lowest=TRUE, labels=FALSE))

## Subset to only include macrophages ##
macs <- c("CD11c+_Mac","CD14+CD11c+_Mac","CD14+_Mac_Mono","CD163+_Mac","CD206+_Mac","CD68+_Mac",
          "FN1+_Mac",'giant_cell')
data_mac <- droplevels(data[data$variable %in% macs,])

## Separate count and density data ##
count_data<-reshape2::dcast(data_mac, sample + log_CFU + burden + CFU_q ~ variable, value.var = "n", fun.aggregate = sum)
density_data<-reshape2::dcast(data_mac, sample + log_CFU + burden + CFU_q ~ variable, value.var = "cell_density", fun.aggregate = sum)

## Append total macrophages ##
count_data <- count_data %>%
  mutate(total_macs = rowSums(select(.,5:12)))

freq_data <- count_data %>% 
  mutate_at(vars(5:12),function(i)i/.$total_macs)

## Melt for easy plotting ##
freq_data.m <- reshape2::melt(freq_data, id.vars = c('sample', 'burden', 'log_CFU','CFU_q'))
density_data.m <- reshape2::melt(density_data, id.vars = c('sample', 'burden', 'log_CFU','CFU_q'))

## Plot ##
plot_data <- freq_data.m[!freq_data.m$variable %in% c('total_macs','gran_CFU'),]
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

## Plot between quartiles ##
plot_data <- freq_data.m[!freq_data.m$variable %in% c('total_macs','gran_CFU'),]
plot_data <- plot_data[plot_data$CFU_q %in% c(1,3),]
ggplot(data = plot_data, aes(x = fct_reorder(as.factor(variable), value,.fun=median,.desc=TRUE), 
                             value, fill=as.factor(CFU_q))) +
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

##..Generate log2(FC) between the high versus low burden grans..##
summary_data <- freq_data[,!names(freq_data) %in% c('sample',
                                                    'log_CFU',
                                                    'CFU_q',
                                                    'total_macs')]
data_summary <- summary_data %>%
  group_by(burden) %>%
  summarize_all(list(mean=mean))

data_summary<-as.data.frame(t(data_summary[,-c(1)]))
names(data_summary)<-c('high','low')
data_summary$pheno_corrected <- rownames(data_summary)
data_summary$pheno_corrected <- stringr::str_remove(data_summary$pheno_corrected, "_mean")
data_summary$high<-as.numeric(as.character(data_summary$high))
data_summary$low<-as.numeric(as.character(data_summary$low))
data_summary$FC<-log2(data_summary$high/data_summary$low)

color_key <- read.csv("./keys/cell_color_key.csv")
plot_populations<-levels(factor(data_summary$pheno_corrected))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

ggplot(data=data_summary[is.finite(data_summary$FC),], 
       aes(x=reorder(pheno_corrected, -FC), y=FC, fill=pheno_corrected)) +
  geom_bar(stat="Identity", width = 1) +
  scale_fill_manual(values = color) +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_minimal() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  labs(y="log2(mean freq high / mean freq low)") 


