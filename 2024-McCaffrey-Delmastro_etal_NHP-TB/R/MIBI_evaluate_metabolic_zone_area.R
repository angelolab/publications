# evaluate_bulk_GLUT1.R
# Author: Erin McCaffrey 
# Date created: 230726
# This script reads in the bulk GLUT1 data produced by MIBI_quantify_bulk_GLUT1.py
# and assesses the relationship between various metrics with bacterial burden.

library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggstatsplot)
library(gtools)

##..Import data..##

setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
data<-read.csv('cell_cohort_data_metabolic_zones.csv')
area_data <- read.csv('./masks/all_samples_mask_data.csv')
signal_data <- read.csv('all_sample_metabolic_signal.csv')
animal_color_key <- read.csv("./keys/animal_color_key.csv")

##..Simplify to sample x area..##

data_per_sample <- data[row.names(unique(data[,c("sample", "glyco_area","IDO1_area")])), 
                        names(data) %in% c("sample", "glyco_area","IDO1_area")]

##..Append area data..##
data_per_sample <- left_join(data_per_sample, area_data, by = c('sample'))

##..Generate ratio and frequency of total..##
data_per_sample$glyco_ratio <- data_per_sample$glyco_area / data_per_sample$IDO1_area
data_per_sample$glyco_freq <- data_per_sample$glyco_area / (data_per_sample$IDO1_area + data_per_sample$glyco_area)
data_per_sample$IDO_freq <- data_per_sample$IDO1_area / (data_per_sample$IDO1_area + data_per_sample$glyco_area)

data_per_sample$glyco_freq_total <- data_per_sample$glyco_area / (data_per_sample$cellular_px + data_per_sample$necrosis_px)
data_per_sample$IDO_freq_total <- data_per_sample$IDO1_area / (data_per_sample$cellular_px + data_per_sample$necrosis_px)
data_per_sample$total_ratio <- data_per_sample$glyco_freq_total / data_per_sample$IDO_freq_total

##..Append granuloma metadata..##
meta_data <- read.csv('./cohort_metadata/study_cohort_metadata.csv')
data_anno <- left_join(data_per_sample, meta_data, by = c('sample'))

##..Append signal data..##
data_anno <- left_join(data_anno, signal_data, by = c('sample'))

##..Create signal ratio..##
data_anno$signal_ratio <- data_anno$glyco_area / data_anno$IDO1_area

##..Plot and compare..##
id_cols <- c('sample', 'log_CFU','burden','Animal_Gran_Code','Animal')
measure_cols <- c('glyco_area','IDO1_area','glyco_ratio','glyco_freq',
                  'IDO_freq','FDG_SUV','glyco_freq_total','IDO_freq_total',
                  'total_ratio','glycozone_signal','IDO1_zone_signal', 
                  'signal_ratio')
data_melted <- melt(data_anno, id.vars = id_cols, measure.vars = measure_cols)

# select channel
plot_data <- data_melted

# plot 1: value distributions 
ggplot(plot_data, aes(reorder(sample, log_CFU), value, fill=burden)) +
  geom_bar(stat="Identity", width = 1) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=90,hjust=1)) + 
  facet_wrap(~variable, scale = "free")

# plot 2: biaxial scatter
ggplot(plot_data, aes(log_CFU, value)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_minimal() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  facet_wrap(~variable, scale = "free")

# plot 3: variables by burden
ggplot(plot_data, aes(burden, value)) +
  geom_boxplot() +
  stat_compare_means(aes(group = burden), method= "t.test", label = "p.format") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~variable, scale = "free")

# plot 4: corr for selected parameter
ggscatterstats(
  data = plot_data[plot_data$variable == 'glyco_freq_total',],
  x    = log_CFU,
  y    = value,
  type = "pearson")

# plot 5: area proportion per granuloma

# get plot data
plot_data <- data_melted[data_melted$variable %in% c('glyco_freq','IDO_freq'),]

# sort the animal-granuloma codes
plot_data$Animal_Gran_Code <- sub("-", ".", plot_data$Animal_Gran_Code, fixed=TRUE)
plot_grans <- levels(factor(plot_data$Animal_Gran_Code))
plot_grans <- mixedsort(plot_grans, scientific = FALSE)
plot_data$Animal_Gran_Code <- factor(plot_data$Animal_Gran_Code, levels = plot_grans)

# plot
ggplot(plot_data, aes(x=Animal_Gran_Code, y=value, fill=variable)) + 
  theme_bw() +
  scale_fill_manual(values=c('#00F222','#0022F2')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("NHP.Granuloma") +
  ylab("Frequency") +
  theme(axis.text.x=element_text(angle=90, hjust=1))

# plot 6: correlation between metabolic zones and necrosis area
data_anno$necrosis_freq <- data_anno$necrosis_px / (data_anno$necrosis_px + data_anno$cellular_px)
data_anno$freq_ratio <- log2(data_anno$glyco_freq / data_anno$IDO_freq)

corr_data <- data_anno[is.finite(data_anno$freq_ratio),]

ggscatterstats(
  data = corr_data,
  x    = freq_ratio,
  y    = necrosis_freq,
  type = "spearman")

plot_animals<-levels(factor(corr_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal_Code<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(corr_data, aes(necrosis_freq, freq_ratio)) +
  geom_point(aes(color = Animal_Code), show.legend = FALSE) +
  scale_color_manual(values = color) +
  theme_bw() 

