# MIBI_evaluate_metabolic_zone_area.R
# Created by: Erin McCaffrey 
# Date created: 230726
#
# Overview: This script reads quantifies and visualizes the proportion of each 
# granuloma's myeloid core made up of the hypoxic and normoxic zones. It visualizes
# the relationship between these attributes and certain metadata. 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggstatsplot)
library(gtools)

##..Step 1: Import data..##

data<-read.csv('cell_cohort_data_metabolic_zones.csv')
area_data <- read.csv('./masks/all_samples_mask_data.csv')
signal_data <- read.csv('all_sample_metabolic_signal.csv')
animal_color_key <- read.csv("./keys/animal_color_key.csv")

##..Step 2: Simplify to sample x area..##

data_per_sample <- data[row.names(unique(data[,c("sample", "glyco_area","IDO1_area")])), 
                        names(data) %in% c("sample", "glyco_area","IDO1_area")]

##..Step 3: Append area data..##

data_per_sample <- left_join(data_per_sample, area_data, by = c('sample'))

##..Step 4: Generate ratio and frequency of total..##

data_per_sample$glyco_ratio <- data_per_sample$glyco_area / data_per_sample$IDO1_area
data_per_sample$glyco_freq <- data_per_sample$glyco_area / (data_per_sample$IDO1_area + data_per_sample$glyco_area)
data_per_sample$IDO_freq <- data_per_sample$IDO1_area / (data_per_sample$IDO1_area + data_per_sample$glyco_area)

data_per_sample$glyco_freq_total <- data_per_sample$glyco_area / (data_per_sample$cellular_px + data_per_sample$necrosis_px)
data_per_sample$IDO_freq_total <- data_per_sample$IDO1_area / (data_per_sample$cellular_px + data_per_sample$necrosis_px)
data_per_sample$total_ratio <- data_per_sample$glyco_freq_total / data_per_sample$IDO_freq_total

##..Step 5: Append granuloma metadata..##

meta_data <- read.csv('./cohort_metadata/study_cohort_metadata.csv')
data_anno <- left_join(data_per_sample, meta_data, by = c('sample'))

##..Step 6: Append signal data..##

data_anno <- left_join(data_anno, signal_data, by = c('sample'))

##..Step 7: Create signal ratio..##

data_anno$signal_ratio <- data_anno$glyco_area / data_anno$IDO1_area

##..Step 8: Plot and compare..##

id_cols <- c('sample', 'log_CFU','burden','Animal_Gran_Code','Animal')
measure_cols <- c('glyco_area','IDO1_area','glyco_ratio','glyco_freq',
                  'IDO_freq','FDG_SUV','glyco_freq_total','IDO_freq_total',
                  'total_ratio','glycozone_signal','IDO1_zone_signal', 
                  'signal_ratio')
data_melted <- melt(data_anno, id.vars = id_cols, measure.vars = measure_cols)

# define plotting data
plot_data <- data_melted

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

# plot correlation between metabolic zones and necrosis area
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

