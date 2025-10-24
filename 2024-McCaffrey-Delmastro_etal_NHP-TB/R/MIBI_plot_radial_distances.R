# MIBI_plot_radial_distances.R
# Date created: 9/19/2023
# Created by: Erin McCaffrey
# 
# Overview: For a selected sample this script plots the radial distribution of 
# selected markers.

library(Hmisc)
library("MetaCyto")
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggalt)
library(ggridges)
library(ggstatsplot)
library(ggpubr)
library(cowplot)
library(introdataviz)

##..Step 1: Import data..##

necrosis_function_data <- read.csv('./spatial_analysis/radial/necrosis_marker_distance_data_all-grans.csv')

##..Step 2: Pre-process..##

# optionally subset by sample
sample_subset <- c('sample46')
data_sub <- necrosis_function_data[necrosis_function_data$sample %in% sample_subset, ]

# normalize the distance to the max distance per sample
data_sub <- data_sub %>%
  group_by(sample) %>%
  mutate(norm_necrosis_distance = necrosis_distance / max(necrosis_distance, na.rm = TRUE)) %>%
  ungroup()

# normalize the signal to 99.9% percentile
max_IDO1 <- quantile(data_sub$IDO1, probs=c(0.999))
max_GLUT1 <- quantile(data_sub$GLUT1, probs=c(0.999))

data_sub$GLUT1 <- data_sub$GLUT1/max_GLUT1
data_sub$IDO1 <- data_sub$IDO1/max_IDO1

# optionally subset further (percentage of each sample)
data_sub_subset <- data_sub %>%
  group_by(sample) %>%
  slice_sample(prop = 0.25) %>%   # keep 25% of rows per sample
  ungroup()

# melt for plotting
data_melted <-melt(data_sub, id.vars = c('sample','norm_necrosis_distance'), measure.vars = c('GLUT1','IDO1'))

##..Step 3: Plot..##

track_plot <- ggplot(data_melted, aes(x=norm_necrosis_distance, 
                                       y = value, color = variable)) + 
  geom_smooth(se=TRUE, linetype="dashed", method = stats::loess) +
  theme_bw() +
  theme(legend.position = 'none') + 
  xlim(0,1) 
ggsave("test.png", track_plot)
track_plot

histo_plot <- ggplot(data_melted_2,aes(x = value, color = variable)) + 
  geom_histogram() +
  theme_bw()
histo_plot
