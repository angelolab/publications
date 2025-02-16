# MIBI_plot_cell_data.R
# Author: Erin McCaffrey 
# Overview: This script generates plots related to cell frequency and counts
# across all granulomas. 

require(dplyr)
library(reshape2)
library(gtools)
library(ggplot2)

##..Import data..##
setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
data <- read.csv('cell_stats_all_samples_meta_data.csv')
data <- droplevels(data[data$category == 'pheno_of_total',])
data <- tibble::rowid_to_column(data, "ID")
metadata <- read.csv("./cohort_metadata/study_cohort_metadata.csv")
animal_color_key <- read.csv("./keys/animal_color_key.csv")
cell_color_key <- read.csv("./keys/cell_color_key.csv")

##..Generate a plot of the total counts per animal..##

# generate sample by total counts
total_counts <- data %>%
  group_by(sample) %>%
  summarize(total_cells = sum(n))

# append the cellular area
cellular_area <- unique(data[,names(data) %in% c('sample','cellular_px')])
total_counts_area <- left_join(total_counts, cellular_area, by = c('sample'))

# generate density 
total_counts_area$density <- total_counts_area$total_cells / total_counts_area$cellular_px

# append animal ID
total_counts_area <- left_join(total_counts_area, 
                               metadata[,names(metadata) %in% c('sample',
                                                                'Animal_Code',
                                                                'Animal_Gran_Code')], 
                               by = c('sample'))

# sort the animals in natural order
plot_animals <- levels(factor(total_counts_area$Animal_Code))
plot_animals <- mixedsort(plot_animals)
total_counts_area$Animal_Code <- factor(total_counts_area$Animal_Code, levels = plot_animals)
total_counts_area <- total_counts_area[order(total_counts_area$Animal_Code),]

# get the plot colors
plot_animals<-levels(factor(total_counts_area$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal_Code<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

# plot
total_cells <- ggplot(total_counts_area, aes(x=reorder(as.factor(Animal_Gran_Code), -total_cells), 
                                y=total_cells, fill = as.factor(Animal_Code))) + 
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values = color) +
  labs(y="Total Cell Count", x="Granuloma") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
total_cells

total_density <- ggplot(total_counts_area, aes(x=reorder(as.factor(Animal_Gran_Code), -density), 
                                             y=density, fill = as.factor(Animal_Code))) + 
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values = color) +
  labs(y="Total Cell Count", x="Granuloma") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
total_density

##..Summarize individual cell frequencies..##

# generate count by cell type
total_subset_counts <- data %>%
  group_by(variable) %>%
  summarize(total = sum(n))

# normalize to total
cell_total <- as.numeric(sum(total_subset_counts$total))
total_subset_counts$frequency <- total_subset_counts$total / cell_total

# sort by descending frequency
total_subset_counts <- total_subset_counts[order(-total_subset_counts$frequency),] 
cell_order <- total_subset_counts$variable
total_subset_counts$variable<- factor(total_subset_counts$variable, levels = cell_order)

# get color key 
plot_populations <- levels(factor(total_subset_counts$variable))
plot_colors <- droplevels(cell_color_key[cell_color_key$Pheno %in% plot_populations,])
plot_colors$Pheno <- factor(plot_colors$Pheno, levels = plot_populations)
plot_colors <- plot_colors[order(plot_colors$Pheno),]
color <- as.vector(plot_colors$Hex)

# plot pie chart
cell_pie <- ggplot(total_subset_counts, aes(x = "", y = frequency, fill = variable)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = color) +
  theme_void()
cell_pie

# add animal granuloma code to the dataset
data <- left_join(data, metadata[,names(metadata) %in% c('sample',
                                                         'Animal_Code',
                                                         'Animal_Gran_Code')],
                  by = c('sample'))

# sort the cell types per the order in the pie chart
data$variable <- factor(data$variable, levels=rev(cell_order))
data <- data[order(data$variable),]

# sort the animal-granuloma codes
data$Animal_Gran_Code <- sub("-", ".", data$Animal_Gran_Code, fixed=TRUE)
plot_grans <- levels(factor(data$Animal_Gran_Code))
plot_grans <- mixedsort(plot_grans, scientific = FALSE)
data$Animal_Gran_Code <- factor(data$Animal_Gran_Code, levels = plot_grans)

# get the color key
plot_populations <- levels(factor(data$variable))
plot_colors <- droplevels(cell_color_key[cell_color_key$Pheno %in% plot_populations,])
plot_colors$Pheno <- factor(plot_colors$Pheno, levels = plot_populations)
plot_colors <- plot_colors[order(plot_colors$Pheno),]
color <- as.vector(plot_colors$Hex)

# generate stacked bar
cell_bar <- ggplot(data, aes(x=Animal_Gran_Code, y=freq_of_total, fill=variable)) + 
  theme_bw() +
  scale_fill_manual(values=color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  theme(axis.text.x=element_text(angle=90, hjust=1))
cell_bar

