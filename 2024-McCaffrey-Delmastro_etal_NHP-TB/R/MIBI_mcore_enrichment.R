# MIBI_metabolic_zone_enrichment.R

library(dplyr)
library(reshape2)
library(ggplot2)
library(forcats)

##..Import data..##

setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
data <- read.csv('cell_cohort_data_metabolic_zones.csv')

##..Drop unassigned cells..##
data <- droplevels(data [!data $pheno_corrected == 'Unassigned',])

##..Preprocess..##

# create a column for myeloid core and not myeloid core
data$mcore <- 0
data[data$glyco_zone == 1 | data$IDO1_zone ==1, ]$mcore <- 1
mcore_cells <- data[data$mcore == 1, ]

##..Get enrichment score..##

# total number of cells
totals_overall <- as.data.frame(table(data$sample, data$pheno_corrected))
colnames(totals_overall) <- c('sample','pheno_corrected','total_overall')

# total number mcore
totals_mcore <- as.data.frame(table(mcore_cells$sample, mcore_cells$pheno_corrected))
colnames(totals_mcore) <- c('sample','pheno_corrected','mcore_total')

# merge
all_data <- merge(totals_overall, totals_mcore, by = c('sample', 'pheno_corrected'), all = T)

# get non-mcore counts
all_data$non_mcore_total <- all_data$total_overall - all_data$mcore_total

# get frequencies
all_data$freq_mcore <- all_data$mcore_total / all_data$total_overall
all_data$freq_non_mcore <- all_data$non_mcore_total / all_data$total_overall

# append the areas
area_data <- read.csv('./masks/all_samples_mask_data.csv')
metabolic_data<-read.csv('cell_cohort_data_metabolic_zones.csv')
metabolic_data_persample <- metabolic_data[row.names(unique(metabolic_data[,c("sample", "glyco_area","IDO1_area")])), 
                                           names(metabolic_data) %in% c("sample", "glyco_area","IDO1_area")]
all_data_area <- merge(all_data, area_data, by = c('sample'), all = T)
all_data_area <- merge(all_data_area, metabolic_data_persample, by = c('sample'), all = T)
all_data_area$mcore_area <- all_data_area$glyco_area + all_data_area$IDO1_area
all_data_area$non_mcore_area <- all_data_area$cellular_px - all_data_area$mcore_area

# get the cell densities
all_data_area$mcore_density <- all_data_area$mcore_total / all_data_area$mcore_area
all_data_area$non_mcore_density <- all_data_area$non_mcore_total/ all_data_area$non_mcore_area

# get the freq densities
all_data_area$mcore_freq_density <- all_data_area$freq_mcore / all_data_area$mcore_area
all_data_area$non_mcore_freq_density <- all_data_area$freq_non_mcore / all_data_area$non_mcore_area

# replace NA with 0
all_data_area[is.na(all_data_area)] <- 0

# calculate enrichment
all_data_area$enrichment <- log2(all_data_area$mcore_density / all_data_area$non_mcore_density)
all_data_area$freq_enrichment <- log2(all_data_area$mcore_freq_density / all_data_area$non_mcore_freq_density)

# remove Inf rows
plot_data <- all_data_area[is.finite(all_data_area$enrichment),]

##..Visualize enrichment score..##

# all cells
color_key <- read.csv("./keys/cell_color_key.csv")
plot_populations<-levels(factor(plot_data$pheno_corrected))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

ggplot(plot_data, aes(fct_reorder(pheno_corrected, enrichment,
                                      .fun=median,.desc=TRUE), enrichment, 
                      fill = pheno_corrected)) +
  geom_boxplot() +
  scale_fill_manual(values = color) +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() +
  labs(x = 'Cell type') + 
  labs(y = 'Enrichment') +
  theme(axis.text.x = element_text(angle=90,hjust=1)) +
  theme(legend.position = "none") 

# just macrophages
macs <- c("CD11c+_Mac","CD14+CD11c+_Mac","CD14+_Mac_Mono","CD163+_Mac","CD206+_Mac","CD68+_Mac",
          "FN1+_Mac", "giant_cell")
mac_data <- plot_data[plot_data$pheno_corrected %in% macs, ]

plot_populations<-levels(factor(mac_data$pheno_corrected))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

ggplot(mac_data, aes(fct_reorder(pheno_corrected, enrichment,
                                  .fun=median,.desc=TRUE), enrichment, 
                      fill = pheno_corrected)) +
  geom_boxplot() +
  scale_fill_manual(values = color) +
  geom_hline(yintercept=c(median(mac_data$enrichment)), linetype="dashed") +
  theme_minimal() +
  theme(legend.position = "none") 




