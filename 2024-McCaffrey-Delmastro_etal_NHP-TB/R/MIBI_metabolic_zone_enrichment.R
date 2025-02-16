# MIBI_mcore_enrichment.R

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

# for now drop the double-zone cells, will run the tiebreaker function later
double_pos_rows <- rownames(data[data$glyco_zone == 1 & data$IDO1_zone ==1, ])
double_neg_rows <- rownames(data[data$glyco_zone == 0 & data$IDO1_zone == 0, ])
# drop_rows <- c(double_pos_rows, double_neg_rows)
data[(row.names(data) %in% double_pos_rows),]$IDO1_zone <- 0
drop_rows <- c(double_neg_rows)
metabolic_cells <- data[!(row.names(data) %in% drop_rows), ]

##..Get enrichment score..##

# total number of cells
totals_overall <- as.data.frame(table(data$sample, data$pheno_corrected))
colnames(totals_overall) <- c('sample','pheno_corrected','total_overall')

# total number metabolic cells
totals <- as.data.frame(table(metabolic_cells$sample, metabolic_cells$pheno_corrected))
colnames(totals) <- c('sample','pheno_corrected','total')

# count of glyco cells
glyco_data <- droplevels(metabolic_cells[metabolic_cells$glyco_zone == 1, ]) 
glyco_totals <- as.data.frame(table(glyco_data$sample, glyco_data$pheno_corrected))
colnames(glyco_totals) <- c('sample','pheno_corrected','total_glyco')

# count of IDO1 cells
IDO1_data <- droplevels(metabolic_cells[metabolic_cells$IDO1_zone == 1, ]) 
IDO1_totals <- as.data.frame(table(IDO1_data$sample, IDO1_data$pheno_corrected))
colnames(IDO1_totals) <- c('sample','pheno_corrected','total_IDO1')

# merge
all_data <- merge(totals_overall, totals, by = c('sample', 'pheno_corrected'), all = T)
all_data <- merge(all_data, glyco_totals, by = c('sample', 'pheno_corrected'), all = T)
all_data <- merge(all_data, IDO1_totals, by = c('sample', 'pheno_corrected'), all = T)

# get frequencies
all_data$freq_glyco_overall <- all_data$total_glyco / all_data$total_overall
all_data$freq_IDO1_overall <- all_data$total_IDO1 / all_data$total_overall
all_data$freq_glyco_zone <- all_data$total_glyco / all_data$total
all_data$freq_IDO1_zone <- all_data$total_IDO1 / all_data$total

# append the areas
area_data <- metabolic_cells[!duplicated(metabolic_cells[7:8]), c(2,7:8)]
all_data_area <- merge(all_data, area_data, by = c('sample'), all = T)

# get the cell densities
all_data_area$glyco_density <- all_data_area$total_glyco / all_data_area$glyco_area
all_data_area$IDO1_density <- all_data_area$total_IDO1 / all_data_area$IDO1_area

# get the freq densities
all_data_area$glyco_freq_density <- all_data_area$freq_glyco_overall / all_data_area$glyco_area
all_data_area$IDO1_freq_density <- all_data_area$freq_IDO1_overall / all_data_area$IDO1_area
all_data_area$glyco_zone_freq_density <- all_data_area$freq_glyco_zone / all_data_area$glyco_area
all_data_area$IDO1_zone_freq_density <- all_data_area$freq_IDO1_zone / all_data_area$IDO1_area

# replace NA with 0
all_data_area[is.na(all_data_area)] <- 0

# calculate enrichment
all_data_area$enrichment <- log2(all_data_area$glyco_density / all_data_area$IDO1_density)
all_data_area$freq_enrichment <- log2(all_data_area$glyco_freq_density / all_data_area$IDO1_freq_density)

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
  theme(axis.text.x = element_text(angle=35,hjust=1)) +
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

# just the means or medians

###...Compare between the zones...##
metabolic_data.m <- reshape2::melt(all_data_area, 
                                   id.vars = c('sample', 'pheno_corrected'),
                                   measure.vars= c('glyco_density',
                                                   'IDO1_density',
                                                   'glyco_freq_density',
                                                   'IDO1_freq_density',
                                                   'glyco_zone_freq_density',
                                                   'IDO1_zone_freq_density'))

plot_data <- metabolic_data.m[metabolic_data.m$variable %in% c('glyco_zone_freq_density',
                                                               'IDO1_zone_freq_density'),]

ggplot(data = plot_data, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.75, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.75)) +
  stat_compare_means(method= "wilcox.test", label = "p.format") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Metabolic Zone') + 
  labs(y = 'Frequency in Zone') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) + 
  facet_wrap(~pheno_corrected, scales = c('free'))
