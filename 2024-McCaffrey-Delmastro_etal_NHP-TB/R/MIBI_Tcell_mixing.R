# MIBI_Tcell_mixing.R
# Date created: 03/21/24
# This script takes the output of the mixing score analysis and evaluates the 
# relationship between mixing between cell types, across zones, and with CFU

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggstatsplot)
library(ggsignif)

##..Read in the data..##
setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2/")
total_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_whole_gran.csv')
mcore_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_mcore-only.csv')
glyco_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_glycozone-only.csv')
ido_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_IDOzone-only.csv')
metabolic_area_data <- read.csv('cell_cohort_data_metabolic_zones.csv')
area_data <- read.csv('./masks/all_samples_mask_data.csv')
meta_data <- read.csv('./cohort_metadata/study_cohort_metadata.csv')

##..First analyze the mixing scores not broken down by cell type..##
total_mixing_sub <- unique(total_mixing[,c(1,2)])
colnames(total_mixing_sub) <- c('sample','total_mixing')

mcore_mixing_sub <- unique(mcore_mixing[,c(1,2)])
colnames(mcore_mixing_sub) <- c('sample','mcore_mixing')

glyco_mixing_sub <- unique(glyco_mixing[,c(1,2)])
colnames(glyco_mixing_sub) <- c('sample','glyco_mixing')

ido_mixing_sub <- unique(ido_mixing[,c(1,2)])
colnames(ido_mixing_sub) <- c('sample','ido_mixing')

##..Merge the area and metadata..##
meta_sub <- meta_data[,c(1,8,9,10)]
mixing_merged <- merge(total_mixing_sub, mcore_mixing_sub) %>% 
  merge(glyco_mixing_sub) %>%
  merge(ido_mixing_sub) %>%
  merge(meta_sub)

##..Melt for easy plotting..##
plot_data <- mixing_merged[,!names(mixing_merged) %in% c('total_mixing','mcore_mixing')]
data.m <- melt(plot_data, id.vars=c('sample','gran_CFU','log_CFU','burden'))

##..Compare overall mixing..##
ggplot(data.m, aes(variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, color = NA, scale = "width") +
  theme_minimal() +
  scale_fill_manual(values=c('#00F422','#0022F4')) +
  stat_summary(fun = median, geom = "point",size=3) 
  # theme(legend.position = 'none')

compare_means(value ~ burden, data.m, method= "t.test")

##..Break down CD4 and CD8 T cells..##
CD4_glyco_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_glycozone-only_CD4T.csv')
CD8_glyco_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_glycozone-only_CD8T.csv')
CD4_IDO1_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_IDOzone-only_CD4T.csv')
CD8_IDO1_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_IDOzone-only_CD8T.csv')

##..Extract the total scores..##
CD4_glyco_sub <- unique(CD4_glyco_mixing[,c(1,2)])
colnames(CD4_glyco_sub) <- c('sample','CD4_glyco_mixing')

CD8_glyco_sub <- unique(CD8_glyco_mixing[,c(1,2)])
colnames(CD8_glyco_sub) <- c('sample','CD8_glyco_mixing')

CD4_IDO1_sub <- unique(CD4_IDO1_mixing[,c(1,2)])
colnames(CD4_IDO1_sub) <- c('sample','CD4_IDO_mixing')

CD8_IDO1_sub <- unique(CD8_IDO1_mixing[,c(1,2)])
colnames(CD8_IDO1_sub) <- c('sample','CD8_IDO_mixing')

##..Merge the area and metadata..##
subset_mixing_merged <- merge(CD4_glyco_sub, CD8_glyco_sub) %>% 
  merge(CD4_IDO1_sub) %>%
  merge(CD8_IDO1_sub) %>%
  merge(meta_sub)

##..Melt for easy plotting..##
# plot_data <- subset_mixing_merged[,!names(subset_mixing_merged) 
#                                   %in% c('CD4_glyco_mixing','CD8_glyco_mixing')]
plot_data <- subset_mixing_merged
data.m <- melt(plot_data, id.vars=c('sample','gran_CFU','log_CFU','burden'))

##..Compare overall mixing..##

my_comparisons <- list(c('CD4_glyco_mixing','CD8_glyco_mixing'),
                 c('CD4_IDO_mixing','CD8_IDO_mixing'),
                 c('CD4_IDO_mixing','CD4_glyco_mixing'),
                 c('CD8_IDO_mixing','CD8_glyco_mixing'))

ggplot(data.m, aes(variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, color = NA, scale = "width") +
  stat_compare_means(comparisons = my_comparisons, method= "wilcox.test") +
  theme_minimal() +
  stat_summary(fun = median, geom = "point",size=3) +
  theme(legend.position = 'none')

compare_means(value ~ variable, data.m, method= "t.test")

##..Look at relationship between mixing and metabolic region area..##

# append area data 
metabolic_area_per_sample <- metabolic_area_data[row.names(unique(metabolic_area_data
                                                                  [,c("sample", "glyco_area","IDO1_area")])), 
                                                 names(metabolic_area_data) %in% c("sample", "glyco_area","IDO1_area")]
mixing_metabolic <- left_join(mixing_merged, metabolic_area_per_sample, by = c('sample'))
mixing_metabolic <- left_join(mixing_metabolic, area_data, by = c('sample'))

# Generate ratios and frequency of total
mixing_metabolic$glyco_ratio <- mixing_metabolic$glyco_area / mixing_metabolic$IDO1_area
mixing_metabolic$glyco_freq <- mixing_metabolic$glyco_area / (mixing_metabolic$IDO1_area + mixing_metabolic$glyco_area)
mixing_metabolic$IDO_freq <- mixing_metabolic$IDO1_area / (mixing_metabolic$IDO1_area + mixing_metabolic$glyco_area)

mixing_metabolic$glyco_freq_total <- mixing_metabolic$glyco_area / (mixing_metabolic$cellular_px + mixing_metabolic$necrosis_px)
mixing_metabolic$IDO_freq_total <- mixing_metabolic$IDO1_area / (mixing_metabolic$cellular_px + mixing_metabolic$necrosis_px)
mixing_metabolic$total_ratio <- mixing_metabolic$glyco_freq_total / mixing_metabolic$IDO_freq_total

# melt
id_cols <- c('sample', 'log_CFU','burden','total_mixing')
measure_cols <- c('mcore_mixing','glyco_mixing','ido_mixing',
                  'glyco_area','IDO1_area','glyco_ratio','glyco_freq','IDO_freq',
                  'glyco_freq_total','IDO_freq_total','total_ratio')
data_melted <- melt(mixing_metabolic, id.vars = id_cols, measure.vars = measure_cols)

# evaluate all relationships with mixing
ggplot(data_melted, aes(total_mixing, value)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_minimal() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  facet_wrap(~variable, scale = "free")


# evaluate individual relatonship
ggplot(mixing_metabolic, aes(total_mixing, glyco_freq)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_minimal() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) 

ggscatterstats(
  data = mixing_metabolic,
  x    = glyco_freq,
  y    = total_mixing,
  type = "pearson")

ggscatterstats(
  data = data_freqs[data_freqs$category == 'pheno_of_totalimmune' &
                      data_freqs$variable == 'CD11c+_Mac',],
  x    = log_CFU,
  y    = freq_of_total,
  type = "spearman")


ggplot(metabolic_summary, aes(bin_num, log_CFU)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_minimal() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) 


ggscatterstats(
  data = metabolic_summary,
  x    = log_CFU,
  y    = bin_num,
  type = "pearson")

