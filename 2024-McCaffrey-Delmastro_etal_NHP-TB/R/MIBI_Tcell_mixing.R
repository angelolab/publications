# MIBI_Tcell_mixing.R
# Cretaed by: Erin McCaffrey
# Date created: 03/21/24
#
# Overview: This script takes the output of the mixing score analysis and evaluates the 
# relationship between mixing between cell types and across zones

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggstatsplot)
library(ggsignif)

##..Step 1: Read in the data..##

total_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_whole_gran.csv')
mcore_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_mcore-only.csv')
glyco_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_glycozone-only.csv')
ido_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_IDOzone-only.csv')
metabolic_area_data <- read.csv('cell_cohort_data_metabolic_zones.csv')
area_data <- read.csv('./masks/all_samples_mask_data.csv')
meta_data <- read.csv('./cohort_metadata/study_cohort_metadata.csv')

##..Step 2: Analyze the mixing scores not broken down by cell type..##
total_mixing_sub <- unique(total_mixing[,c(1,2)])
colnames(total_mixing_sub) <- c('sample','total_mixing')

mcore_mixing_sub <- unique(mcore_mixing[,c(1,2)])
colnames(mcore_mixing_sub) <- c('sample','mcore_mixing')

glyco_mixing_sub <- unique(glyco_mixing[,c(1,2)])
colnames(glyco_mixing_sub) <- c('sample','glyco_mixing')

ido_mixing_sub <- unique(ido_mixing[,c(1,2)])
colnames(ido_mixing_sub) <- c('sample','ido_mixing')

##..Step 3: Merge the area and metadata..##

meta_sub <- meta_data[,c(1,8,9,10)]
mixing_merged <- merge(total_mixing_sub, mcore_mixing_sub) %>% 
  merge(glyco_mixing_sub) %>%
  merge(ido_mixing_sub) %>%
  merge(meta_sub)

##..Step 4: Melt for easy plotting..##

plot_data <- mixing_merged[,!names(mixing_merged) %in% c('total_mixing','mcore_mixing')]
data.m <- melt(plot_data, id.vars=c('sample','gran_CFU','log_CFU','burden'))

##..Step 5: Plot..##

# Compare overall mixing
ggplot(data.m, aes(variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, color = NA, scale = "width") +
  theme_minimal() +
  scale_fill_manual(values=c('#00F422','#0022F4')) +
  stat_summary(fun = median, geom = "point",size=3) 

compare_means(value ~ burden, data.m, method= "t.test")

##..Step 6: Break down CD4 and CD8 T cells..##

CD4_glyco_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_glycozone-only_CD4T.csv')
CD8_glyco_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_glycozone-only_CD8T.csv')
CD4_IDO1_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_IDOzone-only_CD4T.csv')
CD8_IDO1_mixing <- read.csv('./spatial_analysis/mixing/mac_mixing_score_IDOzone-only_CD8T.csv')

##..Step 7: Extract the total scores..##

CD4_glyco_sub <- unique(CD4_glyco_mixing[,c(1,2)])
colnames(CD4_glyco_sub) <- c('sample','CD4_glyco_mixing')

CD8_glyco_sub <- unique(CD8_glyco_mixing[,c(1,2)])
colnames(CD8_glyco_sub) <- c('sample','CD8_glyco_mixing')

CD4_IDO1_sub <- unique(CD4_IDO1_mixing[,c(1,2)])
colnames(CD4_IDO1_sub) <- c('sample','CD4_IDO_mixing')

CD8_IDO1_sub <- unique(CD8_IDO1_mixing[,c(1,2)])
colnames(CD8_IDO1_sub) <- c('sample','CD8_IDO_mixing')

##..Step 7: Merge the area and metadata..##

subset_mixing_merged <- merge(CD4_glyco_sub, CD8_glyco_sub) %>% 
  merge(CD4_IDO1_sub) %>%
  merge(CD8_IDO1_sub) %>%
  merge(meta_sub)

##..Step 8: Melt for easy plotting..##

plot_data <- subset_mixing_merged
data.m <- melt(plot_data, id.vars=c('sample','gran_CFU','log_CFU','burden'))

##..Step 9: Plot..##

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
