# MIBI_buffer_correlation_analysis.R
# Date created: 08/23/24
# Created by: Erin McCaffrey
#
# Overview: This script reads on the buffer features for each sample, which
# comprises both cellular densities, protein expression patterns, and the 
# portion of each buffer belonging to the normoxic versus hypoxic zones. It then
# generates buffer-by-buffer correlations for all features and visualizes those
# in a heatmap. 

# necessary libraries

library(Hmisc)
library(WGCNA)
library(corrplot)
library(pals)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(reshape2)
library(psych)
library(tibble)
library(dplyr)
library(viridis)

###...Step 1: Preprocess...###

# read in the data
buffer_feature_data <- read.csv('all_grans_buffer_features.csv')

# generate a metabolic enrichment
buffer_feature_data$metabolic_enrichment <- log2(buffer_feature_data$glyco_percent / 
                                                   buffer_feature_data$IDO_percent)


###...Step 2: Generate all correlations...###

# define correlation features
cell_features <- colnames(buffer_feature_data[,6:24])
marker_features <- c('Arginase1', 'CD36', 'CHIT1', 'DC.LAMP', 'FTL', 'GLUT1', 'H3K9Ac', 
                                     'H3K27me3', 'HLA.DR', 'HO.1', 'IDO1', 'IFNg', 'IL33', 'iNOS', 
                                     'MMP9', 'pIRF3', 'pS6', 'pSMAD3', 'pSTAT1', 'pSTAT3','Fe','Collagen1')
# marker_features <- colnames(buffer_feature_data[,25:62])
metabolic_features <- colnames(buffer_feature_data[,69:73])

# subset the features for correlation
corr_data <- buffer_feature_data[,c(cell_features, marker_features, metabolic_features)]

# generate big correlation matrix for cellular and marker features
corr_data_mat <- as.matrix(corr_data)
corr <- corAndPvalue(corr_data_mat, method='spearman')

# generate a correlation matrix for the metabolic enrichment with filtration to 
# only include finite values
corr_data_finite <- corr_data[is.finite(corr_data$metabolic_enrichment),]
corr_data_finite_mat <- as.matrix(corr_data_finite)
corr_enrichment <- corAndPvalue(corr_data_finite_mat, method='spearman')

# pull out and append the metabolic enrichment data to the p value and coefficient mats
corrmatrix <- as.data.frame(corr$cor)
corrmatrix_enrichment <- as.data.frame(corr_enrichment$cor)
corrmatrix$metabolic_enrichment <- corrmatrix_enrichment$metabolic_enrichment
corrmatrix[46,] <- corrmatrix_enrichment[46,]

# pull out and adjust p values
pvalues <- as.data.frame(corr$p)
pvalues_enrichment <- as.data.frame(corr_enrichment$p)
pvalues$metabolic_enrichment <- pvalues_enrichment$metabolic_enrichment
pvalues[46,] <- pvalues_enrichment[46,]

pvalues_adj <- pvalues
ltri <- lower.tri(pvalues_adj)
utri <- upper.tri(pvalues_adj)
pvalues_adj[ltri] <- p.adjust(pvalues_adj[ltri], method = "fdr")
pvalues_adj[utri] <- t(pvalues_adj)[utri]

# plot
# corrplot(as.matrix(corrmatrix), order = "hclust", tl.col = "black")
# cross out statistically insignificant, customize other features
custom_col <- colorRampPalette(ocean.curl(255))
corrplot(as.matrix(corrmatrix), 
         method = 'color',
         order = "hclust", 
         tl.col = "black", 
         # p.mat = corr$p, 
         # sig.level = 0.05, 
         # insig = "pch", 
         # pch = 4,
         # pch.cex = 0.5,
         col=custom_col(200))

###...Step 3: Plot individual relationships...###

# plot_data <- buffer_feature_data
plot_data <- buffer_feature_data[is.finite(buffer_feature_data$metabolic_enrichment),]

ggplot(plot_data, aes(x=metabolic_enrichment, y = CD11c._Mac)) + 
  geom_point(aes(color = bin_num)) +
  scale_color_viridis() + 
  scale_y_log10() +
  theme_bw()
  

individual_corr <-corr.test(plot_data$metabolic_enrichment, 
                            plot_data$CD11c._Mac, method = 'spearman')
individual_corr$r
individual_corr$p

###...Step 4: Plot the correlation coefficients for the metabolic score...###

# subset the correlations
metabolic_coefficients <- as.data.frame(corrmatrix_enrichment$metabolic_enrichment)
metabolic_coefficients$feature <- rownames(corrmatrix_enrichment)
names(metabolic_coefficients) <- c('spearman_rho','feature')

# get an adjust the p-values just for these comparisons
metabolic_p <- pvalues_enrichment$metabolic_enrichment
metabolic_p_adj <- p.adjust(metabolic_p, method = "fdr")
metabolic_coefficients$p_adj <- metabolic_p_adj

# sort in descending order
metabolic_coefficients <- metabolic_coefficients %>% 
  arrange(desc(spearman_rho))

# remove the redundant features
remove <- c('metabolic_enrichment','glyco_percent','buffer_glyco_area',
           'buffer_IDO1_area','IDO_percent')
metabolic_coefficients_filtered <- metabolic_coefficients[!metabolic_coefficients$feature %in% remove,]

# optionally remove the non-statistically significant features
metabolic_coeffs_sig <- metabolic_coefficients_filtered[metabolic_coefficients_filtered$p_adj < 0.05,]

# plot as a bar chart
ggplot(metabolic_coeffs_sig, aes(reorder(feature, -spearman_rho), spearman_rho)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

