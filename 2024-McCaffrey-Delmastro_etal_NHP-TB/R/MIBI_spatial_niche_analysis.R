# MIBI_spatial_niche_analysis.R
# Date created: 03/21/24
# This script takes the niche output from Jolene and visualizes the differentially
# abundant niches between controlling and non-controlling granulomas. It also
# visualizes the expression of markers in each niche in a heatmap format.

library(reshape2)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)
library(ggnewscale)
library(gplots)
library(RColorBrewer)
library(pals)
library(tibble)
library(dendextend)
library(UpSetR)
library(ggpubr)
library(ComplexHeatmap)

##..Read in the data..##
setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2/")
data<-read.csv('./spatial_analysis/QUICHE/v2/tb_annotated_table_tb_binary_updated_all-sig.csv')

##..Assign the niches as being a high or low burden-enriched niche..##
data_summary <- data %>%
  group_by(new_label) %>%
  summarize(median = median(logFC))

high_burden_niches <- data_summary[data_summary$median > 0, ]$new_label
low_burden_niches <- data_summary[data_summary$median < 0, ]$new_label

data$enrichment <- 'high'
data[data$new_label %in% low_burden_niches, ]$enrichment <- 'low'

##..Create binary significance level..##
data$sig <- 'True'
data[data$SpatialFDR >= 0.05, ]$sig <- 'False'

##..Optionally filter to only include the top 50 differential..##
high_enriched <- data_summary %>%
  filter(rank(desc(median))<=15)

low_enriched <- data_summary %>%
  filter(rank(desc(median))>83)

diff_niches <- c(high_enriched$new_label, low_enriched$new_label)
# plot_data <- data[data$new_label %in% diff_niches, ]
plot_data <- data

##..Plot as ascending violins..##

ggplot(plot_data, aes(reorder(new_label, logFC, FUN = median), y=logFC, fill=enrichment)) +
  geom_violin(trim = FALSE, color = NA, scale = "width") +
  scale_fill_manual(values = c('#BE1E2D','#1B75BB')) +
  new_scale_color() +
  geom_point(aes(color = sig), fill=NA, alpha = 0.7, stroke = NA,
             position = position_jitter(width = 0.05))  +
  scale_color_manual(values = c('#696969','#000000')) + 
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_text(angle=90,hjust=1)) + 
  labs(y="log2FC") 


##..Append single cell data..##

sc_data <- read.csv('cell_cohort_data_metabolic_zones.csv')
niche_data_sc <- merge(data, sc_data, by = c('tiled_label','sample'))

##..Generate heatmap of niches by cell frequency..##

niche_cell_freq <- niche_data_sc %>% 
  group_by(new_label, .drop = F ) %>%
  count(pheno_corrected) %>%
  mutate(freq_of_total = prop.table(n)) %>%
  rename('variable' = 'pheno_corrected') %>%
  mutate(category = 'pheno_of_total')

niche_cell_hmap <- dcast(niche_cell_freq, new_label ~ variable, value.var = 'freq_of_total')  
niche_cell_hmap[is.na(niche_cell_hmap)] <- 0
niche_cell_hmap.m <- as.matrix(niche_cell_hmap[,2:20])
rownames(niche_cell_hmap.m) <- niche_cell_hmap$new_label

hmap <- heatmap.2(niche_cell_hmap.m[diff_niches,],
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          col = rev(as.vector((brewer.spectral(100)))),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),      
          sepcolor="white",
          colsep=0:ncol(niche_cell_hmap.m),
          rowsep=0:nrow(niche_cell_hmap.m),
          breaks=seq(0,1, length.out=101))

##..Plot binary matrix of the niche cell members..##

# generate binary heatmap
cell_types <- unique(niche_data_sc$pheno_corrected)
cell_types <- cell_types[!cell_types == 'giant_cell']
niche_members <- matrix(nrow = length(cell_types), ncol = length(diff_niches))

colnames(niche_members) <- diff_niches
rownames(niche_members) <- cell_types

for(i in 1:length(diff_niches)) {
  diff_niche <- diff_niches[i]
  niche_cell_types <- strsplit(diff_niche, split = '__')[[1]]
  for(j in 1:length(niche_cell_types)) {
    niche_cell_type <- niche_cell_types[j]
    niche_members[niche_cell_type, diff_niche] <- 1 
  }
}

niche_members[is.na(niche_members)] <- 0

# order to match cell frequency hierarchical order and niche order

# niche order
diff_niches_data_ordered <- data[data$new_label %in% diff_niches,] %>% 
  group_by(new_label) %>% 
  mutate(m = median(logFC)) %>% 
  arrange(desc(m)) %>% select(-m)
diff_niche_order <- rev(unique(diff_niches_data_ordered$new_label))

# cell order
cell_order <- colnames(niche_cell_hmap.m[diff_niches,])[hmap$colInd]
cell_order <- cell_order[!cell_order == 'giant_cell']

niche_members_ordered <- niche_members[cell_order, diff_niche_order]

# plot
ggplot(melt(niche_members_ordered), aes(Var2, Var1, fill = value, alpha = value)) + 
  geom_tile(colour = "gray50") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# make an upset plot
combo_mat <- make_comb_mat(t(niche_members_ordered))
upset(as.data.frame(t(niche_members_ordered)), 
      nsets = 30)

##..Determine niche clustering to choose representative niches..##

dist_mat <- dist(niche_cell_hmap, method="euclidean")
hc <- hclust(dist_mat, method = 'complete')
plot(hc)
dend <- as.dendrogram(hc)
plot(dend)
dend_19 <- color_branches(dend, h = 0.55)
plot(dend_19)
groups <- as.data.frame(cutree(hc, k = 19))
groups$new_label <-rownames(groups)
colnames(groups) <- c('new_label','cluster')

scree_data <- hc$height %>% as_tibble() %>%
  add_column(groups = length(hc$height):1) %>%
  rename(height=value)

ggplot(scree_data, aes(x=groups, y=height)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 19) +
  theme_minimal()

##..Define metabolic score of each niche..##

niche_IDO_count <- niche_data_sc %>% 
  group_by(new_label) %>%
  summarize(n_IDO1_zone = sum(IDO1_zone))

niche_glyco_count <- niche_data_sc %>% 
  group_by(new_label) %>%
  summarize(n_glyco_zone = sum(glyco_zone))

niche_total_count <- niche_data_sc %>%
              group_by(new_label) %>%
              tally()
niche_metabolic_freq <- merge(niche_IDO_count, 
                              merge(niche_glyco_count, niche_total_count),
                              by = c('new_label'))              
niche_metabolic_freq$glyco_freq <- niche_metabolic_freq$n_glyco_zone / niche_metabolic_freq$n
niche_metabolic_freq$IDO_freq <- niche_metabolic_freq$n_IDO1_zone / niche_metabolic_freq$n
niche_metabolic_freq$metabolic_enrichment <- log2(niche_metabolic_freq$n_glyco_zone /
                                                    niche_metabolic_freq$n_IDO1_zone)

# reorder to match the beeswarm plot
data_ordered <- data %>% 
  group_by(new_label) %>% 
  mutate(m = median(logFC)) %>% 
  arrange(desc(m)) %>% select(-m)
niche_order <- rev(unique(data_ordered$new_label))

niche_metabolic_freq <- niche_metabolic_freq %>%
  mutate(across(new_label, factor, levels=niche_order))

# optionally subset
# metabolic_plot_data <- niche_metabolic_freq[niche_metabolic_freq$new_label %in% diff_niches, ]
metabolic_plot_data <- niche_metabolic_freq

# plot
fill <- colorRampPalette(c("white", "#0022F4"))(n=100)
ggplot(metabolic_plot_data, aes(x = new_label, y = 1, fill = IDO_freq)) + 
  scale_fill_gradientn(limits=c(0, max(niche_metabolic_freq$IDO_freq)), colours=fill) +
  theme_bw() +
  geom_tile() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fill <- colorRampPalette(c("white", "#00F422"))(n=100)
ggplot(metabolic_plot_data, aes(x = new_label, y = 1, fill = glyco_freq)) + 
  scale_fill_gradientn(limits=c(0, max(niche_metabolic_freq$glyco_freq)), colours=fill) +
  theme_bw() +
  geom_tile() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##..Compare metabolic enrichment in high versus low burden-enriched niches..##

# append the enrichment information
niche_metabolic_freq$group <- 'high_burden'
niche_metabolic_freq[niche_metabolic_freq$new_label %in% low_burden_niches, ]$group <- "low_burden"

# plot metabolic attributes between groups

ggplot(niche_metabolic_freq, aes(group, y=metabolic_enrichment)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_violin(trim = FALSE, scale = "width") +
  stat_compare_means(method= "wilcox.test") +
  theme_bw() +
  stat_summary(fun = median, geom = "point",size=3) 


# integrate mean log2FC data
niche_metabolic_freq_logFC <- merge(niche_metabolic_freq, data_summary, by = c('new_label'))

# evaluate upper and lower most niches
high <- top_n(niche_metabolic_freq_logFC, 10, median)
low <- top_n(niche_metabolic_freq_logFC, -10, median)
caps <- rbind(high, low)

ggplot(caps, aes(group, y=metabolic_enrichment)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_violin(trim = FALSE, scale = "width") +
  stat_compare_means(method= "wilcox.test") +
  theme_bw() +
  stat_summary(fun = median, geom = "point",size=3) 

  
##..Append phenotype data and export..##

niche_data_sc <- niche_data_sc[,1:12]
niche_data_sc$lineage_category <- niche_data_sc$pheno_corrected
macs <- c("CD11c+_Mac","CD14+CD11c+_Mac","CD14+_Mac_Mono","CD163+_Mac","CD206+_Mac","CD68+_Mac",
                  "FN1+_Mac",'giant_cell')
lymphs <- c("Bcell","CD4+Tcell","CD8+Tcell")
grans <- c("Neutrophil","MastCell")
imm_other <- c("Immune_Other")
non_immune <-c("VIM+Stroma","Endothelial","FN1+Fibro","SMA+Fibro","HLADR+_APC")

niche_data_sc[niche_data_sc$pheno_corrected %in% macs,]$lineage_category <- "macrophage"
niche_data_sc[niche_data_sc$pheno_corrected %in% lymphs,]$lineage_category <- "lymphocyte"
niche_data_sc[niche_data_sc$pheno_corrected %in% grans,]$lineage_category <- "granulocyte"
niche_data_sc[niche_data_sc$pheno_corrected %in% imm_other,]$lineage_category <- "immune-other"
niche_data_sc[niche_data_sc$pheno_corrected %in% non_immune,]$lineage_category <- "non-immune"

write.csv(niche_data_sc, "NHP_TB_QUICHE_subset-annotated.csv", row.names = F)
