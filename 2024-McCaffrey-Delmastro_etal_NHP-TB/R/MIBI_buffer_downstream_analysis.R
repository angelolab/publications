# MIBI_buffer_downstream_analysis.R
# Date created: 06/12/24
# Created by: Erin McCaffrey
#
# Overview: This script conducts downstream analysis for the buffer analysis
# generated in ArcGIS and subsequently processed with MIBI_process_buffer_masks.py
# and MIBI_extract_buffer_features.py. The script does the following:
#
# Step 1: Preprocess - In MIBI_extract_buffer_features.py the buffer masks were
# processed with the cellular data and image data to extract the count and 
# sum total signal from each buffer. In this preprocessing step that data is
# imported and buffers without any cells are removed. The CFU data and multifocal 
# status is also appended. After this, two alternative representations of the 
# data are produced. In one the cellular and protein expression data are area-normalized
# to generate cellular densities and mean expression, respectively. 
#
# Step 2: Individual granuloma analysis - Here for a selected granuloma the data 
# is visualized: 1) Cellular data all together 2) One plot per cell type 3) Marker
# data all together 4) One plot per marker 5) Pseudotime style heatmap
#
# Step 3: Grouped granuloma analysis - In order to bin the analysis, the buffer 
# numer gets normalized to the max buffer number. Both the normalized and 
# unnormalized data gets plotted together and fitted. 
#
# R version 4.2.2
# 
# Packages needed: 
# readxl, dlyr, ggplot2, reshape, ggstatsplot, ggpubr

library(Hmisc)
library(tidyr) 
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggalt)
library(forcats)
library(ggridges)
library(gplots)
library(pals)


library(ggstatsplot)
library(ggpubr)
library(cowplot)
library(introdataviz)

###...Step 1: Preprocess...###

# read in data
setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
buffer_data <- read.csv('buffer_data.csv')
meta_data <- read.csv('./cohort_metadata/study_cohort_metadata.csv')
necrosis_data <- read.csv('./spatial_analysis/radial/necrosis_distance_data.csv')
metabolic_buffer_data <- read.csv('metabolic_buffer_data.csv')

# drop the unassigned samples
buffer_data <- droplevels(buffer_data[,!names(buffer_data) == 'Unassigned'])

# append the metadata for log_CFU and CFU group
meta_data_sub <- meta_data[,names(meta_data) %in% c('sample','log_CFU','burden')]
buffer_data_anno <- left_join(buffer_data, meta_data_sub, by = c('sample'))

# indicate single necrotic foci versus mixed/multifocal
single_foci_samples <- unique(necrosis_data$sample)
buffer_data_anno$histo_group <- 'multifocal'
buffer_data_anno[buffer_data_anno$sample %in% 
                   single_foci_samples, ]$histo_group <- 'unifocal'

# remove the rows with acellular buffers
buffer_data_filtered <- droplevels(buffer_data_anno[!buffer_data_anno$buffer_cell_num == 0,])

# convert remaining NaNs to zero
buffer_data_filtered[is.na(buffer_data_filtered)] <- 0

# turn collagen and pSTAT3 values to NA except for valid samples
pSTAT3_samples <- c('sample4','sample7','sample16','sample17','sample18',
                   'sample24','sample32','sample35','sample41','sample42',
                   'sample45','sample50','sample51','sample53','sample54',
                   'sample55','sample59')
coll_samples <- c('sample16','sample17','sample18','sample32','sample35',
                  'sample41','sample42','sample50','sample51','sample55')

buffer_data_filtered[!buffer_data_filtered$sample %in% pSTAT3_samples,]$pSTAT3 <- NA
buffer_data_filtered[!buffer_data_filtered$sample %in% coll_samples,]$Collagen1 <- NA

# generate area-normalized and cell frequency normalized versions of the data
# area norm
buffer_data_area_norm <- buffer_data_filtered
buffer_data_area_norm[,7:63] <- buffer_data_filtered[,7:63] / buffer_data_filtered$buffer_area
# frequency (and area-norm for proteins)
buffer_data_freq <- buffer_data_filtered
buffer_data_freq[,7:25] <- buffer_data_filtered[,7:25] / buffer_data_filtered$buffer_cell_num
buffer_data_freq[,26:63] <- buffer_data_filtered[,26:63] / buffer_data_filtered$buffer_area

# create a normalized buffer readouts
buffer_data_filtered$buffer_num_norm <- buffer_data_filtered$buffer_num / buffer_data_filtered$total_buffer_num
buffer_data_area_norm$buffer_num_norm <- buffer_data_area_norm$buffer_num / buffer_data_area_norm$total_buffer_num
buffer_data_freq$buffer_num_norm <- buffer_data_freq$buffer_num / buffer_data_freq$total_buffer_num

# bin 
bin_n = 20
buffer_data_filtered$buffer_num_norm_binned <- cut(buffer_data_filtered$buffer_num_norm, breaks = bin_n)
buffer_data_area_norm$buffer_num_norm_binned <- cut(buffer_data_area_norm$buffer_num_norm, breaks = bin_n)
buffer_data_freq$buffer_num_norm_binned <- cut(buffer_data_freq$buffer_num_norm, breaks = bin_n)

# replace the bin variables with numerical codes
bin_num <- 1:20
bin <- levels(buffer_data_area_norm$buffer_num_norm_binned)
bin_key <- data.frame(bin_num, bin)
names(bin_key) <- c('bin_num','buffer_num_norm_binned')

buffer_data_area_norm <- left_join(buffer_data_area_norm, bin_key, by = c('buffer_num_norm_binned'))

# append the metabolic data

# generate proportion of area in zone
metabolic_buffer_data$glyco_percent <- metabolic_buffer_data$buffer_glyco_area / metabolic_buffer_data$buffer_area
metabolic_buffer_data$IDO_percent <- metabolic_buffer_data$buffer_IDO1_area / metabolic_buffer_data$buffer_area
metabolic_buffer_data$enrichment <- log2(metabolic_buffer_data$glyco_percent / metabolic_buffer_data$IDO_percent)

# merge with buffer data 
buffer_data_metabolic_combined <-  left_join(buffer_data_area_norm, 
                                             metabolic_buffer_data, by = c('sample',
                                                                           'total_buffer_num',
                                                                           'buffer_num'))

###...Step 2: Individual granuloma analysis...###

# choose a granuloma
gran <- 'sample33'

# subset the data
gran_data_raw <- buffer_data_filtered[buffer_data_filtered$sample == gran,]
gran_data_area_norm <- buffer_data_metabolic_combined[buffer_data_metabolic_combined$sample == gran,]
gran_data_freq <- buffer_data_freq[buffer_data_freq$sample == gran, ]

# melt for plotting (select different visualizations)
cellular_cols <- names(gran_data_raw[,7:25])
protein_cols <- names(gran_data_raw[,26:63])
cellular_data.m <- melt(gran_data_area_norm, id.vars <- c('sample','buffer_num','bin_num'), 
                        measure.vars <- cellular_cols)
protein_data.m <- melt(gran_data_area_norm, id.vars <- c('sample','buffer_num','bin_num'), 
                       measure.vars <- protein_cols)

# get the color key for the cellular data
color_key <- read.csv("./keys/cell_color_key_buffers.csv")
plot_populations<-levels(factor(cellular_data.m$variable))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

# plot the cellular data in bulk and separately
single_gran_bulk_cell <- ggplot(cellular_data.m, aes(x=buffer_num, y=value, color = variable)) + 
  geom_point() + 
  geom_line() +
  theme_minimal() + 
  scale_color_manual(values = color)
single_gran_bulk_cell

single_gran_individual_cell <- ggplot(cellular_data.m, aes(x=buffer_num, y=value, color = variable)) + 
  geom_point() + 
  geom_line() +
  theme_minimal() + 
  scale_color_manual(values = color) + 
  facet_wrap(~variable, scales=c("free"), ncol = 10) 
single_gran_individual_cell

# create an uncounted version of the data for ggridges
raw_data.m <- melt(gran_data_raw, id.vars <- c('sample','buffer_num'), 
                   measure.vars <- cellular_cols)

uncounted_gran_data <- raw_data.m %>% 
  type.convert(as.is = TRUE) %>% 
  uncount(value)

# reorder the factors 

order <- c('FN1._Mac','FN1.Fibro','Neutrophil','CD11c._Mac','CD68._Mac',
           'CD14.CD11c._Mac','CD14._Mac_Mono','CD206._Mac','giant_cell',
           'CD163._Mac','CD4.Tcell','CD8.Tcell','MastCell','Bcell',
           'Immune_Other','SMA.Fibro','HLADR._APC','VIM.Stroma','Endothelial')

uncounted_gran_data$variable<-factor(uncounted_gran_data$variable, 
                                     levels = order)
uncounted_gran_data<-uncounted_gran_data[order(uncounted_gran_data$variable),]

# get the color key
plot_populations<-levels(factor(uncounted_gran_data$variable))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

# plot
# single_gran_cell_ridges <- ggplot(uncounted_gran_data, aes(x = buffer_num,
#                                                            y = variable,
#                                                            height = ..scaled..,
#                                                            fill = variable)) +
#   geom_density_ridges(stat = "density", alpha=0.8, adjust = 1.5) +
#   scale_fill_manual(values=color) +
#   theme_ridges()
# single_gran_cell_ridges

single_gran_cell_ridges <- ggplot(uncounted_gran_data, aes(x = buffer_num, 
                                                           y = variable, 
                                                           height = stat(density),
                                                           fill = variable)) + 
  geom_density_ridges(stat = "binline", 
                      bins = max(uncounted_gran_data$buffer_num), 
                      draw_baseline = FALSE,
                      scale = 3,
                      alpha = 0.8) +
  xlim(0,23) +
  scale_fill_manual(values=color) +
  theme_ridges() +
  coord_cartesian(clip = "off") +
  theme(legend.position = 'none')
single_gran_cell_ridges

# optionally subset the proteins for plotting
plot_markers <- c('Arginase1', 'CD36', 'CHIT1', 'DC.LAMP', 'FTL', 'GLUT1', 'H3K9Ac', 
                  'H3K27me3', 'HLA.DR', 'HO.1', 'IDO1', 'IFNg', 'IL33', 'iNOS', 
                  'MMP9', 'pIRF3', 'pS6', 'pSMAD3', 'pSTAT1', 'pSTAT3','Fe','Collagen1')
protein_data_sub <- protein_data.m[protein_data.m$variable %in% plot_markers,]


# plot the protein data in bulk and individually
single_gran_bulk_protein <- ggplot(protein_data_sub, aes(x=buffer_num, y=value, color = variable)) + 
  geom_point() + 
  geom_line() +
  theme_minimal() 
single_gran_bulk_protein

single_gran_individual_protein <- ggplot(protein_data_sub, aes(x=buffer_num, y=value)) + 
  geom_point(color = 'grey30') + 
  geom_line(color = 'grey30') +
  theme_minimal() + 
  theme(legend.position = 'None') +
  facet_wrap(~variable, scales=c("free"), ncol = 11) 
single_gran_individual_protein

# choose multiple factors to plot together
factors <- c('buffer_num','bin_num','CD4.Tcell','CD8.Tcell',
             'IDO1','GLUT1','glyco_percent','IDO_percent')
gran_data_factors <- gran_data_area_norm[,names(gran_data_area_norm) %in% factors]

# normalize them for plotting together
gran_data_factors_norm <- gran_data_factors
gran_data_factors_norm$CD4.Tcell<- gran_data_factors_norm$CD4.Tcell/max(gran_data_factors$CD4.Tcell)
gran_data_factors_norm$CD8.Tcell<- gran_data_factors_norm$CD8.Tcell/max(gran_data_factors$CD8.Tcell)
gran_data_factors_norm$GLUT1 <-gran_data_factors_norm$GLUT1/max(gran_data_factors$GLUT1)
gran_data_factors_norm$IDO1 <-gran_data_factors_norm$IDO1/max(gran_data_factors$IDO1)
gran_data_factors_norm$glyco_percent <-gran_data_factors_norm$glyco_percent/max(gran_data_factors$glyco_percent)
gran_data_factors_norm$IDO_percent <-gran_data_factors_norm$IDO_percent/max(gran_data_factors$IDO_percent)

# melt
factor_data.m <- melt(gran_data_factors_norm,
                      id.vars = c('bin_num','buffer_num'),
                      measure.vars = c('CD8.Tcell',
                                       'CD4.Tcell',
                                       'IDO_percent',
                                       'glyco_percent'))

# plot
single_gran_factors <- ggplot(factor_data.m, aes(x=buffer_num, y=value)) + 
  geom_bar(stat = 'identity',
           aes(buffer_num, value, group=variable, fill=variable),
           data = . %>% filter(variable %in% c("IDO_percent" , "glyco_percent"))) +
  scale_fill_manual(values = c('#0022F2','#00F222')) +
  geom_line(aes(buffer_num, value, group=variable, colour=variable),
            data = . %>% filter(variable %in% c("CD8.Tcell" , "CD4.Tcell"))) +
  scale_color_manual(values = c('#ABE3F3','#3BBCA8')) +
  theme_minimal() +
  ggtitle("Sample 40")
single_gran_factors

# create a pseudotime-style heatmap
hmap_cols <- c(plot_populations, plot_markers)
hmap_data <- as.matrix(gran_data_area_norm[,names(gran_data_area_norm) %in% hmap_cols])
rownames(hmap_data) <- gran_data_area_norm$buffer_num
hmap_data.t <- t(hmap_data)

# order for nice visualization
hmap_data.t_reordered <- hmap_data.t[order(apply(hmap_data.t,1,which.max)), ]

# plot
heatmap.2(hmap_data.t_reordered,
          scale = "row",
          Colv = F, Rowv = F,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          col = magma(256),
          # col = cubehelix(256),
          # col = rev(as.vector((brewer.rdbu(100)))),
          # col = rev(as.vector((brewer.spectral(100)))),
          density.info = 'none',
          key.title = '',
          breaks=seq(-1,2.5, length.out=257))


###...Step 3: Grouped granuloma analysis...###

# melt for plotting (select different visualizations)
cellular_data_all.m <- melt(buffer_data_area_norm, 
                            id.vars <- c('sample','buffer_num','buffer_num_norm', 'bin_num',
                                         'buffer_num_norm_binned','burden','log_CFU','histo_group'), 
                            measure.vars <- cellular_cols)
protein_data_all.m <- melt(buffer_data_area_norm, 
                           id.vars <- c('sample','buffer_num','buffer_num_norm', 'bin_num',
                                        'buffer_num_norm_binned','burden','log_CFU','histo_group'),
                           measure.vars <- protein_cols)

# get the color key for the cellular data
plot_populations<-levels(factor(cellular_data_all.m$variable))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

# optionally subset the data

# plot the cellular data in bulk and separately
all_gran_bulk_cell <- ggplot(cellular_data_all.m, aes(x=as.numeric(bin_num), y=value, color = variable)) + 
  # geom_point() + 
  geom_line() +
  theme_minimal() + 
  scale_color_manual(values = color)
all_gran_bulk_cell

all_gran_individual_cell <- ggplot(cellular_data_all.m, aes(x=as.numeric(bin_num), y=value, color = variable)) + 
  geom_point() +
  geom_smooth(method = "loess", color = 'black') +
  # geom_line() +
  theme_minimal() + 
  scale_color_manual(values = color) +
  theme(legend.position = 'None') +
  facet_wrap(~variable, scales=c("free"), ncol = 10) 
all_gran_individual_cell

# create an uncounted version of the data for ggridges
raw_data_all.m <- melt(buffer_data_filtered, 
                       id.vars <- c('sample','buffer_num','buffer_num_norm',
                                    'buffer_num_norm_binned','burden','log_CFU','histo_group'), 
                       measure.vars <- cellular_cols)

uncounted_data_all <- raw_data_all.m %>% 
  type.convert(as.is = TRUE) %>% 
  uncount(value)

# reorder the factors 
order <- c('FN1._Mac','FN1.Fibro','Neutrophil','CD11c._Mac','CD68._Mac',
           'CD14.CD11c._Mac','CD14._Mac_Mono','CD206._Mac','giant_cell',
           'CD163._Mac','CD4.Tcell','CD8.Tcell','MastCell','Bcell',
           'Immune_Other','SMA.Fibro','HLADR._APC','VIM.Stroma','Endothelial',
           'Unassigned')

uncounted_data_all$variable<-factor(uncounted_data_all$variable, 
                                     levels = order)
uncounted_data_all<-uncounted_data_all[order(uncounted_data_all$variable),]

# get the color key
plot_populations<-levels(factor(uncounted_data_all$variable))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

# # plot ridges
# all_gran_cell_ridges <- ggplot(uncounted_data_all, aes(x = buffer_num_norm, 
#                                                            y = variable, 
#                                                            height = stat(density),
#                                                            fill = burden)) + 
#   geom_density_ridges(stat = "binline", 
#                       bins = 30,
#                       draw_baseline = FALSE,
#                       scale = 3,
#                       alpha = 0.8) +
#   scale_fill_manual(values=color) +
#   theme_ridges() +
#   coord_cartesian(clip = "off") 
# all_gran_cell_ridges
# 
# all_gran_cell_ridges <- ggplot(uncounted_data_all, aes(x = buffer_num_norm,
#                                                            y = variable,
#                                                            height = ..scaled..,
#                                                            fill = variable)) +
#   geom_density_ridges(stat = "density", 
#                       alpha=0.8, 
#                       adjust = 1.5,
#                       aes(fill = burden)) +
#   theme_ridges() + 
#   coord_cartesian(clip = "off") 
# all_gran_cell_ridges
# 
# all_gran_cell_density <- ggplot(cellular_data_all.m, aes(buffer_num, fill = burden)) +
#   stat_density(aes(weight = value, y = variable , height = after_stat(density)), 
#                geom = 'density_ridges', position = 'identity', alpha = 0.8) + 
#   theme_minimal() + 
#   coord_cartesian(clip = "off") +
#   facet_wrap(~variable, scales=c("free")) 
# all_gran_cell_density

# optionally subset the proteins for plotting
protein_data_all_sub <- protein_data_all.m[protein_data_all.m$variable %in% plot_markers,]
# 
# # plot the protein data in bulk and individually
# all_gran_bulk_protein <- ggplot(protein_data_all_sub, aes(x=buffer_num, y=value, color = variable)) + 
#   geom_point() + 
#   geom_line() +
#   theme_minimal() 
# all_gran_bulk_protein
# 
all_gran_individual_protein <- ggplot(protein_data_all_sub, aes(x=bin_num, y=value)) +
  geom_point(color = 'grey') +
  geom_smooth(method = "loess", color = 'black') +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~variable, scales=c("free"), ncol = 11)
all_gran_individual_protein

# create a pseudotime-style heatmap
# data_subset <- buffer_data_area_norm[buffer_data_area_norm$burden == 'low', ]
data_subset <- buffer_data_area_norm
hmap_data_all <- as.matrix(data_subset[,names(data_subset) %in% hmap_cols])

# nrows = max(data_subset$buffer_num)
rows <- unique(data_subset$buffer_num_norm_binned)
rows <- rows[order(rows)]
nrows <- length(unique(data_subset$buffer_num_norm_binned))

hm_allclusters <- matrix(, nrow = nrows, ncol = length(hmap_cols))
for(i in 1:max(nrows)) {
  temp_mat <- data_subset[data_subset[,"buffer_num_norm_binned"] == rows[i], hmap_cols]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_allclusters) <- rows
colnames(hm_allclusters) <- hmap_cols
hm_all_clusters.t <- t(hm_allclusters)

# order for nice visualization
hm_all_clusters.t_reordered <- hm_all_clusters.t[order(apply(hm_all_clusters.t,1,which.max)), ]

# plot
heatmap.2(hm_all_clusters.t_reordered,
          scale = "row",
          Colv = F, Rowv = F,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "none",
          trace = "none",
          col = magma(256),
          # col = cubehelix(256),
          # col = rev(as.vector((brewer.rdbu(100)))),
          # col = rev(as.vector((brewer.spectral(100)))),
          density.info = 'none',
          key.title = '',
          breaks=seq(-1,2, length.out=257))

###...Step 4: Metabolic zonation analysis...###

# generate boxplot visualization of zone preference
plot_data <- buffer_data_metabolic_combined[,names(buffer_data_metabolic_combined) %in%
                                              c('sample','glyco_percent','IDO_percent','bin_num')]
plot_data.m <- melt(plot_data,id.vars = c('sample','bin_num'), measure.vars = c('IDO_percent','glyco_percent'))

ggplot(data = plot_data.m, aes(x = as.factor(bin_num), y = value, fill=variable)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c('#0022F4','#00F422')) +
  # theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total Macrophages') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

# generate boxplot visualization of enrichment
ggplot(data = buffer_data_metabolic_combined, aes(x = as.factor(bin_num), y = enrichment)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total Macrophages') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

# generate single heatmap rows for each zone
metabolic_summary <- buffer_data_metabolic_combined %>%
  group_by(bin_num) %>%
  summarise(mean = mean(IDO_percent, na.rm = TRUE))

data_max <- max(metabolic_summary$mean)

# plot
fill <- colorRampPalette(c("white", "#0022F4"))(n=100)
ggplot(metabolic_summary, aes(x = as.factor(bin_num), y = 1, fill = mean)) + 
  scale_fill_gradientn(limits=c(0, data_max), colours=fill) +
  theme_bw() +
  geom_tile()


metabolic_summary <- buffer_data_metabolic_combined %>%
  group_by(sample) %>%
  mutate(max_score = max(CD4.Tcell)) %>% 
  ungroup() %>% 
  filter(CD4.Tcell==max_score)

colnames(metabolic_summary) <- c('sample','CD8.Tcell')

metabolic_summary <- left_join(metabolic_summary, 
                               buffer_data_metabolic_combined,
                               by = c('sample','CD8.Tcell'))

