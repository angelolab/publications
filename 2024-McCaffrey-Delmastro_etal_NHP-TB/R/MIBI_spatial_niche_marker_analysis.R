# MIBI_spatial_niche_marker_analysis.R
# Date created: 03/24/24
# This script takes the niche output from Jolene for functional marker expression.
# It generates the plot Jolene provided broken down by niche and cell subset.
# it also generates one that is a mean of all cell types in a niche. 

library(reshape2)
library(tidyverse)
library(dplyr)
library(gplots)
library(pals)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(data.table)
library(EnhancedVolcano)

##..Read in the data..##

setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2/")
sc_data<-read.csv('cohort_cell_table.csv')
niche_data<-read.csv('./spatial_analysis/QUICHE/v2/tb_annotated_table_tb_binary_updated_all-sig.csv')

##..Combine the sc data and niche data..##
sc_data <- merge(niche_data, sc_data, by = c('tiled_label','sample','centroid.x','centroid.y'))

##..Get the list of niches  ..##
data_summary <- niche_data %>%
  group_by(new_label) %>%
  summarize(median = median(logFC))

high_enriched <- data_summary %>%
  filter(rank(desc(median))<=15)

low_enriched <- data_summary %>%
  filter(rank(desc(median))>83)

diff_niches <- rbind(high_enriched, low_enriched)
diff_niches_ordered <- arrange(diff_niches, median)

##..Adjust single cell data for pSTAT3..##

pSTAT3_pos <- c('sample53','sample59','sample4','sample17','sample55','sample16',
                'sample54','sample45','sample32','sample7','sample18','sample24',
                'sample35','sample42','sample41','sample50','sample51')
coll_pos <- c('sample16','sample17','sample18','sample32','sample35',
                  'sample41','sample42','sample50','sample51','sample55')

sc_data[!sc_data$sample %in% pSTAT3_pos, ]$pSTAT3 <- NA
sc_data[!sc_data$sample %in% coll_pos, ]$Collagen1 <- NA

# ##..Split apart the cell type column..##
# 
# sc_data$pheno_corrected <- sub('.*: ', '', sc_data$pheno_corrected)

##..Filter the expression data to only include the niches..##

# first need to separate the niche info
# hmap_data_splt <- hmap_data %>% separate(pheno_corrected, into=c("niche", "null", "cell_type"), sep= " ")
# hmap_data$niche <- hmap_data_splt$niche
# 
# hmap_summary <- hmap_data %>%
#   group_by(niche) %>%
#   summarise(across(everything(), list(mean)))

##..Subset data by phenotype and/or niche..##

# pheno <- 'Neutrophil'
# niches <- c("CD14+_Mac_Mono__CD206+_Mac__Neutrophil",
#             "CD11c+_Mac__CD8+Tcell__Neutrophil",
#             "CD68+_Mac__FN1+_Mac__Neutrophil",
#             "Bcell__CD14+CD11c+_Mac__Neutrophil",
#             "HLADR+_APC__Neutrophil__VIM+Stroma")
# hmap_data <- sc_data[sc_data$pheno_corrected == pheno &
#                        sc_data$new_label %in% niches,]
# hmap_data <- sc_data[sc_data$pheno_corrected == pheno,]

##..Generate heatmap..##

hmap_data <- sc_data
# markers <- c('MMP9','CHIT1','GLUT1','pIRF3','FTL','IDO1','pS6','iNOS','pSTAT3')
markers <- c('Arginase1', 'CD36', 'CHIT1', 'DC.LAMP', 'FTL', 'GLUT1', 'H3K9Ac', 
                     'H3K27me3', 'HLA.DR', 'HO.1', 'IDO1', 'IFNg', 'IL33', 'iNOS', 
                     'MMP9', 'pIRF3', 'pS6', 'pSMAD3', 'pSTAT1', 'pSTAT3','Fe','Collagen1')
# markers <- names(sc_data[,2:22])
groups <- unique(hmap_data$new_label)
# groups_ordered <- groups[order(match(groups,diff_niches_ordered$new_label))]
groups_ordered <- diff_niches_ordered$new_label

hm_niche_markers <- matrix(, nrow = length(markers), ncol = length(groups_ordered))
for(i in 1:length(groups_ordered)) {
  temp_mat <- hmap_data[hmap_data[,"new_label"] == groups_ordered[i], markers]
  hm_niche_markers[,i] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_niche_markers) <- markers
colnames(hm_niche_markers) <- groups_ordered

hm_niche_markers[is.na(hm_niche_markers)] <- 0

##..Plot..##

custom_pal_div<-c("#0E1949","#1E356C","#31508C","#4272AE","#6A93C6","#98B1DA",
                           "#C8D0EF","#F8F0FE","#F0C5D8","#E19EB0","#D17486",
                           "#BD4B5C","#923346","#691D32","#43071E")
colfunc<-colorRampPalette(custom_pal_div)

heatmap.2(hm_niche_markers,
          scale = "row",
          Colv = F, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "row",
          trace = "none",
          col = colfunc(100),
          # col = rev(as.vector((brewer.spectral(100)))),
          # col = (as.vector((ocean.delta(100)))),
          # col = magma(256),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),      
          sepcolor="black",
          colsep=0:ncol(hm_niche_markers),
          rowsep=0:nrow(hm_niche_markers),
          breaks=seq(-2.5,2.5, length.out=101))
          # breaks=seq(-0.5,2.5, length.out=257))

##..Evaluate effect size of cell-marker combos..##

# function from David Glass
getEffectSize <- function(x, y, hedge=T) {
  # Calculates Hedge's g or Cohen's d effect size
  # Inputs:
  #   x - a numeric vector
  #   y - a numeric vector
  #   hedge - logical, if tru calculate Hedge's g, if false, Cohen's d
  # Ouptut:
  #   the effect size
  x.mean <- mean(x, na.rm = T)
  y.mean <- mean(y, na.rm = T)
  x.sd <- sd(x, na.rm = T)
  y.sd <- sd(y, na.rm = T)
  x.n <- length(x) - 1
  y.n <- length(y) - 1
  if (hedge) pooled.sd <- sqrt((x.n*x.sd^2 + y.n*y.sd^2)/(x.n+y.n))
  else pooled.sd <- sqrt((x.sd^2+y.sd^2)/2)
  return((x.mean-y.mean)/pooled.sd)
}

# reshape the data 
marker_data.m <- reshape2::melt(sc_data, 
                      id.vars = c('new_label','pheno_corrected'),
                      measure.vars = markers)

# drop giant cells because they are too few
marker_data.m <- droplevels(marker_data.m[!marker_data.m$pheno_corrected == 'giant_cell',])

# create marker-cell column
marker_data.m$feature <- paste(marker_data.m$pheno_corrected,"-",marker_data.m$variable)

# remove the symbols from feature column
marker_data.m$feature <- gsub("+", "", marker_data.m$feature, fixed = TRUE)
marker_data.m$feature <- gsub(".", "", marker_data.m$feature, fixed = TRUE)
marker_data.m$feature <- gsub(" ", "", marker_data.m$feature, fixed = TRUE)

# append the niche grouping info
data_summary <- niche_data %>%
  group_by(new_label) %>%
  summarize(median = median(logFC))
high_burden_niches <- data_summary[data_summary$median > 0, ]$new_label
marker_data.m$niche_group <- "low_burden"
marker_data.m[marker_data.m$new_label %in% high_burden_niches,]$niche_group <- "high_burden"

# reshape to have 
# marker_data.cast <- reshape2::dcast(marker_data.m,
#                                     new_label ~ pheno_corrected + variable,
#                                     value.var = "value",
#                                     fun.aggregate = identity)
# marker_data.cast$niche_group <- 'low_burden'
# marker_data.cast[marker_data.cast$new_label %in% high_burden_niches, ]$niche_group <- 'high_burden'

# generate effect size data

# set up input as data table
dt <- as.data.table(marker_data.m)
columns <- unique(marker_data.m$feature)

# Create a dt to store the output.
stats.dt <- data.table(feature=columns,
                       low.high.effect=0,
                       low.high.p.value=0,
                       low.high.adj.p.value=0)

# Generate the effect sizes and p-values for 
for(col in columns) {
  high.v <- unname(unlist(dt[niche_group=="high_burden" & feature == col, "value"]))
  low.v <- unname(unlist(dt[niche_group=="low_burden" & feature == col, "value"]))
  stats.dt[feature==col, `:=`(low.high.effect=getEffectSize(high.v, low.v),
                              low.high.p.value=wilcox.test(high.v, low.v)$p.value)]
}

# Adjust the p-values using FDR method of multiple hypothesis correction
stats.dt$low.high.adj.p.value <- p.adjust(stats.dt$low.high.p.value, method="fdr")

# plot
EnhancedVolcano(stats.dt,
                lab = stats.dt$feature,
                title = 'High versus low burden niches',
                x = 'low.high.effect',
                y = 'low.high.adj.p.value',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 3.0,
                legendLabels=c('NS','Effect Size > 0.5','Adj. p < 0.05',
                               'Adj. p < 0.05 & Effect Size > 0.5'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                parseLabels = TRUE)

# custom plot (source: https://biostatsquid.com/volcano-plots-r-tutorial/)

stats.dt$diffexpressed <- "NO"
stats.dt$diffexpressed[stats.dt$low.high.effect >= 0.5 & 
                         stats.dt$low.high.adj.p.value < 0.05] <- "UP"
stats.dt$diffexpressed[stats.dt$low.high.effect <= -0.5 & 
                         stats.dt$low.high.adj.p.value < 0.05] <- "DOWN"

# set labeling 
down_features <- stats.dt[stats.dt$diffexpressed == "DOWN",]
down_features <- down_features[order(down_features$low.high.effect),]
down_feature_names <- down_features[1:10,]$feature
up_features <- stats.dt[stats.dt$diffexpressed == "UP",]
up_features <- up_features[order(up_features$low.high.effect),]
up_feature_names <- up_features[21:35,]$feature


stats.dt$delabel <- stats.dt$feature
stats.dt[!stats.dt$diffexpressed %in% c("UP","DOWN"),]$delabel <- NA
  


ggplot(data = stats.dt, aes(x = low.high.effect, y = -log10(low.high.adj.p.value), 
                            color = diffexpressed,
                            label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) + 
  scale_color_manual(values = c("#1B75BB", "grey", "#BE1E2D"), 
                     labels = c("Up in Low", "Not significant", "Up in High")) + 
  # coord_cartesian(xlim = c(-2.5, 2.5)) +  
  theme_classic() + 
  ggtitle('Cell-marker Pairs in High versus Low Burden Niches') +
  geom_text_repel(max.overlaps = Inf)

ggplot(data = stats.dt, aes(x = low.high.effect, y = -log10(low.high.adj.p.value), 
                            col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#1B75BB", "grey", "#BE1E2D"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  ggtitle('Cell-marker Pairs in High versus Low Burden Niches') + 
  geom_text_repel(max.overlaps = Inf)


# make a heatmap of the significant features 
diff_marker_pairs <- stats.dt[stats.dt$diffexpressed %in% c('UP','DOWN'),]$label





##..Evaluate log fold-change cell-marker pairs..##

# generate cell-type - marker combo summaries for burden groups separately
markers_low_burden <- marker_data.m[marker_data.m$niche_group == 'low_burden',] %>%
  group_by(pheno_corrected, variable) %>%
  summarize(mean_low = mean(value, na.rm = TRUE))

markers_high_burden <- marker_data.m[marker_data.m$niche_group == 'high_burden',] %>%
  group_by(pheno_corrected, variable) %>%
  summarize(mean_high = mean(value, na.rm = TRUE))

markers_all <- left_join(markers_low_burden, markers_high_burden, 
                         by = c('pheno_corrected', 'variable'))

markers_all$FC <- log2(markers_all$mean_high / markers_all$mean_low)

# convert to heatmap
markers.cast <- reshape2::dcast(markers_all,
                                    pheno_corrected ~ variable,
                                    value.var = "FC")
markers.hmap <- as.matrix(markers.cast[,2:23])
rownames(markers.hmap) <- markers.cast$pheno_corrected

heatmap.2(markers.hmap,
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          col = colfunc(100),
          # col = rev(as.vector((brewer.spectral(100)))),
          # col = (as.vector((ocean.delta(100)))),
          # col = magma(256),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),      
          sepcolor="black",
          colsep=0:ncol(markers.hmap),
          rowsep=0:nrow(markers.hmap),
          breaks=seq(-3,3, length.out=101))

# plot relationsips as violins

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', linewidth = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())

ggplot(marker_data.m, aes(x = niche_group, y = value, fill = pheno_corrected)) +
  geom_violin(aes(fill = niche_group)) +
  scale_fill_manual(values=c('#BE1E2D','#1B75BB')) +
  facet_wrap(~variable, ncol = 6, scales = "free_y") +
  theme +
  stat_compare_means(method= "wilcox.test", label = "p.format") +
  theme(legend.position = 'none') + 
  labs(x="Niche Group") +
  labs(y="Expression") 
