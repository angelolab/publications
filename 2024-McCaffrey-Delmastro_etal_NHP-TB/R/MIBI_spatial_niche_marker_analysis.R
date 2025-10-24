# MIBI_spatial_niche_marker_analysis.R
# Created by: Erin McCaffrey
# Date created: 03/24/24
#
# Overview: This script takes the QUICHE output for functional marker expression.
# It generates the plot Jolene provided broken down by niche and cell subset.
# It also generates one that is a mean of all cell types in a niche. 

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

##..Step 1: Read in the data..##

sc_data<-read.csv('cohort_cell_table.csv')
niche_data<-read.csv('./spatial_analysis/QUICHE/v2/tb_annotated_table_tb_binary_updated_all-sig.csv')

##..Step 2: Combine the sc data and niche data..##
sc_data <- merge(niche_data, sc_data, by = c('tiled_label','sample','centroid.x','centroid.y'))

##..Step 3: Get the list of niches  ..##
data_summary <- niche_data %>%
  group_by(new_label) %>%
  summarize(median = median(logFC))

high_enriched <- data_summary %>%
  filter(rank(desc(median))<=15)

low_enriched <- data_summary %>%
  filter(rank(desc(median))>83)

diff_niches <- rbind(high_enriched, low_enriched)
diff_niches_ordered <- arrange(diff_niches, median)

##..Step 4: Adjust single cell data for pSTAT3 and Coll1..##

pSTAT3_pos <- c('sample53','sample59','sample4','sample17','sample55','sample16',
                'sample54','sample45','sample32','sample7','sample18','sample24',
                'sample35','sample42','sample41','sample50','sample51')
coll_pos <- c('sample16','sample17','sample18','sample32','sample35',
                  'sample41','sample42','sample50','sample51','sample55')

sc_data[!sc_data$sample %in% pSTAT3_pos, ]$pSTAT3 <- NA
sc_data[!sc_data$sample %in% coll_pos, ]$Collagen1 <- NA

##..Step 5: Generate heatmap..##

hmap_data <- sc_data
markers <- c('Arginase1', 'CD36', 'CHIT1', 'DC.LAMP', 'FTL', 'GLUT1', 'H3K9Ac', 
                     'H3K27me3', 'HLA.DR', 'HO.1', 'IDO1', 'IFNg', 'IL33', 'iNOS', 
                     'MMP9', 'pIRF3', 'pS6', 'pSMAD3', 'pSTAT1', 'pSTAT3','Fe','Collagen1')
groups <- unique(hmap_data$new_label)
groups_ordered <- diff_niches_ordered$new_label

hm_niche_markers <- matrix(, nrow = length(markers), ncol = length(groups_ordered))
for(i in 1:length(groups_ordered)) {
  temp_mat <- hmap_data[hmap_data[,"new_label"] == groups_ordered[i], markers]
  hm_niche_markers[,i] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_niche_markers) <- markers
colnames(hm_niche_markers) <- groups_ordered

hm_niche_markers[is.na(hm_niche_markers)] <- 0

##..Step 6: Plot..##

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
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),      
          sepcolor="black",
          colsep=0:ncol(hm_niche_markers),
          rowsep=0:nrow(hm_niche_markers),
          breaks=seq(-2.5,2.5, length.out=101))

##..Step 7: Evaluate effect size of cell-marker combos..##

# function from wise oracle David Glass, PhD
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
