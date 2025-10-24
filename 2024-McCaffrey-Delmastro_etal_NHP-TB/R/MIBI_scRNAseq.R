# NHP_TB_scRNAseq.R
# Created by: Erin McCaffrey
# Date created: 04/24/24
#
# Overview: This script conducts analysis of the scRNAseq data from Peters et al
# using outputs of the pathway enrichment analyses and differential gene
# expresion analysis provided by the original study authors. 

library(EnhancedVolcano)
library(ggplot2)
library(forcats)
library(viridis)
library(tidyverse)
library(reshape2)
library(cowplot)
library(gplots)

##..Step 1: Read in the data..##

deg_data <- read.csv('RM2_vs_RM3.csv')

IDO1_sigdb_data <- read.csv('up_in_RM3_MSigDB.csv')
hypoxia_sigdb_data <- read.csv('RM2_vs_RM3_MSigDB.csv')

IDO1_kegg_data <- read.csv('up_in_RM3_KEGG.csv')
hypoxia_kegg_data <- read.csv('RM2_vs_RM3_KEGG.csv')

metabolic_scores <- read.csv('v2_TableS6_MetabolicPathwayActivities.csv')

##..Step 2: Generate a volcano from the deg data..##

keyvals <- ifelse(deg_data$avg_log2FC < -1 & deg_data$p_val_adj < 0.05, 'blue',
                  ifelse(deg_data$avg_log2FC > 1 & deg_data$p_val_adj < 0.05, 'green','black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'green'] <- 'Glyco'
names(keyvals)[keyvals == 'black'] <- 'Neither'
names(keyvals)[keyvals == 'blue'] <- 'IDO1'

EnhancedVolcano(deg_data,
                lab = deg_data$gene,
                title = 'Glycolytic (RM2) versus IDO1 (RM3) Macrophages',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                legendLabels=c('NS','FC > 1.0','Adj. p < 0.05',
                               'Adj. p < 0.05 & FC > 1.0'),
                # col = c('black', 'black', 'grey', 'blue'),
                colCustom = keyvals,
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                arrowheads = FALSE,
                hline = c(0.05),
                ylab = '-log10(adj p)',
                xlab = 'log2FC',
                xlim = c(-5,5))

##..Step 3: Generate a bubble plot for the pathway data..##

IDO1_kegg_data$log_p <- -log10(IDO1_kegg_data$Adjusted.P.value)

ggplot(IDO1_kegg_data, aes(x = Odds.Ratio, y = fct_reorder(Term, Odds.Ratio,
                                                      .fun=median,.desc=F))) +
  geom_point(aes(color = Combined.Score, size = log_p), alpha = 1) +
  scale_color_viridis() +
  scale_size(range = c(1, 10)) + # Adjust the range of points size
  theme_set(theme_bw() +theme(legend.position = "bottom"))


IDO1_sigdb_data$log_p <- -log10(IDO1_sigdb_data$Adjusted.P.value)
hypoxia_sigdb_data$log_p <- -log10(hypoxia_sigdb_data$Adjusted.P.value)

IDO1_sigdb_data <- IDO1_sigdb_data %>%
  arrange(desc(Odds.Ratio))
IDO1_sigdb_data$Term <- fct_inorder(IDO1_sigdb_data$Term) %>% fct_rev()

hypoxia_sigdb_data <- hypoxia_sigdb_data %>%
  arrange(desc(Odds.Ratio))
hypoxia_sigdb_data$Term <- fct_inorder(hypoxia_sigdb_data$Term) %>% fct_rev()

max_score <- max(hypoxia_sigdb_data$Combined.Score)
max_p <- max(hypoxia_sigdb_data$log_p)

hypoxia_plot <- ggplot(hypoxia_sigdb_data, aes(x = Odds.Ratio, y = Term)) +
  geom_point(aes(color = Combined.Score, size = log_p), alpha = 1) +
  scale_color_viridis() +
  scale_size(range = c(1, 10)) + # Adjust the range of points size
  theme_set(theme_bw() + theme(legend.position = "bottom")) 
hypoxia_plot

IDO1_plot <- ggplot(IDO1_sigdb_data, aes(x = Odds.Ratio, y = Term)) +
  geom_point(aes(color = Combined.Score, size = log_p), alpha = 1) +
  scale_color_viridis() +
  scale_size(range = c(1, 10)) + # Adjust the range of points size
  theme_set(theme_bw() + theme(legend.position = "bottom")) 
IDO1_plot

cowplot::plot_grid(hypoxia_plot, IDO1_plot)

##..Step 4: Plot the metabolic score data..##

# transform the data to have the scores by population
metabolic_scores_cast <- dcast(metabolic_scores, int_level_2 ~ source, value.var = c('score'))
names(metabolic_scores_cast) <- 
  str_replace_all(names(metabolic_scores_cast), c(" " = ".", 
                                                  "/" = "",
                                                  "," = "."))

# plot the score for glycolysis in RM2 v RM3
plot_data <- metabolic_scores_cast[metabolic_scores_cast$int_level_2 %in%
                                     c('Recruited macrophages 2', 
                                       'Recruited macrophages 3'),]

ggplot(plot_data, aes(int_level_2, y=Glycolysis..Gluconeogenesis)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(legend.position = 'none')

# plot heatmap

heatmap_data <- as.matrix(metabolic_scores_cast[, 2:80])
rownames(heatmap_data) <- metabolic_scores_cast$int_level_2

custom_pal_div<-c("#0E1949","#1E356C","#31508C","#4272AE","#6A93C6","#98B1DA",
                           "#C8D0EF","#F8F0FE","#F0C5D8","#E19EB0","#D17486",
                           "#BD4B5C","#923346","#691D32","#43071E")
colfunc<-colorRampPalette(custom_pal_div)

heatmap.2(heatmap_data,
          scale = "col",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "col",
          trace = "none",
          col = colfunc(100),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),      
          sepcolor="black",
          colsep=0:ncol(t(hm_niche_markers)),
          rowsep=0:nrow(t(hm_niche_markers)),
          breaks=seq(-2.5,2.5, length.out=101))

