# Correlate EVT glycoprotein and glycosyltransferase transcript
# Author: Ke Leow
# Date: 04/23/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
sapply(list.files(pattern="[.]R$", path="R_scripts/functions/", full.names=TRUE), source)

#--------------------------------
# Get gene lists for EVT glycoprotein and glycosyltransferase
#--------------------------------
#get EVT glycoproteins
evt_glycoproteins <- read_csv('data/glycoproteomics/evt_HB_25_052725.csv') %>% pull(gene)

#--------------------------------
# Load data - expression of all EVT ROIs
#--------------------------------
data <- read_csv('data/EVT_Nanostring/expression_tab.csv')

#get EVT ROIs
meta <- read_csv('data/EVT_Nanostring/metadata_tab.csv') 
#how many EVT ROIs
meta %>% group_by(SegmentDisplayName) %>% summarise(n=n())
# # Step 1: Filter metadata for EVT segments
# evt_samples <- meta %>%
#   filter(str_detect(SegmentDisplayName, "EVT")) %>%
#   pull(sample_names)
# 
# # Step 2: Subset the expression table by matching column names to EVT samples
# expr_evt <- data[, colnames(data) %in% c("TargetName",evt_samples)]

#--------------------------------
# Load data - expression of all EVT ROIs
#--------------------------------
# Define galectin gene symbols
galectin_genes <- c("LGALS1", "LGALS2", "LGALS3", "LGALS4", 
                    "LGALS7", "LGALS7B", "LGALS8", "LGALS9", "LGALS12", "LGALS13", "LGALS14", "LGALS16")
# write.csv(galectin_genes, "data/glycoproteomics/galectin_genes.csv", row.names = F)

#get galectin interactions
string <- read_tsv("STRING_interactions/EVT_HB_LGALS/string_interactions_short (2).tsv")

# Standardize column to have all interactions as (gene, galectin)
# Swap columns so galectin always in column 'galectin' and partner in 'gene'
string_long <- string %>%
  rename(node1 = `#node1`) %>% 
  mutate(
    galectin = ifelse(node1 %in% galectin_genes, node1, 
                      ifelse(node2 %in% galectin_genes, node2, NA)),
    gene = ifelse(node1 %in% galectin_genes, node2, 
                  ifelse(node2 %in% galectin_genes, node1, NA))
  ) %>%
  filter(!is.na(galectin) & !is.na(gene))  # keep only galectin-gene pairs

# Summarize for each gene
summary_table <- string_long %>%
  group_by(gene) %>%
  summarise(
    num_galectin_interactions = n_distinct(galectin),
    galectins = paste(unique(galectin), collapse = ", "),
    median_score = median(combined_score)
  ) %>%
  arrange(desc(num_galectin_interactions)) %>% 
  filter(gene %in% evt_glycoproteins)

#get galectin interactors - 13 genes
summary_table %>% pull(gene)

# #annotate by median score
# gal_anno <- summary_table %>% 
#   filter(gene %in% map_glyco_transcr$gene_name) %>% 
#   select(gene, median_score) %>% 
#   rename(gal.score = median_score) %>% 
#   column_to_rownames("gene")

#annotate true/false
gal_anno <- summary_table %>% 
  # filter(gene %in% evt_glycoproteins) %>% 
  # mutate(gal_ligand = TRUE) %>% 
  select(gene, num_galectin_interactions) %>% 
  rename(num_gal = num_galectin_interactions) %>% 
  column_to_rownames("gene") 
# %>% 
#   mutate(across(everything(), as.factor))

#--------------------------------
# correlate selected genes
#--------------------------------
#select data for corr
dat1_corr <- subset(data, TargetName %in% evt_glycoproteins) %>% 
  column_to_rownames("TargetName") 

#select data for corr
dat2_corr <- subset(data, TargetName %in% galectin_genes) %>% 
  column_to_rownames("TargetName") 

#run correlation
# corr_outcomes_exprs(t(dat2_corr), t(dat1_corr), method = "pearson")
corr_outcomes_exprs(t(dat2_corr), t(dat1_corr), method = "spearman")

# Define breaks ensuring 0 is centered
# Define symmetric breaks around zero
lim <- 1  # Find the max absolute value in data
breaks <- seq(-lim, lim, length.out = 101)  # Ensure zero is centered

# Create the matrix of significance labels
sig_labels <- ifelse(pvals < 0.001, "***",
                     ifelse(pvals < 0.01, "**",
                            ifelse(pvals < 0.05, "*", "")))

annotation_colors <- list(
  num_gal = colorRampPalette(brewer.pal(n = 7, name ="Purples"))(100)
)

pdf("R_plots/EVT_Nanostring/corr_heatmap_evt_HB_LGALS_112025.pdf", width = 7, height = 3)
pheatmap(coeff,
         annotation_col = gal_anno,
         annotation_colors = annotation_colors,annotation_names_col = FALSE,
         display_numbers = sig_labels,
         number_color = "black",
         na_col = "white",     # color for NA cells
         angle_col = "90",
         breaks = breaks,
         border_color = "grey5",
         # fontsize_row = 5,     # adjust for row labels
         # fontsize_col = 5,     # adjust for column labels
         # show_rownames = FALSE,
         # annotation_colors = annotation_colors,
         color = rev(colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100))
)
dev.off()

