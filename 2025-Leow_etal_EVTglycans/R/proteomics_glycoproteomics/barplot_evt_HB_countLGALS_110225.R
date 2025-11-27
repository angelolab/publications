# Map DE glycoproteins to NS transcriptomics data
# Author: Ke Leow
# Date: 03/15/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)

#--------------------------------
# Load data
#--------------------------------
select_genes <- read_csv('data/glycoproteomics/evt_HB_25_052725.csv') %>% pull(gene)

# Define galectin gene symbols
galectin_genes <- c("LGALS1", "LGALS2", "LGALS3", "LGALS4", 
                    "LGALS7", "LGALS7B", "LGALS8", "LGALS9", "LGALS12", "LGALS13", "LGALS14", "LGALS16")
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
  filter(gene %in% select_genes)

# Barplot of number of galectins each gene interacts with
pdf("R_plots/glycoproteomics/barplot_evtHB_galectin_interactions.pdf", , width = 4, height = 4)
ggplot(summary_table, aes(y = reorder(gene, num_galectin_interactions), 
                          x = num_galectin_interactions)) +
  geom_bar(stat = "identity", fill = "grey", col = "black") +
  labs(
    # title = "Number of Galectins Interacting with Each Gene",
       y = NULL,
       x = NULL) +
  theme_classic(base_size = 16)
dev.off()