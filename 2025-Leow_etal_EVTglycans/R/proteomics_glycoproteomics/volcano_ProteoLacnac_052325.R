# Volcano plot of DeGlycoProteomics data 
# Author: Ke Leow
# Date: 03/02/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggrepel)

#--------------------------------
# Load data 
#--------------------------------
df <- read_csv('data/glycoproteomics/24-02-27_Angelo_Placenta_ProteomicsData.csv') 
lacnac <- read_csv('data/glycoproteomics/24-12-06_E00006_Angelo_LFQ_DB_DP_PL_Updated_PartialGalNew_021025.csv') %>%
  mutate(glycan.class = case_when(
    Partially.Galactosylated.new == TRUE ~ "Partially Galactosylated",
    Polylacnac == TRUE ~ "Polylacnac",
    Tetraantennary == TRUE ~ "Tetraantennary",
    Triantennary == TRUE ~ "Triantennary",
    Biantennary == TRUE ~ "Biantennary",
    `Glycan Classification` == "High Mannose" ~ "High Mannose",
    TRUE ~ "Other Complex/Hybrid"
  ))
lacnac_genes <- lacnac %>% 
  filter(glycan.class %in% c("Polylacnac","Tetraantennary","Partially Galactosylated")) %>% 
  pull(Gene) %>% unique()

#--------------------------------
# Pull comparison pair
#--------------------------------
input = "DB v.s. DP"

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NC"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$`Student's T-test Difference DP_DB` > 1 & df$`Student's T-test q-value DP_DB` < 0.05] <- "DP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$`Student's T-test Difference DP_DB` < -1 & df$`Student's T-test q-value DP_DB` < 0.05] <- "DB"
# flip position of DB and DP on x-axis (FC)
df$`Student's T-test Difference DB_DP` = -df$`Student's T-test Difference DP_DB`


# Filter to keep only the most significant entry per gene - and up in DB
# de_genes <- df %>%
#   group_by(Gene) %>%
#   slice_max(`-Log Student's T-test p-value DP_DB`, n = 1) %>%   # Keeps the entry with the lowest p-value for each gene
#   filter(diffexpressed == "DB") %>% 
#   ungroup() %>% pull(Genes)
# df$delabel <- ifelse(df$Genes %in% de_genes, df$Gene, NA)

# Find common genes (intersection)
common_genes <- intersect(lacnac_genes, df$Genes)

df$deglyco <- ifelse(df$Genes %in% common_genes, "deglyco", NA)
df$deglyco <- factor(df$deglyco)
df$alpha_value <- ifelse(df$Genes %in% common_genes, 1, 0.5)

top_genes_label <- df %>%
  filter(deglyco == "deglyco" & diffexpressed != "NC") %>%
  arrange(desc(`-Log Student's T-test p-value DP_DB`)) %>%
  slice_head(n = 5) %>%
  ungroup()

#volcano - color if detected in lacnac
pdf("R_plots/glycoproteomics/volcano_proteomicsMappedLacnacGlycoproteins_DBDP.pdf", , width = 3, height = 3.5)
ggplot(df, aes(x = `Student's T-test Difference DB_DP`, y = `-Log Student's T-test p-value DP_DB`, col = deglyco, alpha = alpha_value)) +
  geom_point()+
  theme_classic(base_size = 13)+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(legend.position = "none")+
  xlab("Log2(Fold Change)") + 
  ylab("-Log10(p-value)")+
  xlim(-7.5, 7.5) +  # Set equal limits around 0
  scale_color_manual(values = c("#FF7518", "grey")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = 1.3, col = "gray", linetype = 'dashed') + 
  geom_text_repel(
    data = top_genes_label,
    aes(label = Genes), color = "black",
    size = 4, max.overlaps = 10
  ) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank()        # Remove ticks
  )
dev.off()

#--------------------------------
# Check number of differential proteins in DB/DP
#--------------------------------
df <- df %>% 
  mutate(lacnac = ifelse(Genes %in% common_genes & diffexpressed != "NC", "lacnac", NA))

#how many lacnac glycoproteins differential in DB and DP
df %>% filter(lacnac == "lacnac") %>% group_by(diffexpressed) %>% summarise(n=n())

# how many proteins differential
df %>% group_by(diffexpressed) %>% summarise(n=n())

