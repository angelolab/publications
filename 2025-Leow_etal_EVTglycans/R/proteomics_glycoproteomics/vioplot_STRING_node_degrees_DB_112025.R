# Map DE glycoproteins to Andi's transcriptomics data
# Author: Ke Leow
# Date: 02/18/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggpubr)

#--------------------------------
# DB-HB vs DB-other
#--------------------------------
## 1. Read in your STRING node-degree tables
decB_HB <- read_tsv("STRING_interactions/44genes_DB_HB_treated_untreated/string_node_degrees.tsv")
decB_other <- read_tsv("STRING_interactions/104genes_DB_allexcept HB/string_node_degrees.tsv")

## 2. Combine into one long dataframe
## Replace `degree` below with the actual column name that contains node degree
all_deg <- bind_rows(
  decB_HB   %>% transmute(group = "decB-HB",    degree = node_degree),
  decB_other %>% transmute(group = "decB-other", degree = node_degree)
) %>%
  mutate(
    group = factor(group, levels = c("decB-other", "decB-HB"))  # order on x-axis
  )

## 3. Violin + boxplot + stats
ggplot(all_deg, aes(x = group, y = degree)) +
  geom_violin(trim = FALSE, alpha = 0.5, width = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  stat_compare_means(
    method = "t.test",        # or "t.test" if you prefer
    label = "p.format",            # show p-value in formatted form
    comparisons = list(c("decB-other", "decB-HB"))
  ) +
  labs(
    x = NULL,
    y = "STRING node degree",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

#--------------------------------
# DB-HB vs DB-other vs DP-HB vs DP-other
#--------------------------------
decP_other <- read_tsv("STRING_interactions/94genes_DPglyco_minusHB/string_node_degrees.tsv")
decP_HB    <- read_tsv("STRING_interactions/58genes_DP_HB/string_node_degrees.tsv")

## 2. Combine into long df
all_deg <- bind_rows(
  decB_HB    %>% transmute(group = "decB-HB",    degree = node_degree),
  decB_other %>% transmute(group = "decB-other", degree = node_degree),
  decP_HB    %>% transmute(group = "decP-HB",    degree = node_degree),
  decP_other %>% transmute(group = "decP-other", degree = node_degree)
) %>% 
  mutate(group = factor(group, levels = c("decP-other", "decP-HB", "decB-other", "decB-HB")))

## 3. All pairwise t-test comparisons (customize if needed)
my_comparisons <- list(
  c("decP-HB", "decP-other"),
  c("decB-other", "decP-other"),
  c("decB-HB", "decP-other"),
  c("decB-HB", "decP-HB"),
  c("decB-HB", "decB-other")
)

## 4. Violin + boxplot + t-test significance
pdf("R_plots/glycoproteomics/vioplot_STRING_numEdges_DBDP_HB_other.pdf", , width = 4, height = 3)
ggplot(all_deg, aes(x = group, y = degree)) +
  geom_violin(trim = FALSE, alpha = 0.5, width = 0.8) +
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.85) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif",
    vjust = 0.5,
    hide.ns = TRUE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
dev.off()
