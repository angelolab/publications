# Volcano plot of DE glycans
# Adapted from: https://biostatsquid.com/volcano-plots-r-tutorial/
# Author: Ke Leow
# Date: 03/17/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggrepel)

#--------------------------------
# Load data 
#--------------------------------
df <- read_csv("data/DBDP_transition/MALDI/library_matched/DE_pval_glycanTypes_DBDP.csv") %>% 
  #get -log10(pval)
  # mutate(neglog10_pval = -log10(DB_DP_P.Value)) %>% 
  mutate(neglog10_pval = -log10(DB_DP_adj.P.Val)) %>% 
  filter(composition != "highlybranched") %>% 
  mutate(composition = ifelse(composition == "fucosylated_sialylated", "fucosyl_sialyl", composition))

#match with EVT glycan type annotations
selected_types <- read_csv("data/EVT_Nanostring/Z_score_tab_nopartialgal_type_list_K_4.csv") %>% 
  pull(gly_type_nodes)

df_filter <- df %>% 
  filter(composition %in% selected_types)

#--------------------------------
# Plot Volcano
#--------------------------------
max_axis = max(df$DB_DP_logFC)

pdf("R_plots/DBDP_transition/volcano_DEglycanTypes_limma_alldonors_10types.pdf", width = 4, height = 4.5)
ggplot(df_filter, aes(x = DB_DP_logFC, y = neglog10_pval, label = composition)) +
  geom_point(size = 3, col = "grey20")+
  geom_text_repel(max.overlaps = Inf)+ # To show all labels 
  theme_classic(base_size = 13)+
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-max_axis, max_axis) +  # Set equal limits around 0
  labs(x = "Log2(Fold Change)",
       y = "-Log10(adj. p-value)"
  )+ 
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank()        # Remove ticks
  )
  # xlab("log2fc") + 
  # ylab("-log10(pval)")
# theme(legend.position = "none")+
  # scale_color_manual(values = c("#5C0B18", "#2D4A72", "grey")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  # geom_vline(xintercept = c(-logfc, logfc), col = "gray", linetype = 'dashed') +
  # geom_hline(yintercept = 1.5, col = "gray", linetype = 'dashed') +
  # geom_text_repel(max.overlaps = Inf) # To show all labels 
dev.off()

