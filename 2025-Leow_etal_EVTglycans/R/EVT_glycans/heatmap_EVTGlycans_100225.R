# Heatmap of EVT glycans and genes
# Author: Ke Leow
# Date: 03/20/25

#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

#--------------------------------
# Load data
#--------------------------------
#glycan annotations
glycan_anno <- read_csv("data/glycan_peaklist_paperAnnotations_032025.csv") %>% 
  mutate(Branching = Class_Branching)

#glycan+gene zscore data - with kmeans clusters
data <- read_csv("data/EVT_Nanostring/combined_gly_enz_heatmap_nopartialgal_type_list_K_4.csv") %>% 
  rename(
    PV   = mean_by_stage_norm1,
    pEVT = mean_by_stage_norm2,
    iEVT = mean_by_stage_norm3,
    eEVT = mean_by_stage_norm4
  ) %>% 
  mutate(k_cluster_new = case_when(
    k_cluster == 3 ~ 2,
    k_cluster == 2 ~ 3,
    TRUE ~ k_cluster
  )) %>%
  arrange(k_cluster_new)%>% 
  inner_join(glycan_anno, by = join_by(gene_name == composition)) %>% 
  column_to_rownames("gene_name")

#--------------------------------
# Plot heatmap for glycans
#--------------------------------
# data_sub <- data %>% filter(Feature == "Glycan")
data_int <- data %>% select_at(vars(contains("V", ignore.case=FALSE))) %>% t()
data_pheno <- data %>% select(Sialyl, Fucose, Branching)

# data_pheno$Cluster <- as.factor(data_pheno$Cluster)
data_pheno$Branching <- factor(data_pheno$Branching, levels = c("High Mannose", "2","3","4","4P"))
# data_pheno$Branches <- as.factor(data_pheno$Branches)
data_pheno$Sialyl <- as.factor(data_pheno$Sialyl)
data_pheno$Fucose <- as.factor(data_pheno$Fucose)

# rename branching labels
levels(data_pheno$Branching)[levels(data_pheno$Branching) == "2"] <- "Bi"
levels(data_pheno$Branching)[levels(data_pheno$Branching) == "3"] <- "Tri"
levels(data_pheno$Branching)[levels(data_pheno$Branching) == "4"] <- "Tetra"
levels(data_pheno$Branching)[levels(data_pheno$Branching) == "4P"] <- "PolyLacNAc"

# Use the "YlOrRd" color palette from RColorBrewer
blue_colors <- brewer.pal(5, "Blues")
purple_colors <- brewer.pal(5, "Purples")
red_colors <- brewer.pal(4, "Reds")
# set3_colors <- brewer.pal(4, "Set3")

# Define custom annotation colors
annotation_colors <- list(
  # Cluster = c(`1` = set3_colors[1], `2` = set3_colors[2],`3` = set3_colors[3],`4` = set3_colors[4]), 
  # Class = c(`Complex/Hybrid` = "cornflowerblue", `High Mannose` = "darkolivegreen1"),
  # Class_Branching = c(`Complex/Hybrid` = blue_colors[5], `High Mannose` = "darkolivegreen1"),
  Branching = c(`High Mannose` = "darkolivegreen1",`Bi` = blue_colors[2], `Tri` = blue_colors[3], `Tetra` = blue_colors[4], `PolyLacNAc` = blue_colors[5]),
  Sialyl = c(`0` = "white", `1` = purple_colors[3], `2` = "purple"),
  Fucose = c(`0` = "white", `1` = red_colors[2], `2` = red_colors[3], `3` = red_colors[4])
)

#get cluster breaks by cumulative sum of number genes in each cluster e.g.c(10,20,30,40)
data %>% group_by(k_cluster_new) %>% summarise(n=n()) %>% pull(n)
# compute where cluster changes
gaps_col <- cumsum(rle(data$k_cluster_new)$lengths)
gaps_col <- gaps_col[-length(gaps_col)]

# Define breaks ensuring 0 is centered
# Define symmetric breaks around zero
lim <- max(abs(data_int))  # Find the max absolute value in data
# lim = 2 #manual input
breaks <- seq(-lim, lim, length.out = 101)  # Ensure zero is centered

pdf("R_plots/MALDI_IF_EVT_glycans/heatmap_EVTglycans_wGlyLabel_100225.pdf", , width = 11.5, height = 4)
pheatmap(data_int,
         annotation_col = data_pheno,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         angle_col = 90,
         # show_colnames=FALSE,
         border_color = "grey5",
         gaps_col=gaps_col,
         breaks = breaks,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         annotation_colors = annotation_colors)
dev.off()

# #--------------------------------
# # Plot heatmap for genes
# #--------------------------------
# data_sub <- data %>% filter(Feature == "Gene")
# 
# data_int <- data_sub %>% select_at(vars(contains("V", ignore.case=FALSE))) 
# data_pheno <- data_sub %>% select(Cluster)
# data_pheno$Cluster <- as.factor(data_pheno$Cluster)
# 
# pheatmap(data_int,
#          annotation_row = data_pheno,
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          # show_rownames=FALSE,
#          color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                    "RdBu")))(100))
# 
