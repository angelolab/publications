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
glycan_anno <- read_csv("data/MALDI_IF_EVT_glycans/glycan_types_enzymes_051524.csv") %>% 
  rename(gene_name = 1) %>% 
  mutate(Fucose = ifelse(fucosylated == 1, "Fucose", NA),
         Sialyl = ifelse(sialylated == 1, "Sialyl", NA)) %>% 
  mutate(Branching = case_when(
    highMannose == 1 ~ "High Mannose",
    # bisecting == 1 ~ "Bisecting",
    biantennary == 1  ~ "Bi",
    triantennary == 1 | tetraantennary == 1  ~ "Tri/Tetra",
    polylacnac == 1 ~ "PolyLacNAc",
    TRUE ~ NA_character_
    ) 
  ) %>% 
  mutate(Branching = case_when(
    str_detect(gene_name, "GALT") ~ NA_character_,
    TRUE ~ Branching
  ))%>%
  mutate(gene_name = str_replace(gene_name, "MGAT4a", "MGAT4A"))

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
  inner_join(glycan_anno, by = join_by(gene_name)) %>% 
  column_to_rownames("gene_name")

#--------------------------------
# Plot heatmap for enzymes
#--------------------------------
# data_sub <- data %>% filter(Feature == "Glycan")
data_int <- data %>% select_at(vars(contains("V", ignore.case=FALSE))) %>% t()
data_pheno <- data %>% select(Sialyl, Fucose, Branching)

# data_pheno$Cluster <- as.factor(data_pheno$Cluster)
# data_pheno$Class <- as.factor(data_pheno$Class)
data_pheno$Branching <- factor(data_pheno$Branching, levels = c("High Mannose", "Bi","Tri/Tetra","PolyLacNAc"))
data_pheno$Sialyl <- as.factor(data_pheno$Sialyl)
data_pheno$Fucose <- as.factor(data_pheno$Fucose)

# Use the "YlOrRd" color palette from RColorBrewer
blue_colors <- brewer.pal(5, "Blues")
purple_colors <- brewer.pal(5, "Purples")
red_colors <- brewer.pal(4, "Reds")
# set3_colors <- brewer.pal(4, "Set3")

# Define custom annotation colors
annotation_colors <- list(
  # Cluster = c(`1` = set3_colors[1], `2` = set3_colors[2],`3` = set3_colors[3],`4` = set3_colors[4]), 
  # Class = c(`Complex/Hybrid` = "cornflowerblue", `High Mannose` = "darkolivegreen1"),
  # Class = c(`Complex/Hybrid` = blue_colors[5], `High Mannose` = "darkolivegreen1"),
  # Branching = c(`2` = blue_colors[2], `3` = blue_colors[3], `4` = blue_colors[4], `4P` = blue_colors[5]),
  Branching = c(`High Mannose` = "darkolivegreen1",`Bi` = blue_colors[2], `Tri/Tetra` = blue_colors[4], `PolyLacNAc` = blue_colors[5]),
  Sialyl = c(`Sialyl` = purple_colors[5]),
  Fucose = c(`Fucose` = red_colors[4])
)

#get cluster breaks by cumulative sum of number genes in each cluster e.g.c(10,20,30,40)
data %>% group_by(Cluster) %>% summarise(n=n()) %>% pull(n)
# compute where cluster changes
gaps_col <- cumsum(rle(data$k_cluster_new)$lengths)
gaps_col <- gaps_col[-length(gaps_col)]

# Define breaks ensuring 0 is centered
# Define symmetric breaks around zero
lim <- max(abs(data_int))  # Find the max absolute value in data
breaks <- seq(-lim, lim, length.out = 101)  # Ensure zero is centered

pdf("R_plots/MALDI_IF_EVT_glycans/heatmap_EVTenzymes_100225.pdf", , width = 5.5, height = 3)
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
