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
#glycan+gene zscore data - with kmeans clusters
glycan_data <- read_csv("data/MALDI_IF_EVT_glycans/EVT_glycans_deg_kmeans_120624.csv") %>% 
  filter(Feature == "Glycan")
evt_glycans <- glycan_data %>% pull(gene_name)
  
glycopeptide_data <- read_csv("data/glycoproteomics/24-12-06_E00006_Angelo_LFQ_DB_DP_PL_Updated_PartialGalNew_021025.csv")

evt_hb <- read_csv("data/glycoproteomics/evt_HB_25_052725.csv") %>% pull(gene)

#select genes ordered in NS data - updated to combine FV+VCT Sept 2025
map_glyco_transcr_new <- read_csv("data/EVT_Nanostring/EVT_DEG_heatmap_DESeq2_090525.csv") %>% 
  rename(
    PV   = mean_by_stage_norm1,
    pEVT = mean_by_stage_norm2,
    iEVT = mean_by_stage_norm3,
    eEVT = mean_by_stage_norm4
  )%>% 
  subset(gene_name %in% evt_hb) 

#--------------------------------
# Plot heatmap - EVTA/I genes with NS data
#--------------------------------
#plot heatmap 
#using updated NS data
plotheatmap <- map_glyco_transcr_new %>% 
  # filter(k_cluster == 2 | k_cluster ==3)%>% 
  column_to_rownames("gene_name") %>% 
  select_at(vars(contains("V"))) 
# row.names(plotheatmap) <- map_glyco_transcr$gene_name

# Define color palette
color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# Define breaks ensuring 0 is centered
# Define symmetric breaks around zero
lim <- max(abs(plotheatmap))  # Find the max absolute value in data
breaks <- seq(-lim, lim, length.out = 101)  # Ensure zero is centered

# pdf("R_plots/glycoproteomics/heatmap_evt_HB_NS_selected_101925.pdf", width = 3, height = 3)
pheatmap(plotheatmap,
         # cellwidth = 30,
         # annotation_col = gal_anno,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE, cluster_rows = FALSE,
         angle_col = "90",
         border_color = "grey5",
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         breaks = breaks
)
dev.off()

#--------------------------------
# Get logFC and pval (DB vs DP treated) for EVT-HB glycopeptides
#--------------------------------
evt_hb_selected <- map_glyco_transcr_new %>% pull(gene_name)
#for duplicates, keep glycopeptide w largest treated_DB_vs_DP_limma_logFC
evt_hb_glycopeptide_filtered <- glycopeptide_data %>%
  filter(Gene %in% evt_hb_selected) %>% 
  group_by(Peptide.Sequence, composition) %>%
  filter((treated_DB_vs_DP_limma_logFC > 1 & `treated_DB_vs_DP_limma_adj.P.Val` < 0.05) | 
           (untreated_DB_vs_DP_limma_logFC > 1 & `untreated_DB_vs_DP_limma_adj.P.Val` < 0.05)) %>% #DB glycopeptides
  select(Gene, Peptide.Sequence, composition, treated_DB_vs_DP_limma_logFC, treated_DB_vs_DP_limma_adj.P.Val, Partially.Galactosylated.new) %>%
  ungroup() 

# 2. Pivot to make columns = glycans (ordered), rows = glycopeptides
heatmap_data <- evt_hb_glycopeptide_filtered %>%
  select(Gene, Peptide.Sequence, composition, treated_DB_vs_DP_limma_logFC) %>%
  pivot_wider(
    names_from = composition,
    values_from = treated_DB_vs_DP_limma_logFC,
    values_fill = list(treated_DB_vs_DP_limma_logFC = NA_real_)
  )

#get glycan annotations
glycan_anno <- read_csv("data/glycan_peaklist_paperAnnotations_032025.csv") %>% 
  mutate(Branching = Class_Branching)

#glycan+gene zscore data - with kmeans clusters
data <- read_csv("data/EVT_Nanostring/combined_gly_enz_heatmap_nopartialgal_type_list_K_4.csv") %>% 
  mutate(Feature = if_else(startsWith(gene_name, "H"), "Glycan", "Gene")) %>% 
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
  arrange(k_cluster_new)%>% filter(Feature == "Glycan") %>% 
  rename(composition = gene_name, cluster = k_cluster_new) %>% 
  full_join(glycan_anno) %>% select(composition, cluster, Sialyl, Fucose, Branching) 

glycan_anno_df <- evt_hb_glycopeptide_filtered %>% 
  mutate(Partial.Gal = ifelse(Partially.Galactosylated.new, TRUE, NA)) %>% 
  select(composition, Partial.Gal) %>% unique() %>% 
  left_join(.,data) %>% 
  column_to_rownames("composition") %>% 
  mutate(across(everything(), as.factor)) 

#--------------------------------
# Plot heatmap for each protein - sum glycopeptides
#--------------------------------
#get total glycan intensities
gene_totals <- heatmap_data %>%
  group_by(Gene) %>%
  summarise(across(-Peptide.Sequence, ~sum(.x, na.rm = TRUE)))

gene_matrix <- gene_totals %>%
  column_to_rownames("Gene") %>%
  as.matrix()
# #group heatmap by gene order
gene_matrix <- gene_matrix[evt_hb_selected, rownames(glycan_anno_df)]

# convert to logical matrix: TRUE if >0, FALSE if ==0
binary_matrix <- as.matrix(gene_matrix > 0)
binary_matrix_num <- binary_matrix * 1 # convert to numeric (0 = FALSE, 1 = TRUE)

glycan_anno_df$cluster <- as.factor(glycan_anno_df$cluster)
glycan_anno_df$Branching <- factor(glycan_anno_df$Branching, levels = c("High Mannose", "2","3","4","4P"))
glycan_anno_df$Sialyl <- as.factor(glycan_anno_df$Sialyl)
glycan_anno_df$Fucose <- as.factor(glycan_anno_df$Fucose)

# rename branching labels
levels(glycan_anno_df$Branching)[levels(glycan_anno_df$Branching) == "2"] <- "Bi"
levels(glycan_anno_df$Branching)[levels(glycan_anno_df$Branching) == "3"] <- "Tri"
levels(glycan_anno_df$Branching)[levels(glycan_anno_df$Branching) == "4"] <- "Tetra"
levels(glycan_anno_df$Branching)[levels(glycan_anno_df$Branching) == "4P"] <- "PolyLacNAc"

# Use the "YlOrRd" color palette from RColorBrewer
blue_colors <- brewer.pal(5, "Blues")
purple_colors <- brewer.pal(5, "Purples")
red_colors <- brewer.pal(4, "Reds")
# set3_colors <- brewer.pal(4, "Set3")

# Define custom annotation colors
annotation_colors <- list(
  cluster =c("1" = "#FF8811FF",  "3" = "#D44D5CFF", "4" = "#046E8FFF"),
  # Class = c(`Complex/Hybrid` = "cornflowerblue", `High Mannose` = "darkolivegreen1"),
  # Class_Branching = c(`Complex/Hybrid` = blue_colors[5], `High Mannose` = "darkolivegreen1"),
  Branching = c(`High Mannose` = "darkolivegreen1",`Bi` = blue_colors[2], `Tri` = blue_colors[3], `Tetra` = blue_colors[4], `PolyLacNAc` = blue_colors[5]),
  Sialyl = c(`0` = "white", `1` = purple_colors[3], `2` = "purple"),
  Fucose = c(`0` = "white", `1` = red_colors[2], `2` = red_colors[3], `3` = red_colors[4])
)

# pdf("R_plots/glycoproteomics/heatmap_evt_HB_NS_cluster4_glycoforms.pdf", width = 8, height = 4)
pheatmap(binary_matrix_num,
         annotation_col = glycan_anno_df,
         color = c("white", "black"),
         legend = FALSE,
         cluster_rows = FALSE,
         # cluster_cols = FALSE,
         angle_col = 90,
         border_color = "grey5",
         annotation_colors = annotation_colors
         # labels_col = glycan_anno_df$composition_label
         # gaps_col=c(17, 29, 59),
         # show_rownames=FALSE
)
dev.off()