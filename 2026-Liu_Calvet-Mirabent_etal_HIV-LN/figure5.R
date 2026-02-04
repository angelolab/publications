# Code for Figure 5, Supplementary Figure 9
# Author: Candace Liu

library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(nlme)
library(emmeans)
library(VennDiagram)
library(grid)

cell_tab = fread("../data/tables/cell_table_size_normalized.csv")
functional_tab = fread("../data/tables/nimbus_binarized.csv")
feature_tab = fread("../data/tables/feature_tab.csv")
cd8t_me_cell_tab = fread("../data/spatial_analysis/neighborhood_analysis_cd8t/cell_meta_cluster_radius50_k4_counts/cell_table_size_normalized_kmeans_nh.csv")
tissue_area = fread("../data/tables/tissue_area.csv")
metadata = fread("../data/tables/metadata.csv")


# Function for getting pvalues between viremia groups
pval_viremia_groups <- function(feature, dat) {
  formula = as.formula(paste(feature,"~status_with_viremia"))
  model = lme(formula, random= ~1|sample_id, data = dat, na.action=na.omit)
  emm = emmeans(model, "status_with_viremia")
  
  pairwise_results = pairs(emm, adjust = "tukey")
  pairwise_df = as.data.frame(pairwise_results)
  return(pairwise_df)
}



### Fig 5a: Heatmap of cell densities across viremia groups
keep_features = c("B_density", "CD4T_density","CD8T_density","CD11c_density","CD11c_CD14_density","CD11c_CD68_density")
compare = feature_tab[, lapply(.SD, median), by = .(status_with_viremia), .SDcols = keep_features]
compare_norm = compare[, lapply(.SD, function(x) x / quantile(x[x!=0], 0.999)), .SDcols = keep_features]
compare_norm = cbind(compare[,c("status_with_viremia")], compare_norm)
compare_dt = data.frame(compare_norm[,..keep_features])
rownames(compare_dt) = compare_norm$status_with_viremia

pdf("statuswithviremia_cell_density_heatmap.pdf", height=3, width=6)
pheatmap(compare_dt,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#512DA8"))(100),
         breaks = seq(0, 1, length.out = 101))
dev.off()



### Fig 5b: Functional markers in all CD11c+ cells
cd11c_celltypes = c("CD11c","CD11c_CD14","CD11c_CD68")
keep_markers = c("Caspase1","NLRP3","Galectin9","CD86")
keep_functional_tab = functional_tab[cell_meta_cluster %in% cd11c_celltypes]
cd11c_counts = keep_functional_tab[,.N,by=.(fov)]

functional_counts = keep_functional_tab[, lapply(.SD, sum, na.rm = TRUE), by = fov, .SDcols = keep_markers]
functional_counts = functional_counts[cd11c_counts[,c("fov","N")], on=.(fov)]
functional_counts[, (keep_markers) := lapply(.SD, function(x) x / N), .SDcols = keep_markers]
functional_counts_melt = melt(functional_counts, id.vars=c("fov"), measure.vars=keep_markers, variable.name="marker", value.name="prop")
functional_counts_melt = functional_counts_melt[metadata, on=.(fov)]

ggplot(functional_counts_melt, aes(x=status_with_viremia, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  facet_wrap(~marker, scales = "free_y") +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of all CD11c+ cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_blank())
ggsave("cd11c_functional_statuswithviremia.pdf", height=5, width=5)

# One marker (with pvalue)
one_marker = "CD86"
one_marker_dat = functional_counts_melt[marker==one_marker]
pairwise_df_out = pval_viremia_groups("prop", one_marker_dat)
comparisons = data.frame(
  comparison = pairwise_df_out$contrast,
  p_value = pairwise_df_out$p.value,
  estimate = pairwise_df_out$estimate,
  p_formatted = ifelse(pairwise_df_out$p.value < 0.001, "p < 0.001",
                       ifelse(pairwise_df_out$p.value < 0.01, paste("p =", round(pairwise_df_out$p.value, 3)),
                              paste("p =", round(pairwise_df_out$p.value, 2))))
)
comp_list = list()
for(i in 1:nrow(pairwise_df_out)) {
  groups = trimws(strsplit(pairwise_df_out$contrast[i], " - ")[[1]])
  comp_list[[i]] = groups
}
ggplot(one_marker_dat, aes(x=status_with_viremia, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of all CD11c+ cells",
       title = one_marker) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(
    data.frame(
      group1 = sapply(comp_list, `[`, 1),
      group2 = sapply(comp_list, `[`, 2),
      p.adj = pairwise_df_out$p.value,
      p.adj.signif = comparisons$p_formatted
    ),
    label = "p.adj.signif",
    y.position = c(max(one_marker_dat$prop, na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )



### Fig 5c: NLRP3+Caspase1+ cells
subset_cells = functional_tab[Caspase1==1 & NLRP3==1]
cell_counts = subset_cells[,.N,by=.(fov)]
cell_counts = cell_counts[tissue_area, on=.(fov)]
cell_counts[is.na(cell_counts)] = 0
cell_counts[,density := N/tissue_area/0.64]
cell_counts = cell_counts[metadata, on=.(fov)]
ggplot(cell_counts, aes(x=status_with_viremia, y=density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = "NLRP3+Caspase1+ cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("caspase1_nlrp3_statuswithviremia.pdf", height=4, width=3)

# Add p-value
pairwise_df_out = pval_viremia_groups("density", cell_counts)
comparisons = data.frame(
  comparison = pairwise_df_out$contrast,
  p_value = pairwise_df_out$p.value,
  estimate = pairwise_df_out$estimate,
  p_formatted = ifelse(pairwise_df_out$p.value < 0.001, "p < 0.001",
                       ifelse(pairwise_df_out$p.value < 0.01, paste("p =", round(pairwise_df_out$p.value, 3)),
                              paste("p =", round(pairwise_df_out$p.value, 2))))
)
comp_list = list()
for(i in 1:nrow(pairwise_df_out)) {
  groups = trimws(strsplit(pairwise_df_out$contrast[i], " - ")[[1]])
  comp_list[[i]] = groups
}
ggplot(cell_counts, aes(x=status_with_viremia, y=density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = "NLRP3+Caspase1+ cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(
    data.frame(
      group1 = sapply(comp_list, `[`, 1),
      group2 = sapply(comp_list, `[`, 2),
      p.adj = pairwise_df_out$p.value,
      p.adj.signif = comparisons$p_formatted
    ),
    label = "p.adj.signif",
    y.position = c(max(cell_counts$density, na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )
ggsave("caspase1_nlrp3_statuswithviremia_withpval.pdf", height=4, width=3)



### Supp Fig 9a-b: Boxplot of features compared between viremia groups
feature = "B_density"
ggplot(feature_tab, aes(x=status_with_viremia, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = feature) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithviremia.pdf"), height=4, width=3.75)

# Add p-value
pairwise_df_out = pval_viremia_groups(feature, feature_tab)
comparisons = data.frame(
  comparison = pairwise_df_out$contrast,
  p_value = pairwise_df_out$p.value,
  estimate = pairwise_df_out$estimate,
  p_formatted = ifelse(pairwise_df_out$p.value < 0.001, "p < 0.001",
                       ifelse(pairwise_df_out$p.value < 0.01, paste("p =", round(pairwise_df_out$p.value, 3)),
                              paste("p =", round(pairwise_df_out$p.value, 2))))
)
comp_list = list()
for(i in 1:nrow(pairwise_df_out)) {
  groups = trimws(strsplit(pairwise_df_out$contrast[i], " - ")[[1]])
  comp_list[[i]] = groups
}
ggplot(feature_tab, aes(x=status_with_viremia, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = feature) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(
    data.frame(
      group1 = sapply(comp_list, `[`, 1),
      group2 = sapply(comp_list, `[`, 2),
      p.adj = pairwise_df_out$p.value,
      p.adj.signif = comparisons$p_formatted
    ),
    label = "p.adj.signif",
    y.position = c(max(feature_tab[[feature]], na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )
ggsave(paste0(feature, "_statuswithviremia_withpval.pdf"), height=4, width=3.75)


feature = "Ki67pos_B_prop"
ggplot(feature_tab, aes(x=status_with_viremia, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of B cells",
       title = "Ki67+ B cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithviremia.pdf"), height=4, width=3.75)

# Get p-value
pairwise_df_out = pval_viremia_groups(feature, feature_tab)
comparisons = data.frame(
  comparison = pairwise_df_out$contrast,
  p_value = pairwise_df_out$p.value,
  estimate = pairwise_df_out$estimate,
  p_formatted = ifelse(pairwise_df_out$p.value < 0.001, "p < 0.001",
                       ifelse(pairwise_df_out$p.value < 0.01, paste("p =", round(pairwise_df_out$p.value, 3)),
                              paste("p =", round(pairwise_df_out$p.value, 2))))
)
comp_list = list()
for(i in 1:nrow(pairwise_df_out)) {
  groups = trimws(strsplit(pairwise_df_out$contrast[i], " - ")[[1]])
  comp_list[[i]] = groups
}
ggplot(feature_tab, aes(x=status_with_viremia, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of B cells",
       title = "Ki67+ B cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(
    data.frame(
      group1 = sapply(comp_list, `[`, 1),
      group2 = sapply(comp_list, `[`, 2),
      p.adj = pairwise_df_out$p.value,
      p.adj.signif = comparisons$p_formatted
    ),
    label = "p.adj.signif",
    y.position = c(max(feature_tab[[feature]], na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )
ggsave(paste0(feature, "_statuswithviremia_withpval.pdf"), height=4, width=3.75)



### Supp Fig 9c: CD8T ME distribution between viremia groups
nh_counts = cd8t_me_cell_tab[,.N,by=.(kmeans_neighborhood, fov)]
nh_counts = nh_counts[!is.na(kmeans_neighborhood)]
nh_counts_wide = dcast(nh_counts, fov ~ kmeans_neighborhood, value.var = "N", fill = 0)
nh_counts = melt(nh_counts_wide, id.vars="fov", value.name = "N", variable.name = "kmeans_neighborhood")

all_cd8t_counts = cd8t_me_cell_tab[, .(total_cd8t=.N), by=.(fov)]
nh_counts = nh_counts[all_cd8t_counts, on=.(fov)]
nh_counts[, prop := N/total_cd8t]
nh_counts = nh_counts[metadata, on=.(fov)]

nh_num = 4
nh_counts_sub = nh_counts[kmeans_neighborhood==nh_num]
ggplot(nh_counts_sub, aes(x=status_with_viremia, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  labs(y="Proportion of CD8+ T cells") +
  theme_bw() +
  labs(title = paste0("Cluster ", nh_num)) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0("cd8t_me",nh_num,"_proportions_statuswithviremia.pdf"), height=4, width=3.75)

# Add p-value
pairwise_df_out = pval_viremia_groups("prop", nh_counts_sub)
comparisons = data.frame(
  comparison = pairwise_df_out$contrast,
  p_value = pairwise_df_out$p.value,
  estimate = pairwise_df_out$estimate,
  p_formatted = ifelse(pairwise_df_out$p.value < 0.001, "p < 0.001",
                       ifelse(pairwise_df_out$p.value < 0.01, paste("p =", round(pairwise_df_out$p.value, 3)),
                              paste("p =", round(pairwise_df_out$p.value, 2))))
)
comp_list = list()
for(i in 1:nrow(pairwise_df_out)) {
  groups = trimws(strsplit(pairwise_df_out$contrast[i], " - ")[[1]])
  comp_list[[i]] = groups
}
ggplot(nh_counts_sub, aes(x=status_with_viremia, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  labs(y="Proportion of CD8+ T cells") +
  theme_bw() +
  labs(title = paste0("Cluster ", nh_num)) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(
    data.frame(
      group1 = sapply(comp_list, `[`, 1),
      group2 = sapply(comp_list, `[`, 2),
      p.adj = pairwise_df_out$p.value,
      p.adj.signif = comparisons$p_formatted
    ),
    label = "p.adj.signif",
    y.position = c(max(feature_tab[[feature]], na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )
ggsave(paste0("cd8t_me",nh_num,"_proportions_statuswithviremia_withpval.pdf"), height=4, width=3.75)



### Supp Fig 9e: Venn diagram of NLRP3+Caspase1+ cells
v = venn.diagram(
  x = list(
    Caspase1 = which(functional_tab$Caspase1 == 1),
    NLRP3 = which(functional_tab$NLRP3 == 1)
  ),
  filename = NULL,
  category.names = c("Caspase1", "NLRP3"),
  fill = c("lightblue", "lightgreen"),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.pos = 20,     # Position categories at 0 degrees
  cat.dist = 0.025,  # Distance of category names from the edge
  euler.d = TRUE,  # Use Euler diagram style
  scaled = TRUE 
)
grid.newpage()
grid.draw(v)

total_num_inflammasome = nrow(functional_tab[Caspase1==1 | NLRP3==1])

# No labels
v = venn.diagram(
  x = list(
    Caspase1 = which(functional_tab$Caspase1 == 1),
    NLRP3 = which(functional_tab$NLRP3 == 1)
  ),
  filename = "inflammasome_venn_nolabels.pdf",
  category.names = c("Caspase1", "NLRP3"),
  fill = c("lightblue", "lightgreen"),
  cat.cex = 0,       # Category text size = 0 (hides category names)
  cex = 0,          # Area label text size = 0 (hides overlap numbers)
  label.col = "transparent"
)



### Supp Fig 9f: Breakdown of NLRP3+Caspase1+ cells
keep_phenos = c("CD11c_CD68","CD11c","CD14_CD68_CD163","CD11c_CD14")
subset_cells = functional_tab[Caspase1==1 & NLRP3==1]
cell_counts = subset_cells[,.N,by=.(cell_meta_cluster)]
total_counts = functional_tab[,.(total_count = .N),by=.(cell_meta_cluster)]
cell_counts = total_counts[cell_counts, on=.(cell_meta_cluster)]
cell_counts[, prop := N/total_count]
ggplot(cell_counts[cell_meta_cluster %in% keep_phenos], aes(x=reorder(cell_meta_cluster, -prop), y=prop)) +
  geom_bar(stat="identity") +
  labs(y="Proportion of cell type",
       title="NLRP3+Caspase1+ cells") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("caspase1_nlrp3_bar_prop.pdf", height=4, width=2.5)



### Supp Fig 9g: Comparison of NLRP3+Caspase1+ cells in all CD11c+ cells
cd11c_celltypes = c("CD11c","CD11c_CD14","CD11c_CD68")
cd11c_cell_tab = functional_tab[cell_meta_cluster %in% cd11c_celltypes]
cd11c_counts = cd11c_cell_tab[,.(total_cd11c = .N),by=.(fov)]

subset_cells = cd11c_cell_tab[Caspase1==1 & NLRP3==1]
cell_counts = subset_cells[,.N,by=.(fov)]
cell_counts = cell_counts[cd11c_counts, on=.(fov)]
cell_counts[is.na(cell_counts)] = 0
cell_counts[,prop := N/total_cd11c]
cell_counts = cell_counts[metadata, on=.(fov)]
ggplot(cell_counts, aes(x=status_with_viremia, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of all CD11c+ cells",
       title = "NLRP3+Caspase1+ in all CD11c+ cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("caspase1_nlrp3_cd11c_statuswithviremia.pdf", height=4, width=3.75)



### Supp Fig 9h: Breakdown of Caspase1+ only cells
subset_cells = functional_tab[Caspase1==1 & NLRP3==0]
# Only keep top 6
top_phenos = names(sort(table(subset_cells$cell_meta_cluster), decreasing = TRUE)[1:6])
subset_cells_top = subset_cells[subset_cells$cell_meta_cluster %in% top_phenos, ]
ggplot(subset_cells_top, aes(x = reorder(cell_meta_cluster, cell_meta_cluster, function(x) -length(x)))) +
  geom_bar() +
  labs(y = "Number of cells",
       title = "Caspase1+ only") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("caspase1_only_bar.pdf", height=4, width=3)



### Supp Fig 9i: Caspase1+ CD8T+ between viremia groups
subset_cells = functional_tab[cell_meta_cluster=="CD8T"]
cd8t_counts = subset_cells[,.(total_cd8t = .N), by=.(fov)]
cell_counts = subset_cells[Caspase1==1 & NLRP3==0, .N, by=.(fov)]
cell_counts = cell_counts[cd8t_counts, on=.(fov)]
cell_counts[is.na(cell_counts)] = 0
cell_counts = metadata[cell_counts, on=.(fov)]
cell_counts[, prop := N/total_cd8t]
ggplot(cell_counts, aes(x=status_with_viremia, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of CD8+ T cells",
       title = "Caspase1+ only") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave("Caspase1only_CD8T_proportion_statuswithviremia.pdf", height=5, width=4)


