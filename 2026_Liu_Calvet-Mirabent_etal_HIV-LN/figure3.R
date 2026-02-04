# Code for Figure 3, Supplementary Figure 6
# Author: Candace Liu

library(data.table)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(nlme)
library(emmeans)

cell_tab = fread("../data/tables/cell_table_size_normalized.csv")
functional_tab = fread("../data/tables/nimbus_binarized.csv")
tissue_area = fread("../data/tables/tissue_area.csv")
metadata = fread("../data/tables/metadata.csv")
cd8t_me_cell_tab = fread("../data/spatial_analysis/neighborhood_analysis_cd8t/cell_meta_cluster_radius50_k4_counts/cell_table_size_normalized_kmeans_nh.csv")
cd8t_me_counts_tab = fread("../data/spatial_analysis/neighborhood_analysis_cd8t/cell_meta_cluster_radius50_k4_counts/nh_input_features.csv")
cd8t_follicle_dist_tab = fread("../data/tables/follicle_edge_distance.csv")
quiche_cell_tab = fread("../data/spatial_analysis/quiche/quiche_cell_tab.csv")
quiche_sig_niches = fread("../data/spatial_analysis/quiche/quiche_sig_niches_median.csv")
cd8t_me_colors = fread("../data/tables/kmeans_nh_colors_k4.csv")
enrichment_score_table_dir = "../data/spatial_analysis/cell_cell_enrichment_20um/tables"
enrichment_score_table_func_dir = "../data/spatial_analysis/cell_cell_distances_functional"

hivpos_niches = quiche_sig_niches[logFC > 0]$quiche_niche_neighborhood
hivneg_niches = quiche_sig_niches[logFC < 0]$quiche_niche_neighborhood

# Function for getting pvalues between p24 groups
pval_p24_groups <- function(feature, dat) {
  formula = as.formula(paste(feature,"~status_with_p24"))
  model = lme(formula, random= ~1|sample_id, data = dat, na.action=na.omit)
  emm = emmeans(model, "status_with_p24")
  
  pairwise_results = pairs(emm, adjust = "tukey")
  pairwise_df = as.data.frame(pairwise_results)
  return(pairwise_df)
}



### Supp Fig 6b: CD8T enrichment scores
one_file_name = "CD8T_APC"
one_tab = fread(file.path(enrichment_score_table_dir, paste0(one_file_name,".csv")))
one_tab = one_tab[metadata, on=.(fov)]
one_tab_melt = melt(one_tab,
                    id.vars = c("pheno1", "pheno2", "fov", "z", "sample_id"),
                    measure.vars = c("status", "status_with_viremia", "status_with_p24"),
                    variable.name = "status_type",
                    value.name = "group")
ggplot(one_tab_melt[status_type=="status_with_p24"], aes(x=group, y=z)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.6) +
  ggtitle(paste(unique(one_tab_melt$pheno1), ":", unique(one_tab_melt$pheno2))) +
  theme_bw() +
  labs(y = "Enrichment score") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0("es_",one_file_name,".pdf")), height=3, width=2.5)



### Fig 3c: Count of pos/neg QUICHE niches per group
hivpos_quiche_cell_tab = quiche_cell_tab[quiche_niche_neighborhood %in% hivpos_niches]
hivpos_niche_count = hivpos_quiche_cell_tab[,.N,by=.(fov)]
hivpos_niche_count = metadata[hivpos_niche_count, on=.(fov)]
hivpos_niche_count = tissue_area[hivpos_niche_count, on=.(fov)]
hivpos_niche_count[, density := N/tissue_area]
ggplot(hivpos_niche_count, aes(x=status_with_p24, y=density/0.64)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1) +
  theme_bw() +
  labs(y="Cell density (cells/mm2)",
       title="Significant niches in HIV+ group") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("quiche_hivpos.pdf", height=4, width=4)

# Add p-value
pairwise_df_out = pval_p24_groups("density", hivpos_niche_count)
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
ggplot(hivpos_niche_count, aes(x=status_with_p24, y=density/0.64)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = "Significant niches in HIV+ group") +
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
    y.position = c(max(hivpos_niche_count$density/0.64, na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )
ggsave("quiche_hivpos_withpval.pdf", height=4, width=4)

hivneg_quiche_cell_tab = quiche_cell_tab[quiche_niche_neighborhood %in% hivneg_niches]
hivneg_niche_count = hivneg_quiche_cell_tab[,.N,by=.(fov)]
hivneg_niche_count = metadata[hivneg_niche_count, on=.(fov)]
hivneg_niche_count = tissue_area[hivneg_niche_count, on=.(fov)]
hivneg_niche_count[, density := N/tissue_area]
ggplot(hivneg_niche_count, aes(x=status_with_p24, y=density/0.64)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1) +
  theme_bw() +
  labs(y="Cell density (cells/mm2)",
       title="Significant niches in HIV- group") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("quiche_hivneg.pdf", height=4, width=4)



### Supp Fig 6e: Count of pos QUICHE niches - split up p24+ into individual donors
hivpos_quiche_cell_tab = quiche_cell_tab[quiche_niche_neighborhood %in% hivpos_niches]
hivpos_niche_count = hivpos_quiche_cell_tab[,.N,by=.(fov)]
hivpos_niche_count = metadata[hivpos_niche_count, on=.(fov)]
hivpos_niche_count = tissue_area[hivpos_niche_count, on=.(fov)]
hivpos_niche_count[, density := N/tissue_area]

hivpos_niche_count[, group := status_with_p24]
hivpos_niche_count[grepl("slide3", fov), group := "hiv_pos_p24pos_diagnosis"]
hivpos_niche_count[grepl("slide4", fov), group := "hiv_pos_p24pos_stoppedart"]

ggplot(hivpos_niche_count, aes(x=group, y=density/0.64)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.6) +
  theme_bw() +
  labs(y="Cell density (cells/mm2)",
       title="Significant niches in HIV+ group") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("quiche_hivpos_statuswithp24_split2donors.pdf", height=3, width=3)



### Fig 3d: Functional marker differences between quiches
keep_markers = c("Caspase1","CD45RO","HLADR","Ki67","GranzymeB")
func_quiche_cell_tab = functional_tab[quiche_cell_tab, on=.(fov, label, cell_meta_cluster)]
func_quiche_cell_tab[, group := "other"]
func_quiche_cell_tab[quiche_niche_neighborhood %in% hivpos_niches, group := "hivpos_sig"]
func_quiche_cell_tab[quiche_niche_neighborhood %in% hivneg_niches, group := "hivneg_sig"]

cd8_only = func_quiche_cell_tab[cell_meta_cluster=="CD8T"]
func_counts = cd8_only[, lapply(.SD, sum, na.rm = TRUE), by = .(group), .SDcols = keep_markers]
num_total = cd8_only[,.N,by=.(group)]
func_counts = func_counts[num_total, on=.(group)]
func_counts[, (keep_markers) := lapply(.SD, function(x) x/N), .SDcols = keep_markers]
func_counts[, group := factor(group, levels = c("hivpos_sig","hivneg_sig", "other"))]
func_counts = func_counts[order(group)]
func_counts_norm = func_counts[, lapply(.SD, function(x) x / quantile(x[x!=0], 0.999)), .SDcols = keep_markers]
df = data.frame(func_counts_norm[,..keep_markers])
rownames(df) = func_counts$group

pdf("quiche_cd8t_func.pdf", height=4, width=4)
pheatmap(df,
         cluster_rows = FALSE,
         color = colorRampPalette(c("white", "#512DA8"))(100),
         breaks = seq(0.5, 1, length.out = 101))
dev.off()



### Supp Fig 6f: Caspase1-/+ CD8T distances
file00 = fread(file.path(enrichment_score_table_func_dir,"CD8T_0_CD8T_0_Caspase1.csv"))
file00[, group := paste0(pheno1, '_', pheno2)]
file01 = fread(file.path(enrichment_score_table_func_dir,"CD8T_0_CD8T_1_Caspase1.csv"))
file01[, group := paste0(pheno1, '_', pheno2)]
file11 = fread(file.path(enrichment_score_table_func_dir,"CD8T_1_CD8T_1_Caspase1.csv"))
file11[, group := paste0(pheno1, '_', pheno2)]

combined = rbindlist(list(file00, file01, file11))
combined = combined[metadata, on=.(fov)]

ggplot(combined, aes(x=group, y=mean*(800/2048))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(y = "Mean distance (um)",
       title = "Mean distance between Caspase1+/- CD8+ T cells") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("caspase1_cd8t_distances.pdf", height=5, width=5)



### Fig 3e: Heatmaps of CD8T MEs
# Blue-white-red colors
rwb_cols = colorRampPalette(c("royalblue4","white","red4"))(99)

## Cell composition (input features)
cd8t_me_counts_tab[, Other := Other + Immune_other + Unassigned]
cd8t_me_counts_tab[, kmeans_neighborhood := paste0("Cluster",kmeans_neighborhood)]
cluster_names = cd8t_me_counts_tab$kmeans_neighborhood
cd8t_me_counts_tab[, c("kmeans_neighborhood", "Immune_other", "Unassigned") := NULL]
mat_dat = data.frame(cd8t_me_counts_tab)
rownames(mat_dat) = cluster_names
mat_dat = scale(mat_dat) #z-score columns
mat_dat = pmin(mat_dat, 3) #cap at 3
# Make colors
mat_colors = cd8t_me_colors$color
names(mat_colors) = cd8t_me_colors$kmeans_neighborhood
mat_colors = list(pheno = mat_colors)
# Annotation bar
mat_col = data.frame(pheno = cluster_names)
rownames(mat_col) = cluster_names
# Make heatmap
range = 2
breaks = seq(-range,range,length.out=100)

pdf("cd8t_kmeans_heatmap_input_features.pdf", height=8, width=10)
pheatmap(mat_dat,
         color = rwb_cols,
         cluster_rows = FALSE,
         breaks = breaks,
         annotation_row = mat_col,
         annotation_colors = mat_colors,
         main = "Z-score by column")
dev.off()

## Functional marker positivity
keep_markers = c("Biotin","Caspase1","CD45RO","CD69","CXCR5","Galectin9","Glut1","GranzymeB","ICOS","IFNg","Ki67","Lag3","NLRP3","PD1","TCF1TCF7","TIGIT","TIM3","Vimentin")
functional_tab_with_nh = cd8t_me_cell_tab[,c("fov","label","kmeans_neighborhood")][functional_tab[cell_meta_cluster=="CD8T"], on=.(fov,label)]
functional_tab_with_nh = functional_tab_with_nh[metadata, on=.(fov)]
functional_tab_with_nh = functional_tab_with_nh[!is.na(kmeans_neighborhood)]
count_tab = functional_tab_with_nh[, lapply(.SD, sum, na.rm = TRUE), by = .(kmeans_neighborhood), .SDcols=keep_markers]
total_cd8t = functional_tab_with_nh[, .(total_cd8t = .N), by=.(kmeans_neighborhood)]
count_tab = count_tab[total_cd8t, on=.(kmeans_neighborhood)]
count_tab[, (keep_markers) := lapply(.SD, function(x) x/total_cd8t), .SDcols = keep_markers]
count_tab = count_tab[order(kmeans_neighborhood)]
mat_dat = data.frame(count_tab[,..keep_markers])
rownames(mat_dat) = count_tab$kmeans_neighborhood
mat_dat = scale(mat_dat) #z-score columns
mat_dat = pmin(mat_dat, 3) #cap at 3
range = 2
breaks = seq(-range,range,length.out=100)
pdf("cd8t_kmeans_heatmap_functional_percent.pdf", height=8, width=10)
pheatmap(mat_dat,
         color = rwb_cols,
         breaks = breaks,
         cluster_rows = FALSE)
dev.off()



### Supp Fig 6h: CD8T ME inside vs outside follicle
functional_tab_with_features = cell_tab[,c("fov","label","in_follicle_mask")][functional_tab, on=.(fov,label)]
functional_tab_with_features = cd8t_me_cell_tab[,c("fov","label","kmeans_neighborhood")][functional_tab_with_features[cell_meta_cluster=="CD8T"], on=.(fov,label)]
functional_tab_with_features = functional_tab_with_features[!is.na(kmeans_neighborhood)]
in_follicle_mask = functional_tab_with_features[,.N,by=.(kmeans_neighborhood, in_follicle_mask)]
ggplot(in_follicle_mask, aes(x = kmeans_neighborhood, y=N, fill=in_follicle_mask)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "#B07AA1"),
                    labels = c("FALSE" = "Outside follicle", "TRUE" = "Inside follicle")) +
  theme_bw() +
  labs(y="Number of cells") +
  theme(axis.title.x = element_blank())
ggsave("cd8t_kmeans_follicle_breakdown.pdf", height=3, width=5)



### Supp Fig 6i: CD8T ME distance from follicle edge
follicle_dist_tab = cd8t_follicle_dist_tab[cd8t_me_cell_tab[,c("fov","label","cell_meta_cluster","kmeans_neighborhood")], on=.(fov, label, cell_meta_cluster)]
follicle_dist_tab = follicle_dist_tab[metadata, on=.(fov)]
follicle_dist_tab = cell_tab[,c("fov","label","in_follicle_mask")][follicle_dist_tab,, on=.(fov,label)]
follicle_dist_tab = follicle_dist_tab[!is.na(kmeans_neighborhood)]
follicle_dist_tab$kmeans_neighborhood = factor(follicle_dist_tab$kmeans_neighborhood, levels=c("1","2","3","4"))
# Remove FOVs that don't have any follicles
follicle_dist_tab = follicle_dist_tab[!is.na(dist_to_gc)]
ggplot(follicle_dist_tab[in_follicle_mask==FALSE], aes(x=kmeans_neighborhood, y=dist_to_gc*(800/2048))) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width", fill = "gray") +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y="Distance to follicle edge (um)")
ggsave("cd8t_kmeans_dist_to_follicle.pdf", height=3, width=3)



### Supp Fig 6j: CD8T ME distribution between p24 groups
nh_counts = cd8t_me_cell_tab[,.N,by=.(kmeans_neighborhood, fov)]
nh_counts = nh_counts[!is.na(kmeans_neighborhood)]
nh_counts_wide = dcast(nh_counts, fov ~ kmeans_neighborhood, value.var = "N", fill = 0)
nh_counts = melt(nh_counts_wide, id.vars="fov", value.name = "N", variable.name = "kmeans_neighborhood")

all_cd8t_counts = cd8t_me_cell_tab[, .(total_cd8t=.N), by=.(fov)]
nh_counts = nh_counts[all_cd8t_counts, on=.(fov)]
nh_counts[, prop := N/total_cd8t]
nh_counts = nh_counts[metadata, on=.(fov)]

ggplot(nh_counts, aes(x=status_with_p24, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.2) +
  facet_wrap(~kmeans_neighborhood, scales = "free_y") +
  labs(y="Proportion of CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        strip.background = element_blank())
ggsave("cd8t_me_proportions.pdf", height=4.5, width=4.5)



### Supp Fig 6c, d, k: Look at combined spatial metrics for CD8T
# Get order of niches
quiche_sig_niches_names = quiche_sig_niches[logFC > 0][order(logFC)]$quiche_niche_neighborhood
mean_quiche_cell_tab = quiche_cell_tab[, lapply(.SD, mean, na.rm=TRUE), by=.(quiche_niche_neighborhood), .SDcols=c("logFC")]
quiche_order = mean_quiche_cell_tab[!is.na(logFC)][order(-logFC)]$quiche_niche_neighborhood
quiche_order = quiche_order[quiche_order %in% quiche_sig_niches_names]

# Combine CD8T spatial metrics into one table
new_tab = cd8t_me_cell_tab[,c("fov","label","cell_meta_cluster","kmeans_neighborhood")]
new_tab = cell_tab[cell_meta_cluster=="CD8T",c("fov","label","cell_meta_cluster","in_follicle_mask")][new_tab, on=.(fov,label,cell_meta_cluster)]
new_tab = quiche_cell_tab[cell_meta_cluster=="CD8T",c("fov","label","cell_meta_cluster","quiche_niche_neighborhood")][new_tab, on=.(fov,label,cell_meta_cluster)]
new_tab$kmeans_neighborhood = factor(new_tab$kmeans_neighborhood)

# Percentage of cells in each quiche niche inside/outside
keep_tab = new_tab[quiche_niche_neighborhood %in% quiche_order]
keep_tab$quiche_niche_neighborhood = factor(keep_tab$quiche_niche_neighborhood, levels=quiche_order)
ggplot(keep_tab, aes(x=quiche_niche_neighborhood, fill=in_follicle_mask)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "#B07AA1"),
                    labels = c("FALSE" = "Outside follicle", "TRUE" = "Inside follicle")) +
  labs(x = "QUICHE niche neighborhood", y = "Proportion of cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("quiche_niche_inside_outside.pdf", height=4, width=8)

# Percentage of cells in each quiche niche inside/outside - split by p24 status
keep_tab = metadata[keep_tab, on=.(fov)]
ggplot(keep_tab, aes(x=quiche_niche_neighborhood, fill=in_follicle_mask)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "#B07AA1"),
                    labels = c("FALSE" = "Outside follicle", "TRUE" = "Inside follicle")) +
  facet_wrap(~status_with_p24) +
  labs(x = "QUICHE niche neighborhood", y = "Proportion of cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Which niches are inside the follicles
inside_follicles = cell_tab[,c("fov","label","cell_meta_cluster","in_follicle_mask")][quiche_cell_tab, on=.(fov,label,cell_meta_cluster)]
inside_follicles = inside_follicles[in_follicle_mask==TRUE & cell_meta_cluster=="CD8T"]
top = inside_follicles[, .N, by = quiche_niche_neighborhood][order(-N)][1:10]
top_names = top$quiche_niche_neighborhood
ggplot(inside_follicles[quiche_niche_neighborhood %in% top_names], aes(x = reorder(quiche_niche_neighborhood, quiche_niche_neighborhood, function(x) -length(x)))) +
  geom_bar() +
  labs(x = "QUICHE niche neighborhood", y = "Number of cells", title = "QUICHEs inside the follicles") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
ggsave("quiche_niche_inside_breakdown.pdf", height=4, width=5)


# Which niches are in which kmeans neighborhood
me_colors = setNames(cd8t_me_colors$color, cd8t_me_colors$phenotype_num)
ggplot(keep_tab, aes(x=quiche_niche_neighborhood, fill=kmeans_neighborhood)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = me_colors) +
  labs(x = "QUICHE niche neighborhood", y = "Proportion of cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("quiche_niche_kmeans.pdf", height=4, width=8)
