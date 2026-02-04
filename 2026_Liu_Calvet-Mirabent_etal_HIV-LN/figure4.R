# Code for Figure 4, Supplementary Figure 8
# Author: Candace Liu

library(data.table)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(nlme)
library(ggpubr)

follicle_cell_tab = fread("../data/tables/cell_table_follicles.csv")
follicle_regions = fread("../data/tables/follicle_regions.csv")
functional_tab = fread("../data/tables/nimbus_binarized.csv")
exp_cell_tab = fread("../data/tables/cell_table_size_normalized.csv")
metadata = fread("../data/tables/metadata.csv")
metadata[,sample_id := as.character(sample_id)]

all_markers = colnames(follicle_regions)
all_markers = all_markers[!all_markers %in% c("fov","follicle_label","follicle_size","follicle_p24")]
follicle_regions = metadata[follicle_regions, on=.(fov)]

color_mapping = c("3"="#A3BE8C", "4"="#B48EAD")

# Function for getting pvalues
pval_p24_follicles <- function(dat, feature, column_name) {
  formula = as.formula(paste(feature, paste0("~",column_name)))
  model = lme(formula, random= ~1|sample_id, data = dat, na.action=na.omit)
  p_value = summary(model)$tTable[paste0(column_name,'TRUE'),'p-value']
  return(p_value)
}



### Fig 4b: Cell density differences between follicles
breakdown = follicle_cell_tab[,.N,by=.(fov, cell_meta_cluster, follicle_mask_label, follicle_mask_size, follicle_mask_p24)]
breakdown[,density := N/follicle_mask_size/0.1526*1e6]
breakdown = metadata[breakdown, on=.(fov)]
median_data = breakdown[grepl("slide3|slide4", fov)] %>%
  group_by(sample_id, follicle_mask_p24, cell_meta_cluster) %>%
  summarize(median_value = median(density, na.rm = TRUE), .groups = "drop")

one_cell_pheno = "CD4T"
sub_breakdown = breakdown[cell_meta_cluster==one_cell_pheno]
sub_breakdown = sub_breakdown[grepl("slide3|slide4", fov)]
ggplot(sub_breakdown, aes(x=follicle_mask_p24, y=density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, aes(color=sample_id)) +
  geom_line(data = filter(median_data, cell_meta_cluster == one_cell_pheno),
            aes(y = median_value, group = sample_id, color = sample_id),
            linewidth = 1) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(y = "Cell density (cells/mm2)",
       title = paste0(one_cell_pheno, " density")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5,),
        legend.position = "none")
ggsave(paste0("p24follicle_",one_cell_pheno,"_density.pdf"), height=3.5, width=2.75)

# Add p-value
feature = "density"
pval = pval_p24_follicles(sub_breakdown, feature, 'follicle_mask_p24')
ggplot(sub_breakdown, aes(x=follicle_mask_p24, y=density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, aes(color=sample_id)) +
  geom_line(data = filter(median_data, cell_meta_cluster == one_cell_pheno),
            aes(y = median_value, group = sample_id, color = sample_id),
            linewidth = 1) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(y = "Cell density (cells/mm2)",
       title = paste0(one_cell_pheno, " density")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5,),
        legend.position = "none") +
  stat_pvalue_manual(
    data.frame(
      group1 = "FALSE",
      group2 = "TRUE", 
      p = round(pval,2),
      y.position = max(sub_breakdown[[feature]]) * 1.1
    ),
    label = "p = {p}")
ggsave(paste0("p24follicle_",one_cell_pheno,"_density_withpval.pdf"), height=3.5, width=2.75)



### Fig 4c: Mean expression in follicles
follicle_regions_melt = melt(follicle_regions, id.vars=c("fov","follicle_label","follicle_size","follicle_p24","status","sample_id"), measure.vars=all_markers)
median_data = follicle_regions_melt[grepl("slide3|slide4", fov)] %>%
  group_by(sample_id, follicle_p24, variable) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = "drop")

one_feature = "CD86"
sub_dat = follicle_regions_melt[variable==one_feature]
sub_dat = sub_dat[grepl("slide3|slide4", fov)]
ggplot(sub_dat, aes(x=follicle_p24, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, aes(color=sample_id)) +
  geom_line(data = filter(median_data, variable == one_feature),
            aes(y = median_value, group = sample_id, color = sample_id),
            linewidth = 1) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(y = "Mean expression",
       title = paste0(one_feature, " expression")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5,),
        legend.position = "none")
ggsave(paste0("p24follicle_",one_feature,"_mean_exp.pdf"), height=3.5, width=2.75)

# Add p-value
feature = "value"
pval = pval_p24_follicles(sub_dat, feature, 'follicle_p24')
ggplot(sub_dat, aes(x=follicle_p24, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, aes(color=sample_id)) +
  geom_line(data = filter(median_data, variable == one_feature),
            aes(y = median_value, group = sample_id, color = sample_id),
            linewidth = 1) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(y = "Mean expression",
       title = paste0(one_feature, " expression")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_pvalue_manual(
    data.frame(
      group1 = "FALSE",
      group2 = "TRUE", 
      p = round(pval,2),
      y.position = max(sub_dat[[feature]]) * 1.1
    ),
    label = "p = {p}")
ggsave(paste0("p24follicle_",one_feature,"_mean_exp_withpval.pdf"), height=3.5, width=2.75)



### Supp Fig 8a,c: Breakdown of Ki67+ and CD86+ cells within follicles
keep_follicle_tab = follicle_cell_tab[grepl("slide3|slide4", fov)]
combined = functional_tab[metadata, on=.(fov)]
combined = combined[keep_follicle_tab, on=.(fov,label,cell_meta_cluster)]

one_marker = "CD86"
count_tab = combined[get(one_marker)==1, .N, by=.(fov, cell_meta_cluster, follicle_mask_p24, follicle_mask_size, sample_id)]
count_tab[, density := N/follicle_mask_size/0.1526*10e6]

# Only keep top phenos and remove infrequent points
keep_phenos = count_tab %>%
  group_by(cell_meta_cluster) %>%
  summarize(median_density = median(density, na.rm = TRUE), n = n()) %>%
  filter(n >= 3) %>%
  arrange(desc(median_density)) %>%
  slice_head(n =5) %>%
  pull(cell_meta_cluster)

count_tab = count_tab %>%
  filter(cell_meta_cluster %in% keep_phenos) %>%
  mutate(cell_meta_cluster = factor(cell_meta_cluster, 
                                    levels = keep_phenos))

ggplot(count_tab, aes(x = cell_meta_cluster, y = density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7, aes(color=sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_mapping) +
  labs(y = "Cell density (cells/mm2)",
       title = paste0(one_marker,"+ cells")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") 
ggsave(paste0("p24follicle_", one_marker, "_breakdown_bar.pdf"), height=3, width=3)



### Supp Fig 8b, d: Expression of Ki67 and CD86 in certain cell types
keep_follicle_tab = follicle_cell_tab[grepl("slide3|slide4", fov)]
combined = exp_cell_tab[keep_follicle_tab, on=.(fov,label,cell_meta_cluster)]

mean_tab = combined[,lapply(.SD,mean), .SDcols = all_markers, by = .(fov,cell_meta_cluster, follicle_mask_p24)]
mean_tab = metadata[mean_tab, on=.(fov)]

one_marker = "Ki67"
one_pheno = "CD11c_CD68"

median_data = mean_tab %>%
  group_by(sample_id, follicle_mask_p24, cell_meta_cluster) %>%
  summarize(median_value = median(get(one_marker), na.rm = TRUE), .groups = "drop")

subset_mean_tab = mean_tab[cell_meta_cluster==one_pheno]
subset_median_tab = median_data %>%
  filter(cell_meta_cluster == one_pheno)

ggplot(subset_mean_tab, aes(x=follicle_mask_p24, y=get(one_marker))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7, aes(color=sample_id)) +
  geom_line(data = subset_median_tab,
            aes(y = median_value, group = sample_id, color = sample_id),
            linewidth = 0.7) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(y = "Mean expression",
       title = paste0(one_marker,"+ ", one_pheno)) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(paste0("p24follicle_", one_marker, "_", one_pheno,"_exp.pdf"), height=3, width=2.75)



### Fig 4d-e, Supp Fig 8f: Distance analysis - cell density around p24+ cells
dist_count = fread("../data/p24pos_cells_dist_counts_0_500.csv")
unique_fov_index = unique(paste0(dist_count$fov,"-",dist_count$index_cell_label))
unique_cell_meta_cluster = unique(dist_count$cell_meta_cluster)
unique_window = unique(dist_count$window)
all_combinations = CJ(
  fov_index = unique_fov_index,
  cell_meta_cluster = unique_cell_meta_cluster,
  window = unique_window
)
all_combinations[, c("fov","index_cell_label") := tstrsplit(fov_index, "-", fixed = TRUE)]
all_combinations[, index_cell_label := as.numeric(index_cell_label)]

dist_count = merge(
  all_combinations[,c("fov","index_cell_label","cell_meta_cluster","window")],
  dist_count,
  by = c("fov", "index_cell_label", "cell_meta_cluster", "window"),
  all.x = TRUE
)
dist_count[is.na(count), count := 0]

# Get area of concentric circles for normalization
window_areas = data.table(window=unique(dist_count$window))
window_size = unique(diff(window_areas$window))
window_areas[, window_area := (pi * window ^2) - (pi * (window-window_size)^2)]
dist_count = dist_count[window_areas, on=.(window)]
dist_count[,count_norm := count / window_area]

mean_by_fov = dist_count[, .(mean_count = mean(count_norm)), by = .(fov, cell_meta_cluster, window)]
mean_by_fov = mean_by_fov[, slide := gsub("^slide(.*?)_.*$", "\\1", fov)]
summary_stats = mean_by_fov[, .(
  mean_across_fovs = mean(mean_count),
  se = sd(mean_count) / sqrt(.N),  # standard error
  lower_ci = mean(mean_count) - 1.96 * sd(mean_count) / sqrt(.N),  # 95% CI lower bound
  upper_ci = mean(mean_count) + 1.96 * sd(mean_count) / sqrt(.N)   # 95% CI upper bound
), by = .(cell_meta_cluster, window, slide)]

one_pheno = "CD11c"
ggplot(summary_stats[cell_meta_cluster==one_pheno], aes(x = window, y = mean_across_fovs, group = slide, color=slide, fill=slide)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Shaded confidence interval
  geom_line(linewidth = 1) +  # Mean line
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  theme_bw() +
  labs(x = "Window",
       y = "Mean normalized count",
       title = one_pheno) +
  theme(plot.title = element_text(hjust = 0.5,),
        legend.position = "none")
ggsave(paste0("p24cell_distance_",one_pheno,".pdf"), height=3, width=3)

keep_features = c("B","CD11c_CD68","CD11c_CD14","CD14_CD68_CD163","CD4T","Foxp3","NK","SMA","Endothelial")
subset_summary_stats = summary_stats[cell_meta_cluster %in% keep_features]
subset_summary_stats$cell_meta_cluster = factor(subset_summary_stats$cell_meta_cluster, levels = keep_features)
ggplot(subset_summary_stats, aes(x = window, y = mean_across_fovs, group = slide, color=slide, fill=slide)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Shaded confidence interval
  geom_line(linewidth = 1) +  # Mean line
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  facet_wrap(~cell_meta_cluster, scales = "free") +
  theme_bw() +
  labs(x = "Window",
       y = "Mean normalized count") +
  theme(legend.position = "none",
        strip.background = element_blank())
ggsave("p24cell_distance_selected.pdf", height=6, width=7)



### Supp Fig 8g: Functional markers in CD8T cells
p24_cd8_tab = fread("../data/p24pos_cells_dist_counts_0_500_cd8t_func.csv")

unique_fov_index = unique(paste0(p24_cd8_tab$fov,"-",p24_cd8_tab$index_cell_label))
unique_marker = unique(p24_cd8_tab$marker)
unique_window = unique(p24_cd8_tab$window)
all_combinations = CJ(
  fov_index = unique_fov_index,
  marker = unique_marker,
  window = unique_window
)
all_combinations[, c("fov","index_cell_label") := tstrsplit(fov_index, "-", fixed = TRUE)]
all_combinations[, index_cell_label := as.numeric(index_cell_label)]
p24_cd8_tab = merge(
  all_combinations[,c("fov","index_cell_label","marker","window")],
  p24_cd8_tab,
  by = c("fov", "index_cell_label", "marker", "window"),
  all.x = TRUE
)
p24_cd8_tab[is.na(prop), prop := 0]

mean_by_fov = p24_cd8_tab[, .(mean_prop = mean(prop)), by = .(fov, marker, window)]
mean_by_fov = mean_by_fov[, slide := gsub("^slide(.*?)_.*$", "\\1", fov)]
summary_stats = mean_by_fov[, .(
  mean_across_fovs = mean(mean_prop),
  se = sd(mean_prop) / sqrt(.N),  # standard error
  lower_ci = mean(mean_prop) - 1.96 * sd(mean_prop) / sqrt(.N),  # 95% CI lower bound
  upper_ci = mean(mean_prop) + 1.96 * sd(mean_prop) / sqrt(.N)   # 95% CI upper bound
), by = .(marker, window, slide)]

subset_markers = c("CXCR5","Caspase1","CD45RO","Glut1","GranzymeB","HLADR","Ki67","PD1","TIGIT")
ggplot(summary_stats[marker %in% subset_markers], aes(x = window, y = mean_across_fovs, group = slide, color=slide, fill=slide)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  # Shaded confidence interval
  geom_line(linewidth = 1) +  # Mean line
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  facet_wrap(~marker, scales = "free") +
  theme_bw() +
  labs(x="Window",
       y="Proportion of CD8+ T cells in window") +
  theme(strip.background = element_blank(),
        legend.position = "none")
ggsave("p24cell_distance_cd8t_func_selected.pdf", height=6, width=7) 



### Fig 4f: Enrichment score
enrichment_score_tab = fread("../data/spatial_analysis/p24_cell_enrichment_20um/p24pos_enrichment_scores.csv")
enrichment_score_tab = metadata[enrichment_score_tab, on=.(fov)]
# Context-dependent for FDCs
enrichment_score_tab_context = fread("../data/spatial_analysis/p24_cell_enrichment_20um_context_dependent_follicle/p24pos_enrichment_scores.csv")
enrichment_score_tab_context = metadata[enrichment_score_tab_context, on=.(fov)]
enrichment_score_tab_context[, pheno := paste0(pheno,"_context")]
enrichment_score_tab_context = rbind(enrichment_score_tab, enrichment_score_tab_context)

select_phenos = c('FDC_context','CD11c_CD68','Tfh','B','CD4T','CD8T')
select_dt = enrichment_score_tab_context[pheno %in% select_phenos]
select_dt[, pheno := factor(pheno, levels=select_phenos)]
select_dt = select_dt[order(pheno)]
ggplot(select_dt, aes(x=pheno, y=z)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, aes(color=sample_id)) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(x = "Cell phenotype",
       y = "Enrichment score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave("../data/plots/paper/p24_enrichment_selected_context.pdf", height=3, width=4)



### Supp Fig 8h: Enrichment score between Caspase1+/- CD8T and p24+ cells
file_path =  "../data/spatial_analysis/p24_cell_enrichment_20um/p24pos_CD8T_Caspase1_enrichment_scores.csv"
one_tab = fread(file_path)
one_tab = metadata[one_tab, on=.(fov)]

# Add p-value
model = lme(z ~ pheno, random= ~1|sample_id, data = one_tab, na.action=na.omit)
p_value = summary(model)$tTable['phenoCD8T_1','p-value']
ggplot(one_tab, aes(x=pheno, y=z)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, aes(color=sample_id)) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(y = "Enrichment score",
       title = "Enrichment around p24+ cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5,),
        legend.position = "none") +
  stat_pvalue_manual(
    data.frame(
      group1 = "CD8T_0",
      group2 = "CD8T_1", 
      p = round(p_value,2),
      y.position = max(one_tab$z) * 1.2
    ),
    label = "p = {p}")
ggsave("caspase1_cd8t_p24pos_withpval.pdf", height=3.5, width=3)
