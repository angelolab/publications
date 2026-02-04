# Code for Figure 2, Supplementary Figure 5, Supplementary Figure 7
# Author: Candace Liu

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(UpSetR)
library(nlme)
library(emmeans)

cell_tab = fread("../data/tables/cell_table_size_normalized.csv")
cell_tab_follicle = fread("../data/tables/cell_table_follicles.csv")
functional_tab = fread("../data/tables/nimbus_binarized.csv")
feature_tab = fread("../data/tables/feature_tab.csv")
mixed_effect_model_tab = fread("../data/tables/mixed_effect_model_output.csv")
tissue_area = fread("../data/tables/tissue_area.csv")
metadata = fread("../data/tables/metadata.csv")


# Function for getting pvalues between p24 groups
pval_p24_groups <- function(feature) {
  formula = as.formula(paste(feature,"~status_with_p24"))
  model = lme(formula, random= ~1|sample_id, data = feature_tab, na.action=na.omit)
  emm = emmeans(model, "status_with_p24")
  
  pairwise_results = pairs(emm, adjust = "tukey")
  pairwise_df = as.data.frame(pairwise_results)
  return(pairwise_df)
}


### Fig 2a: Mixed effect model volcano plot
all_features = colnames(feature_tab)
all_features = setdiff(all_features, c("fov","status","sample_id","status_with_viremia","status_with_viremia","status_with_p24"))
feature_tab_withmeta = metadata[,c("fov","sample_id")][feature_tab, on=c("fov","sample_id")]
feature_tab_withmeta[,fov:=NULL]
long_dt = melt(feature_tab_withmeta,
               id.vars = "status",
               measure.vars = all_features,
               variable.name = "feature_name", 
               value.name = "value")
stats_dt = long_dt[, .(
  n = sum(!is.na(value)),
  mean = mean(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE)
), by = .(feature_name, status)]

stats_wide = dcast(stats_dt, feature_name ~ status, 
                   value.var = c("n", "mean", "sd"))
results_dt_es = stats_wide[, .(
  feature_name = feature_name,
  cohens_d = (mean_hiv_pos - mean_hiv_neg) / 
    sqrt(((n_hiv_pos - 1) * sd_hiv_pos^2 + (n_hiv_neg - 1) * sd_hiv_neg^2) / 
           (n_hiv_pos + n_hiv_neg - 2))
)]
results_dt_es[, feature_name := trimws(feature_name)]
results_dt_es = results_dt_es[mixed_effect_model_tab, on=.(feature_name)]
results_dt_es$group = ifelse(results_dt_es$p_value < 0.05 & results_dt_es$cohens_d > 0.5, "hivpos_sig",
                             ifelse(results_dt_es$p_value < 0.05 & results_dt_es$cohens_d < -0.5, "hivneg_sig",
                                    "other"))
results_dt_es[, to_label := FALSE]
results_dt_es[feature_name=="CD4T_density", to_label := TRUE]
results_dt_es[feature_name=="CD8T_density", to_label := TRUE]

ggplot(results_dt_es, aes(x = cohens_d, y = -log10(p_value))) +
  geom_point(aes(color = group), shape=19) +
  scale_color_manual(values = c("hivpos_sig" = "#BF616A", "hivneg_sig" = "#5E81AC", "other" = "darkgray")) +
  theme_bw() +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
  labs(x = "Standardized effect size (Cohen's d)", y = "-log10(p-value)") +
  geom_text_repel(
    data = results_dt_es[to_label == TRUE],
    aes(label = feature_name),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.5,
    segment.alpha = 0.8,
    min.segment.length = 0,
    max.overlaps = 20
  ) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave("volcano_plot_cohensd.pdf", height=5, width=6)



### Fig 2b: Boxplots of CD4/CD8 ratio
keep_tab = feature_tab[,c("fov","sample_id","status","CD4T_CD8T_ratio")]
keep_tab[, median_ratio := median(log2(CD4T_CD8T_ratio)), by = sample_id]
keep_tab[, sample_id := factor(sample_id, levels = unique(sample_id[order(median_ratio)]))]
ggplot(keep_tab, aes(x=log2(CD4T_CD8T_ratio), y=sample_id, color=status)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(height=0.05) +
  geom_vline(xintercept=0, linetype="dashed", color='gray', linewidth=1) +
  scale_color_manual(values = c("hiv_pos" = "#BF616A", "hiv_neg" = "#5E81AC")) +
  theme_bw() +
  labs(x = "log2(# CD4+ T cells / # CD8+ T cells)", y = "Donor") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
ggsave("cd4_cd8_ratio.pdf", height=8, width=7)



### Fig 2c, Supp Fig 5b-f, Supp Fig 7c-d, h: Boxplot of features compared between p24 groups
feature = "CD8T_density"
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = "CD8+ T cell density") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithp24.pdf"), height=4, width=3.75)

# Add p-value (density)
feature = "CD4T_density"
pairwise_df_out = pval_p24_groups(feature)
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
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Density (cells/mm2)",
       title = "Caspase1+ CD8+ T cells") +
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
ggsave(paste0(feature, "_statuswithp24_withpval.pdf"), height=4, width=3.5)

ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Density (cells/mm2)",
       title = "Caspase1+ CD8+ T cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithp24.pdf"), height=4, width=3.5)



# Add p-value (proportion)
feature = "Caspase1pos_CD8T_prop"
#feature = "GranzymeBpos_CD8T_prop"
pairwise_df_out = pval_p24_groups(feature)
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
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of CD8+ T cells",
       title = "Caspase-1+ CD8+ T cells") +
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
ggsave(paste0(feature, "_statuswithp24_withpval.pdf"), height=4, width=3.75)


feature = "Ki67pos_CD8T_density"
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Density (cells/mm2)",
       title = "Ki67+ CD8+ T cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithp24.pdf"), height=3, width=2.75)



### Fig 2d: Breakdown of number of functional markers expressed in CD8T
functional_markers = c("Biotin","Caspase1","CD45RO","CD69","Galectin9","Glut1","GranzymeB","ICOS","IFNg","Ki67","Lag3","NLRP3","PD1","TCF1TCF7","TIGIT","TIM3")
keep_tab = functional_tab[cell_meta_cluster=="CD8T"]
keep_tab[, num_markers := rowSums(.SD), .SDcols = functional_markers]
keep_tab$num_markers_grouped = factor(
  ifelse(keep_tab$num_markers > 3, ">3", as.character(keep_tab$num_markers)),
  levels = c("0", "1", "2", "3", ">3")
)
keep_tab_count = keep_tab[,.N,by=.(fov,num_markers_grouped)]
keep_tab_count = keep_tab_count[CJ(fov = unique(fov), num_markers_grouped = unique(num_markers_grouped)), on = .(fov, num_markers_grouped)] #fill with 0's for missing
keep_tab_count[is.na(N), N := 0]
keep_tab_count = keep_tab_count[metadata, on=.(fov)]
total_cd8t_counts = keep_tab[,.(total_cd8t = .N), by=.(fov)]
keep_tab_prop = keep_tab_count[total_cd8t_counts, on=.(fov)]
keep_tab_prop[, prop := N/total_cd8t]
keep_tab_mean = keep_tab_prop[,lapply(.SD, mean), .SDcols=c("prop"), by=.(status_with_p24, num_markers_grouped)]
ggplot(keep_tab_mean, aes(x = status_with_p24, y = prop, fill = num_markers_grouped)) +
  geom_col() +
  scale_fill_viridis_d() +
  labs(y="Mean proportion of CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank())
ggsave("naive_cd8t_breakdown_statuswithp24.pdf", height=3, width=4.5)



### Fig 2e-f: Boxplots of functional marker+ CD8T (proportion of total CD8T)
feature = "Caspase1pos_CD8T_prop"
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of CD8+ T cells",
       title = "Caspase1+ CD8+ T cells") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithp24.pdf"), height=4, width=3.75)



### Supp Fig 5g: Boxplots of double-positive CD8T (proportion of total CD8T)
keep_markers = c("Caspase1","GranzymeB")
sub_cell_tab = functional_tab[cell_meta_cluster=="CD8T"]
sub_cell_tab_keep = sub_cell_tab[rowSums(sub_cell_tab[, ..keep_markers] == 1) == length(keep_markers)]
sub_cell_counts = sub_cell_tab_keep[, .N, by=.(fov)]
totals = sub_cell_tab[,.(totals = .N), by=.(fov)]
sub_cell_counts = sub_cell_counts[totals, on=.(fov)]
sub_cell_counts[is.na(sub_cell_counts)] = 0
sub_cell_counts[, prop := N/totals]
sub_cell_counts = metadata[sub_cell_counts, on=.(fov)]

ggplot(sub_cell_counts, aes(x=status_with_p24, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.7) +
  labs(y = "Proportion of CD8+ T cells",
       title = paste0(paste(keep_markers, collapse = "+"),"+ CD8+ T cells")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0("cd8t_functional_proportion_",paste(keep_markers, collapse = "_"),".pdf"), height=4, width=3.5)

# Add p-value (proportion)
formula = as.formula(paste("prop","~status_with_p24"))
model = lme(formula, random= ~1|sample_id, data = sub_cell_counts, na.action=na.omit)
emm = emmeans(model, "status_with_p24")
pairwise_results = pairs(emm, adjust = "tukey")
pairwise_df_out = as.data.frame(pairwise_results)

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
ggplot(sub_cell_counts, aes(x=status_with_p24, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Proportion of CD8+ T cells",
       title = paste0(paste(keep_markers, collapse = "+"),"+ CD8+ T cells")) +
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
    y.position = c(max(sub_cell_counts$prop, na.rm = TRUE) * c(1.1, 1.3, 1.2))
  )
ggsave(paste0("cd8t_functional_proportion_",paste(keep_markers, collapse = "_"),"_withpval.pdf"), height=4, width=3.5)



### Supp Fig 5h: Boxplots of double-positive CD8T (density)
keep_markers = c("Caspase1","GranzymeB")
sub_cell_tab = functional_tab[cell_meta_cluster=="CD8T"]
sub_cell_tab_keep = sub_cell_tab[rowSums(sub_cell_tab[, ..keep_markers] == 1) == length(keep_markers)]
sub_cell_counts = sub_cell_tab_keep[, .N, by=.(fov)]
sub_cell_counts = sub_cell_counts[tissue_area, on=.(fov)]
sub_cell_counts[is.na(sub_cell_counts)] = 0
sub_cell_counts[, density := N/tissue_area]
sub_cell_counts = metadata[sub_cell_counts, on=.(fov)]

ggplot(sub_cell_counts, aes(x=status_with_p24, y=density/0.64)) + #0.64 because tissue area was calculated as percentage of image, image is 800umx800um
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.7) +
  labs(y = "Density (cells/mm2)",
       title = paste0(paste(keep_markers, collapse = "+"),"+ CD8+ T cells")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0("cd8t_functional_density_",paste(keep_markers, collapse = "_"),".pdf"), height=4, width=3.5)



### Supp Fig 5j: GranzymeB+Caspase1- cells
sub_cell_tab = functional_tab[cell_meta_cluster=="CD8T"]
sub_cell_tab_keep = sub_cell_tab[GranzymeB==1 & Caspase1==0]
sub_cell_counts = sub_cell_tab_keep[, .N, by=.(fov)]
totals = sub_cell_tab[,.(totals = .N), by=.(fov)]
sub_cell_counts = sub_cell_counts[totals, on=.(fov)]
sub_cell_counts[is.na(sub_cell_counts)] = 0
sub_cell_counts[, prop := N/totals]
sub_cell_counts = metadata[sub_cell_counts, on=.(fov)]

ggplot(sub_cell_counts, aes(x=status_with_p24, y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.7) +
  labs(y = "Proportion of CD8+ T cells",
       title = "GranzymeB+Caspase1- CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("cd8t_granpos_caspneg_proportion.pdf", height=4, width=3.5)



### Supp Fig 5k-l: Plot multiple features
#keep_features = c("Ki67pos_CD8T_prop", "CD45ROpos_CD8T_prop","TIGITpos_CD8T_prop","TIM3pos_CD8T_prop","Galectin9pos_CD8T_prop","PD1pos_CD8T_prop","ICOSpos_CD8T_prop","IFNgpos_CD8T_prop", "CD69pos_CD8T_prop", "TCF1TCF7pos_CD8T_prop")
keep_features = c("Ki67pos_CD8T_density", "CD45ROpos_CD8T_density","TIGITpos_CD8T_density","TIM3pos_CD8T_density","Galectin9pos_CD8T_density","PD1pos_CD8T_density","ICOSpos_CD8T_density","IFNgpos_CD8T_density", "CD69pos_CD8T_density", "TCF1TCF7pos_CD8T_density")
keep_feature_tab_cols = c("fov","status_with_p24",keep_features)
keep_feature_tab_melt = melt(feature_tab[,..keep_feature_tab_cols], id.vars=c("fov","status_with_p24"))
ggplot(keep_feature_tab_melt, aes(x=status_with_p24, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.3) +
  facet_wrap(~variable, scales="free_y", nrow=2) +
  theme_bw() +
  #labs(y = "Proportion of CD8+ T cells") +
  labs(y = "Density (cells/mm2)") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())
ggsave("cd8t_func_density_statuswithp24.pdf", height=4, width=8)



### Supp Fig 5m: Upset plot of double positive CD8T
one_group = "hiv_pos_p24pos"
keep_markers = c("Biotin","Caspase1","CD45RO","CD69","Galectin9","Glut1","GranzymeB","ICOS","IFNg","Ki67","Lag3","NLRP3","PD1","TCF1TCF7","TIGIT","TIM3")
cd8t_func_phenostrings = functional_tab[cell_meta_cluster=="CD8T"][metadata, on=.(fov)]
cd8t_func_phenostrings = cd8t_func_phenostrings[status_with_p24 == one_group]
cd8t_func_phenostrings$pheno_string = do.call(paste0, cd8t_func_phenostrings[,..keep_markers])

# Function for getting the binary string for a phenotype
get_binary_string <- function(phenos) {
  binned_phenos = keep_markers %in% phenos
  binned_phenos = binned_phenos*1
  return_string = paste(binned_phenos, collapse="")
  return(return_string)
}

# Count each combination of phenos
count_table = cd8t_func_phenostrings[,.N, by=.(pheno_string)]
total_cd8t = sum(count_table$N)
count_table[, perc_total_cd8t := N/total_cd8t]

# Only keep degree 2
combs = combn(keep_markers, 2, simplify = FALSE)
# Only keep non-empty groups
set_list = list()
for(marker in keep_markers) {
  set_list[[marker]] = which(cd8t_func_phenostrings[[marker]] == 1)
}
intersection_data = fromList(set_list)
non_empty_combs = list()
for(comb in combs) {
  condition = rep(TRUE, nrow(intersection_data))
  for(marker in keep_markers) {
    if(marker %in% comb) {
      condition = condition & intersection_data[[marker]] == 1
    } else {
      condition = condition & intersection_data[[marker]] == 0
    }
  }
  if(sum(condition) > 100) {
    non_empty_combs = c(non_empty_combs, list(comb))
  }
}

pdf(paste0("cd8t_upset_",one_group,".pdf"), height=4.5, width=6)
upset(cd8t_func_phenostrings, sets = keep_markers, intersections = non_empty_combs)
dev.off()
# Make bar graph as % of total CD8T
new_dt = as.data.table(do.call(rbind, non_empty_combs))
colnames(new_dt) = c("marker1","marker2")
new_dt[, pheno_string := sapply(non_empty_combs, get_binary_string)]

keep_count_table = count_table[pheno_string %in% new_dt$pheno_string]
keep_count_table = keep_count_table[new_dt, on=.(pheno_string)]
keep_count_table = keep_count_table[order(perc_total_cd8t, decreasing=TRUE)]
keep_count_table$marker_combo = paste(keep_count_table$marker1, "+", keep_count_table$marker2)

ggplot(keep_count_table, aes(x = reorder(marker_combo, perc_total_cd8t, decreasing=TRUE), y = perc_total_cd8t)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y="% of CD8+ T cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank())
ggsave(paste0("cd8t_upset_",one_group,"_barplot.pdf"), height=3, width=6)


# Plot all 3 groups together
keep_order_tab = keep_count_table[,c("pheno_string","marker1","marker2","marker_combo")] #from p24+ group
desired_order = keep_count_table$marker_combo
  
cd8t_func_phenostrings = functional_tab[cell_meta_cluster=="CD8T"][metadata, on=.(fov)]
cd8t_func_phenostrings$pheno_string = do.call(paste0, cd8t_func_phenostrings[,..keep_markers])

count_table = cd8t_func_phenostrings[,.N, by=.(pheno_string, status_with_p24)]
total_cd8t = cd8t_func_phenostrings[,.(total_cd8t = .N), by=.(status_with_p24)]
count_table = count_table[total_cd8t, on=.(status_with_p24)]
count_table[, perc_total_cd8t := N/total_cd8t]

merged_count_table = count_table[keep_order_tab, on=.(pheno_string)]
merged_count_table$marker_combo = factor(merged_count_table$marker_combo, levels = desired_order)
ggplot(merged_count_table, aes(x = marker_combo, y = perc_total_cd8t, fill = status_with_p24)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("hiv_neg" = "#5E81AC", "hiv_pos_p24neg" = "#FFBE7D", "hiv_pos_p24pos" = "#BF616A"),
    labels = c("HIV-", "HIV+ p24-", "HIV+ p24+")
  ) +
  labs(
    y = "Percentage of total CD8+ T cells",
    fill = "Group"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("cd8t_upset_allgroups_barplot.pdf", height=3, width=9)



### Supp Fig 5n: Distribution of checkpoint markers in CD8T
checkpoint_markers = c("Biotin","Galectin9","Lag3","PD1","TIGIT","TIM3")
keep_tab = functional_tab[cell_meta_cluster=="CD8T"]
keep_tab[, num_markers := rowSums(.SD), .SDcols=checkpoint_markers]
keep_tab$num_markers_grouped = factor(
  ifelse(keep_tab$num_markers > 1, ">1", as.character(keep_tab$num_markers)),
  levels = c("0", "1", ">1")
)
co_exp_counts = table(keep_tab$num_markers_grouped)
total_num_cd8t = sum(co_exp_counts)
co_exp_counts = co_exp_counts / total_num_cd8t

df = as.data.frame(co_exp_counts)
names(df) = c("group","frequency")
ggplot(df, aes(x=group, y=frequency)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(x="Number of checkpoint markers co-expressed",
       y="Proportion of CD8+ T cells")
ggsave("cd8t_checkpoint.pdf", height=3, width=4)

# Broken down by group
checkpoint_markers = c("Biotin","Galectin9","Lag3","PD1","TIGIT","TIM3")
keep_tab = functional_tab[cell_meta_cluster=="CD8T"]
keep_tab[, num_markers := rowSums(.SD), .SDcols = checkpoint_markers]
keep_tab$num_markers_grouped = factor(
  ifelse(keep_tab$num_markers > 1, ">1", as.character(keep_tab$num_markers)),
  levels = c("0", "1", ">1")
)
keep_tab_count = keep_tab[,.N,by=.(fov,num_markers_grouped)]
keep_tab_count = keep_tab_count[CJ(fov = unique(fov), num_markers_grouped = unique(num_markers_grouped)), on = .(fov, num_markers_grouped)] #fill with 0's for missing
keep_tab_count[is.na(N), N := 0]
keep_tab_count = keep_tab_count[metadata, on=.(fov)]
total_cd8t_counts = keep_tab[,.(total_cd8t = .N), by=.(fov)]
keep_tab_prop = keep_tab_count[total_cd8t_counts, on=.(fov)]
keep_tab_prop[, prop := N/total_cd8t]
keep_tab_mean = keep_tab_prop[,lapply(.SD, mean), .SDcols=c("prop"), by=.(status_with_p24, num_markers_grouped)]
ggplot(keep_tab_mean, aes(x = status_with_p24, y = prop, fill = num_markers_grouped)) +
  geom_col() +
  scale_fill_viridis_d() +
  labs(y="Mean proportion of CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank())
ggsave("cd8t_checkpoint_statuswithp24.pdf", height=2, width=4)



### Supp Fig 5o: Upset plot of checkpoint markers in CD8T
checkpoint_markers = c("Biotin","Galectin9","Lag3","PD1","TIGIT","TIM3")
cd8t_func_phenostrings = functional_tab[cell_meta_cluster=="CD8T"]
cd8t_func_phenostrings$pheno_string = do.call(paste0, cd8t_func_phenostrings[,..checkpoint_markers])

# Function for getting the binary string for a phenotype
get_binary_string <- function(phenos) {
  binned_phenos = checkpoint_markers %in% phenos
  binned_phenos = binned_phenos*1
  return_string = paste(binned_phenos, collapse="")
  return(return_string)
}

# Count each combination of phenos
count_table = cd8t_func_phenostrings[,.N, by=.(pheno_string)]
total_cd8t = sum(count_table$N)
count_table[, perc_total_cd8t := N/total_cd8t]

# Only keep degree 2
combs = combn(checkpoint_markers, 2, simplify = FALSE)
# Only keep non-empty groups
set_list = list()
for(marker in checkpoint_markers) {
  set_list[[marker]] = which(cd8t_func_phenostrings[[marker]] == 1)
}
intersection_data = fromList(set_list)
non_empty_combs = list()
for(comb in combs) {
  condition = rep(TRUE, nrow(intersection_data))
  for(marker in checkpoint_markers) {
    if(marker %in% comb) {
      condition = condition & intersection_data[[marker]] == 1
    } else {
      condition = condition & intersection_data[[marker]] == 0
    }
  }
  if(sum(condition) > 0) {
    non_empty_combs = c(non_empty_combs, list(comb))
  }
}

pdf("cd8t_checkpoint_upset.pdf", height=3.2, width=8)
upset(cd8t_func_phenostrings, sets=checkpoint_markers, intersections = non_empty_combs)
dev.off()
# Make bar graph as % of total CD8T
new_dt = as.data.table(do.call(rbind, non_empty_combs))
colnames(new_dt) = c("marker1","marker2")
new_dt[, pheno_string := sapply(non_empty_combs, get_binary_string)]

keep_count_table = count_table[pheno_string %in% new_dt$pheno_string]
keep_count_table = keep_count_table[new_dt, on=.(pheno_string)]
keep_count_table = keep_count_table[order(perc_total_cd8t, decreasing=TRUE)]
keep_count_table$marker_combo = paste(keep_count_table$marker1, "+", keep_count_table$marker2)

ggplot(keep_count_table, aes(x = reorder(marker_combo, perc_total_cd8t, decreasing=TRUE), y = perc_total_cd8t)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y="% of CD8+ T cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank())
ggsave("cd8t_checkpoint_upset_barplot.pdf", height=2, width=4.5)

# Plot all 3 groups together (for combs > 100 counts)
keep_order_tab = keep_count_table[N > 100, c("pheno_string","marker1","marker2","marker_combo")] #order from all cells
desired_order = keep_count_table$marker_combo

cd8t_func_phenostrings = functional_tab[cell_meta_cluster=="CD8T"][metadata, on=.(fov)]
cd8t_func_phenostrings$pheno_string = do.call(paste0, cd8t_func_phenostrings[,..checkpoint_markers])

count_table = cd8t_func_phenostrings[,.N, by=.(pheno_string, status_with_p24)]
total_cd8t = cd8t_func_phenostrings[,.(total_cd8t = .N), by=.(status_with_p24)]
count_table = count_table[total_cd8t, on=.(status_with_p24)]
count_table[, perc_total_cd8t := N/total_cd8t]

merged_count_table = count_table[keep_order_tab, on=.(pheno_string)]
merged_count_table$marker_combo = factor(merged_count_table$marker_combo, levels = desired_order)
ggplot(merged_count_table, aes(x = marker_combo, y = perc_total_cd8t, fill = status_with_p24)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("hiv_neg" = "#5E81AC", "hiv_pos_p24neg" = "#FFBE7D", "hiv_pos_p24pos" = "#BF616A"),
    labels = c("HIV-", "HIV+ p24- LN", "HIV+ p24+ LN")
  ) +
  labs(
    y = "Percentage of total CD8+ T cells",
    fill = "Group"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("cd8t_checkpoint_upset_allgroups_barplot.pdf", height=3, width=6.5)



### Supp Fig 5p-q, 7e: Co-expression of markers with GranzymeB+/Caspase1+ CD8T cells
keep_markers = c("Biotin","Caspase1","CD45RO","CD69","Galectin9","Glut1","GranzymeB","ICOS","IFNg","Ki67","Lag3","NLRP3","PD1","TCF1TCF7","TIGIT","TIM3")

# Functional marker+ as percentage of GranzymeB+ cells
keep_cells = functional_tab[cell_meta_cluster=="CD8T" & GranzymeB==1]
sub_total = keep_cells[, .N, by=.(fov)]
keep_cells_freq = keep_cells[, lapply(.SD,sum), by=.(fov), .SDcols=keep_markers]
keep_cells_freq = sub_total[keep_cells_freq, on=.(fov)]
keep_cells_freq = metadata[keep_cells_freq, on=.(fov)]
keep_cells_freq[, (keep_markers) := .SD / N, .SDcols = keep_markers]
keep_cells_freq_long = melt(keep_cells_freq, measure.vars = keep_markers)
ggplot(keep_cells_freq_long[variable != "GranzymeB"], aes(x=reorder(variable, value, FUN=median, decreasing=TRUE), y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.3) +
  labs(y="Proportion of GranzymeB+ CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("cd8t_granzymeB_coexpression.pdf", height=3, width=6)

ggplot(keep_cells_freq_long[variable != "GranzymeB"], aes(x=reorder(variable, value, FUN=median, decreasing=TRUE), y=value + 0.0001)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.3) +
  scale_y_continuous(trans = "log2") +
  labs(y="Proportion of GranzymeB+ CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Function marker+ as percentage of Caspase1+ cells
keep_cells = functional_tab[cell_meta_cluster=="CD8T" & Caspase1==1]
sub_total = keep_cells[, .N, by=.(fov)]
keep_cells_freq = keep_cells[, lapply(.SD,sum), by=.(fov), .SDcols=keep_markers]
keep_cells_freq = sub_total[keep_cells_freq, on=.(fov)]
keep_cells_freq = metadata[keep_cells_freq, on=.(fov)]
keep_cells_freq[, (keep_markers) := .SD / N, .SDcols = keep_markers]
keep_cells_freq_long = melt(keep_cells_freq, measure.vars = keep_markers)
ggplot(keep_cells_freq_long[variable != "Caspase1"], aes(x=reorder(variable, value, FUN=median, decreasing=TRUE), y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.3) +
  labs(y="Proportion of Caspase1+ CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("cd8t_caspase1_coexpression.pdf", height=3, width=6)

ggplot(keep_cells_freq_long[variable != "Caspase1"], aes(x=reorder(variable, value, FUN=median, decreasing=TRUE), y=value + 0.0001)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.3) +
  scale_y_continuous(trans = "log2") +
  labs(y="Proportion of Caspase1+ CD8+ T cells") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



### Supp Fig 5r: Look at any checkpoint with Caspase1 OR GranzymeB
new_tab = functional_tab[cell_meta_cluster == "CD8T"]
new_tab = cell_tab[,c("fov","label","in_follicle_mask")][new_tab, on=.(fov,label)]
new_tab[, group := "None"]
new_tab[GranzymeB==1 | Caspase1==1, group := "GranzymeB+ or Caspase1+"]
new_tab[Biotin == 1 | Galectin9 == 1 | Lag3 == 1 | PD1 == 1 | TIGIT == 1 | TIM3 == 1, group := "Checkpoint+"]
new_tab[(GranzymeB==1 | Caspase1==1) & (Biotin == 1 | Galectin9 == 1 | Lag3 == 1 | PD1 == 1 | TIGIT == 1 | TIM3 == 1), group := "Both"]

ggplot(new_tab[group!="None"], aes(x = reorder(group, group, function(x) -length(x)))) +
  geom_bar() +
  theme_bw() +
  theme(axis.title.x = element_blank())
ggsave("cd8t_caspgran_checkpoint_barplot.pdf", height=3, width=3)

ggplot(new_tab[group!="None"], aes(x = reorder(group, group, function(x) -length(x)), fill=in_follicle_mask)) +
  geom_bar(position="fill") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "#B07AA1"),
                    labels = c("FALSE" = "Outside follicle", "TRUE" = "Inside follicle")) +
  theme_bw() +
  labs(y = "Proportion") +
  theme(axis.title.x = element_blank()) 
ggsave("cd8t_caspgran_checkpoint_barplot_follicle.pdf", height=3, width=4)



### Fig 2g-h: Boxplot of CD8T density in follicles
#feature = "CD8T_density_follicle"
feature = "CD8T_density_extrafollicular"
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.7) +
  labs(y = "Cell density (cells/mm2)",
       title = "Extrafollicular") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(feature, "_statuswithp24.pdf"), height=4, width=3.75)

# Add p-value
pairwise_df_out = pval_p24_groups(feature)
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
ggplot(feature_tab, aes(x=status_with_p24, y=.data[[feature]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, size=0.7) +
  theme_bw() +
  labs(x = NULL,
       y = "Cell density (cells/mm2)",
       title = "Extrafollicular") +
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
ggsave(paste0(feature, "_statuswithp24_withpval.pdf"), height=4, width=3.75)



### Fig 2j: Heatmap of functional marker proportion in CD8T
keep_markers = c("CXCR5","GranzymeB","PD1","IFNg","TIGIT","NLRP3","Ki67","CD45RO","TIM3","Caspase1","ICOS","Lag3","Glut1","Galectin9","TCF1TCF7")

combined_tab = functional_tab[cell_tab_follicle, on=c("fov","label","cell_meta_cluster")]
combined_cd8t = combined_tab[cell_meta_cluster == 'CD8T']
# Get total number of CD8T cells in each compartment
total_cd8t_counts = combined_cd8t[,.N,by=.(fov, follicle_mask_label)]
# Get sum of functional marker+ CD8T cells in each compartment
counts_fov = combined_cd8t[,lapply(.SD, sum), by=.(fov, follicle_mask_label, in_follicle_mask), .SDcols=keep_markers]
# Get proportions
counts_fov = counts_fov[total_cd8t_counts, on=.(fov, follicle_mask_label)]
counts_fov[, (keep_markers) := lapply(.SD, function(x) x/N), .SDcols = keep_markers]
counts_fov = metadata[counts_fov, on=.(fov)]

# Get median of each group
compare = counts_fov[, lapply(.SD, median), by = .(status_with_p24, in_follicle_mask), .SDcols = keep_markers]
compare[, group := paste0(status_with_p24,"_",in_follicle_mask)]
# 99.9th percentile normalization
compare_norm = compare[, lapply(.SD, function(x) x / quantile(x[x!=0], 0.999)), .SDcols = keep_markers]
# Make heatmap
compare_df = data.frame(compare_norm[,..keep_markers])
rownames(compare_df) = compare$group
compare_df = t(compare_df)

pdf("cd8t_functional_inside_outside_follicle_heatmap.pdf", height=6, width=4)
pheatmap(compare_df,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#512DA8"))(100),
         breaks = seq(0, 1, length.out = 101))
dev.off()
