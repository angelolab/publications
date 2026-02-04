# Code for Figure 1, Supplementary Figure 4 
# Author: Candace Liu

library(data.table)
library(ggplot2)
library(pheatmap)
library(viridis)

cell_tab = fread("../data/tables/cell_table_size_normalized.csv")
metadata = fread("../data/tables/metadata.csv")
cell_colors = fread("../data/hiv_colors.csv")
tissue_area = fread("../data/tables/tissue_area.csv")

# Make color mappings
cell_color_mapping = setNames(cell_colors$color, cell_colors$phenotype)
sample_mapping = c("#F19A8A","#99B3D9","#B2D1AB","#C4A1CA","#FEC58C","#F2F2B2","#DBCEA3","#F9C0D2","#DEDEDE","#8FD4B3","#E6F5C9")
names(sample_mapping) = 1:10
status_mapping = c("hiv_pos"="#BF616A", "hiv_neg"="#5E81AC")
color_mapping = c(cell_color_mapping, sample_mapping, status_mapping)



### Fig 1c: Cell phenotype heatmap
# Blue-white-red colors
rwb_cols = colorRampPalette(c("royalblue4","white","red4"))(99)
# Make colors for heatmaps
mat_colors = cell_colors$color
names(mat_colors) = cell_colors$phenotype
mat_colors = list(pheno = mat_colors)

markers = c("Calprotectin","CD3e","CD4","CD8a","CD11c","CD14","CD20",
            "CD21","CD31","CD45","CD56","CD68","CD163","CXCR5","FOXP3",
            "HLADR","MastCellTryptase","PD1","SMA")
clusters_copy = copy(cell_tab)
setnames(clusters_copy,"cell_meta_cluster","phenotype")
# Combine immune_other, unassigned, and other for visualization
clusters_copy[phenotype == "Immune_other", phenotype := "Other"]
clusters_copy[phenotype == "Unassigned", phenotype := "Other"]

mean_dat = clusters_copy[, lapply(.SD, mean), by = phenotype, .SDcols=markers]
mat_dat = data.frame(mean_dat[,..markers])
rownames(mat_dat) = mean_dat$phenotype
# Z-score and cap
mat_dat = scale(mat_dat)
mat_dat = pmin(mat_dat, 3)
# Annotations
mat_col = data.frame(pheno = mean_dat$phenotype)
rownames(mat_col) = mean_dat$phenotype
# Make heatmap
breaks = seq(-3, 3, length.out=100)
pdf("cells_heatmap.pdf", height=8, width=10)
pheatmap(mat_dat,
         color = rwb_cols,
         breaks = breaks,
         annotation_row = mat_col,
         annotation_colors = mat_colors,
         main = "Cell phenotypes")
dev.off()



### Sup Fig 4a: Cell phenotype heatmap (no normalization)
mat_dat = data.frame(mean_dat[,..markers])
rownames(mat_dat) = mean_dat$phenotype
# Annotations
mat_col = data.frame(pheno = mean_dat$phenotype)
rownames(mat_col) = mean_dat$phenotype
# Make heatmap
pdf("cells_heatmap_noz.pdf", height=8, width=10)
pheatmap(mat_dat,
         color = viridis(99),
         breaks = seq(0,0.003,length.out=100),
         annotation_row = mat_col,
         annotation_colors = mat_colors,
         main = "Cell phenotypes")
dev.off()



### Supp Fig 4b: Cell phenotype quantification
# Count number of each cell type in each fov
count_tab = cell_tab[,.N, by=.(fov,cell_meta_cluster)]
count_tab = merge(count_tab, metadata, on='fov')
count_tab = count_tab[order(sample_id)]
count_tab$sample_id = factor(count_tab$sample_id)
count_tab$fov = factor(count_tab$fov, levels=unique(metadata$fov))
# Stacked bar plot (proportion)
ggplot(count_tab, aes(x = fov, y = N, fill=cell_meta_cluster)) +
  geom_bar(position = "fill", stat="identity") +
  scale_fill_manual(values = color_mapping) +
  geom_tile(data = count_tab, aes(x = fov, y = -0.02, fill = sample_id), height = 0.02) +
  geom_tile(data = count_tab, aes(x = fov, y = -0.05, fill = status), height = 0.02) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.08,NA)) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("celltype_proportions.pdf", height=8, width=18)

# Pie chart
all_counts = cell_tab[,.N, by=.(cell_meta_cluster)]
ggplot(all_counts, aes(x = "", y = N, fill = cell_meta_cluster)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = color_mapping) +
  coord_polar("y", start = 0) +
  theme_void()
ggsave("celltype_proportions_pie.pdf", height=8, width=8)

