# Liu and Bosse et al., Reproducible, high-dimensional imaging in archival human tissue by Multiplexed Ion Beam Imaging by Time-of-Flight (MIBI-TOF)
# Create cell tSNE in Supplementary Figure 7C

library(data.table)
library(Rtsne)
library(ggplot2)
library(stringr)
library(ggrastr)

set.seed(1024)

# Read in cell table that has been clustered
cell_dat = fread("cell_table_size_normalized_clusters.csv")
cell_dat[,slide := str_extract(fov,"(?<=Slide)\\d+")]
cell_dat[,core := str_extract(fov,"(?!.*_)\\w+")]

# Points to include in analysis
good_fovs = c('R1C2', 'R1C3', 'R1C5', 'R1C7', 'R1C9', 'R1C10',
             'R2C10', 'R2C11', 'R2C12', 'R3C2', 'R3C4', 'R6C7', 'R6C10',
             'R6C11', 'R7C6', 'R7C7', 'R7C10', 'R8C1', 'R8C10', 'R8C11')
cell_dat = cell_dat[core %in% good_fovs]

# List of markers
markers = c("CD11c","CD20","CD3","CD31","CD4","CD45","CD56","CD68","CD8","HLA DR","HLA class 1 A, B, and C, Na-K-ATPase alpha1","PANCK","PAX5","Vimentin","beta-tubulin","dsDNA")
marker_dat = cell_dat[,..markers]

# 99.9% normalization
quantile_norm <- function(x) {
  q = quantile(x,0.999)
  if (q!=0) {
    return(x/q)
  } else {
    return(x/max(x))
  }
 }
marker_dat = marker_dat[, lapply(.SD, function(x) quantile_norm(x))]

# Colors for plotting
pheno_colors = data.table(phenotype=c("B_cell","Epithelial","Endothelial","T_cell","Macrophage","Fibroblast","DC","NK_cell","Other"),
                          color=c("#A0CBE8","#F1CE63","#8CD17D","#FF9D9A","#B07AA1","#499894","#9D7660","#FFBE7D","#808080"))
core_colors = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#FBB4AE","#B3CDE3","#CCEBC5")

# tSNE
tsne_out = Rtsne(as.matrix(marker_dat))
tsne_dims = data.table(tsne_out$Y)
colnames(tsne_dims) = c("tsne1","tsne2")
all_dat = cbind(cell_dat[,c("phenotype","label","fov","slide","core")],tsne_dims)

# Color by phenotype
ggplot(all_dat, aes(x=tsne1, y=tsne2, color=phenotype)) +
  geom_point_rast(size=0.1) +
  scale_color_manual(breaks=pheno_colors$phenotype, values=pheno_colors$color) +
  xlab("TSNE1") +
  ylab("TSNE2") +
  theme_void()
ggsave("cell_tsne_pheno.pdf",height=8,width=9)

# Color by slide
ggplot(all_dat, aes(x=tsne1, y=tsne2, color=slide)) +
  geom_point_rast(size=0.1, alpha=0.5) +
  scale_color_brewer(palette="Set2") +
  xlab("TSNE1") +
  ylab("TSNE2") +
  theme_void()
ggsave("cell_tsne_slide.pdf",height=8,width=9)

# Color by ROI
ggplot(all_dat, aes(x=tsne1, y=tsne2, color=core)) +
  geom_point_rast(size=0.1, alpha=0.5) +
  scale_color_manual(breaks=good_fovs, values=core_colors) +
  xlab("TSNE1") +
  ylab("TSNE2") +
  theme_void()
ggsave("cell_tsne_core.pdf",height=8,width=9)


