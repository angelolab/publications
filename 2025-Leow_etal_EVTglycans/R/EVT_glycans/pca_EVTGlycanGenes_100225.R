# PCA of EVT glycans and genes
# Author: Ke Leow
# Date: 05/13/24

#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
#load PCA function - 
calculatePC <- function(exprs.t, pheno) {
  variances<-apply(exprs.t, 2, var)  # calculates the variances
  rejected<-(variances==0 | is.na(variances))
  exprs.t2<-exprs.t[, !rejected]
  exprs.model<<-prcomp(exprs.t2, center=TRUE, scale.=TRUE)
  #write.table(exprs.model$x, file=scoresOutfile, sep="\t", col.names=NA, quote=F)
  #write.table(exprs.model$rotation, file=loadingsOutfile, sep="\t", col.names=NA, quote=F)
  v<<-round((exprs.model$sdev^2)/sum(exprs.model$sdev^2)*100,digits=1)
  #plot(v, type="b", ylab="Variation", xlab="Principal component")
  #text(v, labels=v, cex= 0.7, pos=3)
  scores<-exprs.model$x
  data.melt<-data.frame(PC1=scores[,1], PC2=scores[,2], PC3=scores[,3], PC4=scores[,4], pheno)
  data.melt
}
#--------------------------------
# Load data
#--------------------------------
#glycan+gene zscore data - with kmeans clusters
data <- read_csv("data/EVT_Nanostring/combined_gly_enz_heatmap_nopartialgal_type_list_K_4.csv") %>%
  mutate(Feature = if_else(startsWith(gene_name, "H"), "Glycan", "Gene")) %>% 
  rename(
    PV   = mean_by_stage_norm1,
    pEVT = mean_by_stage_norm2,
    iEVT = mean_by_stage_norm3,
    eEVT = mean_by_stage_norm4
  )%>% 
  mutate(k_cluster_new = case_when(
    k_cluster == 3 ~ 2,
    k_cluster == 2 ~ 3,
    TRUE ~ k_cluster
  )) %>%
  arrange(k_cluster_new)%>% 
  column_to_rownames("gene_name")

data_int <- data %>% select_at(vars(contains("V", ignore.case=FALSE))) 

data_pheno <- data[,c("Feature", "k_cluster_new")]
data_pheno$Feature <- factor(data_pheno$Feature, levels = c("Glycan", "Gene"))
data_pheno$k_cluster_new <- as.factor(data_pheno$k_cluster_new)
data_pheno$gene_name <- rownames(data_pheno)


# #check heatmap
# my_gene_col <- data_pheno[,c("Feature", "Cluster")]
# 
# pheatmap(data_int,
#          annotation_row = my_gene_col,
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          # show_rownames=FALSE,
#          color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                    "RdBu")))(100))

#--------------------------------
# Plot PCA 
#--------------------------------
pc <- calculatePC(data_int, data_pheno)
loadings <- exprs.model$rotation



#scatter plot of scores
pdf("R_plots/MALDI_IF_EVT_glycans/pca_EVTglycanGenes_100225.pdf", width = 4, height = 4)
pc %>% 
  ggplot(aes(x=PC1, y=PC2, colour=k_cluster_new, shape = Feature)) + 
  scale_color_manual(values = c("1" = "#FF8811FF", "2" = "#462255FF", "3" = "#D44D5CFF", "4" = "#046E8FFF"))+
  scale_shape_manual(values=c(1, 17))+
  geom_point(size=3)+ 
  xlab(paste("PC1 (", v[1], "%)", sep="")) +
  ylab(paste("PC2 (", v[2], "%)", sep="")) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims))+
  # geom_label_repel(aes(label = gene_name), na.rm = TRUE,
  #                  box.padding = 0.4, point.padding = 0.3,
  #                  segment.color = 'grey50', max.overlaps = 20,
  #                  fill = "white", color = "black", size = 3) +
  theme_classic(base_size = 18)+ 
  theme(
    # text=element_text(size=18),
        legend.position = "none"
        )
dev.off()
# save_plot("figure/pca_everyGly_glyRatio.svg", width=10, height=8)
# ggsave(file="figures/pca_everyGly_glyRatio.svg", width=4, height=4)

#--------------------------------
# Plot PCA - color by glycan types
#--------------------------------
#load gene and glycan annotations
glycan_annotation <- read_csv("data/MALDI_IF_EVT_glycans/glycan_peaklist_typeAnnotations_051524.csv")
gene_annotation <- read_csv("data/MALDI_IF_EVT_glycans/glycan_types_enzymes_051524.csv")%>%
  rename(gene_name = 1)%>%
  mutate(gene_name = str_to_upper(gene_name))

#combine gene and glycan annotations
glycan_gene_anno <- glycan_annotation %>%
  select(-mz) %>% rename(gene_name = 1)  %>%
  rbind(gene_annotation)

#merge with pc data
pc_anno <- inner_join(pc, glycan_gene_anno)

# Use the "YlOrRd" color palette from RColorBrewer
blue_colors <- brewer.pal(5, "Blues")
purple_colors <- brewer.pal(5, "Purples")
red_colors <- brewer.pal(4, "Reds")


#####################
#color by type - biantennary
selected_glycan <- pc_anno %>%
    filter(biantennary == 1 & Feature =="Glycan") %>%
    pull(gene_name)

type= "biantennary"
selected_gene = c("MGAT1", "MGAT2") #add gene
pc_anno_update <- pc_anno %>% 
  mutate(annotate = ifelse(gene_name %in% c(selected_glycan, selected_gene), type, NA))

pc_anno_update$label_genes <- ifelse(pc_anno_update$gene_name %in% selected_gene,
                                     pc_anno_update$gene_name, NA)
# Find overall range across PC1 and PC2
lims <- max(abs(pc[,1:4]))

pdf("R_plots/MALDI_IF_EVT_glycans/pca_EVTglycanGenes_biantennary_100225.pdf", width = 4, height = 4)
pc_anno_update %>%
  filter(annotate == type) %>%
  ggplot(aes(x=PC1, y=PC2, shape = Feature), colour="black") +
  # scale_color_manual(values = col_decidua)+
  scale_shape_manual(values=c(1, 17))+
  geom_point(size=4)+
  scale_color_manual(values = c(blue_colors[3])) +
  xlab(paste("PC1 (", v[1], "%)", sep="")) +
  ylab(paste("PC2 (", v[2], "%)", sep="")) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims))+
  geom_text_repel(aes(label = label_genes),
                  box.padding = 2,      # Adjust distance between the label and the point
                  point.padding = 0.5,    # Adjust padding between text and the point
                  nudge_x = 0.4,          # Move labels slightly along x-axis
                  nudge_y = 0.1,
                  size=8
                  ) +
  theme_classic()+
  theme(text=element_text(size=18), 
        legend.position = "none"
  )
# +
  # ggtitle(type)  # Add plot title
dev.off()


#####################
#color by type - Tetra and Lacnac
selected_glycan <- pc_anno %>%
  filter((tetraantennary == 1 | polylacnac == 1) & Feature =="Glycan") %>%
  pull(gene_name)

type= "highly branched"
selected_gene = c("MGAT5", "B3GNT7", "B3GNT2") #add gene
pc_anno_update <- pc_anno %>% 
  mutate(annotate = ifelse(gene_name %in% c(selected_glycan, selected_gene), type, NA))

pc_anno_update$label_genes <- ifelse(pc_anno_update$gene_name %in% selected_gene,
                                     pc_anno_update$gene_name, NA)

pdf("R_plots/MALDI_IF_EVT_glycans/pca_EVTglycanGenes_highlybranched_100225.pdf", width = 4, height = 4)
pc_anno_update %>%
  filter(annotate == type) %>%
  ggplot(aes(x=PC1, y=PC2, shape = Feature), colour="black") +
  # scale_color_manual(values = col_decidua)+
  scale_shape_manual(values=c(1, 17))+
  geom_point(size=4)+
  scale_color_manual(values = c(blue_colors[5])) +
  xlab(paste("PC1 (", v[1], "%)", sep="")) +
  ylab(paste("PC2 (", v[2], "%)", sep="")) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims))+
  geom_text_repel(aes(label = label_genes),
                  box.padding = 1,      # Adjust distance between the label and the point
                  point.padding = 0.5,    # Adjust padding between text and the point
                  nudge_x = 0.1,          # Move labels slightly along x-axis
                  nudge_y = 0.1,
                  size=8
                  ) +
  theme_classic()+
  theme(text=element_text(size=18), 
        legend.position = "none"
  )
# +
#   ggtitle(type)  # Add plot title
dev.off()

# #####################
# #color by type - core fucose
# selected_glycan <- read_csv('data/musc_ta511_endoF3_PNGase_AAXL/summary_glycan_extracted_coreFucose_051325.csv') %>%
#   pull(composition)
# 
# type= "core fucose"
# selected_gene = "FUT8"#add gene
# pc_anno_update <- pc_anno %>%
#   mutate(annotate = ifelse(gene_name %in% c(selected_glycan, selected_gene), type, NA))
# 
# #summarize how many fucosylated
# pc_anno_update %>% group_by(fucosylated) %>% summarise(n=n())
# #summarize how many core fucosylated
# pc_anno_update %>% group_by(annotate) %>% summarise(n=n())
# 
# pc_anno_update$label_genes <- ifelse(pc_anno_update$gene_name %in% selected_gene,
#                                      pc_anno_update$gene_name, NA)
# 
# pdf("R_plots/MALDI_IF_EVT_glycans/pca_EVTglycanGenes_coreFucose.pdf", width = 4, height = 4)
# pc_anno_update %>%
#   filter(annotate == type) %>%
#   ggplot(aes(x=PC1, y=PC2, shape = Feature), colour="black") +
#   # scale_color_manual(values = col_decidua)+
#   scale_shape_manual(values=c(1, 17))+
#   geom_point(size=3)+
#   scale_color_manual(values = c(red_colors[4])) +
#   xlab(paste("PC1 (", v[1], "%)", sep="")) +
#   ylab(paste("PC2 (", v[2], "%)", sep="")) +
#   geom_text_repel(aes(label = label_genes),
#                   box.padding = 0.5,      # Adjust distance between the label and the point
#                   point.padding = 0.5,    # Adjust padding between text and the point
#                   nudge_x = 0.8,          # Move labels slightly along x-axis
#                   nudge_y = 0.1,
#                   size = 6
#   ) +
#   theme_classic()+
#   theme(text=element_text(size=18),
#         legend.position = "none"
#   )+
#   ggtitle(type)  # Add plot title
# dev.off()


# #####################
# #color by type - 2,3 sialyl
# selected_glycan <- read_csv('data/musc_ta511_endoF3_PNGase_AAXL/summary_glycan_extracted_SA_051325.csv') %>% 
#   filter(SA_23) %>% 
#   pull(AAXL_parent)
# 
# type= "2,3 sialyl"
# selected_gene = c("ST3GAL4", "ST3GAL6") #add gene
# pc_anno_update <- pc_anno %>% 
#   mutate(annotate = ifelse(gene_name %in% c(selected_glycan, selected_gene), type, NA))
# 
# #summarize how many sialyl
# pc_anno_update %>% group_by(sialylated) %>% summarise(n=n())
# #summarize how many core fucosylated
# pc_anno_update %>% group_by(annotate) %>% summarise(n=n())
# 
# pc_anno_update$label_genes <- ifelse(pc_anno_update$gene_name %in% selected_gene,
#                                      pc_anno_update$gene_name, NA)
# 
# pdf("R_plots/MALDI_IF_EVT_glycans/pca_EVTglycanGenes_23sialyl.pdf", width = 4, height = 4)
# pc_anno_update %>%
#   filter(annotate == type) %>%
#   ggplot(aes(x=PC1, y=PC2, shape = Feature), colour="black") +
#   # scale_color_manual(values = col_decidua)+
#   scale_shape_manual(values=c(1, 17))+
#   geom_point(size=3)+
#   scale_color_manual(values = c(purple_colors[5])) +
#   xlab(paste("PC1 (", v[1], "%)", sep="")) +
#   ylab(paste("PC2 (", v[2], "%)", sep="")) +
#   geom_text_repel(aes(label = label_genes),
#                   box.padding = 1,      # Adjust distance between the label and the point
#                   point.padding = 0.5,    # Adjust padding between text and the point
#                   nudge_x = 0.8,          # Move labels slightly along x-axis
#                   nudge_y = 0.1,
#                   size=6
#   ) +
#   theme_classic()+
#   theme(text=element_text(size=18), 
#         legend.position = "none"
#   )+
#   ggtitle(type)  # Add plot title
# dev.off()
