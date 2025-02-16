## MIBI_plot_cell_heatmap.R ##
# Author: Erin McCaffrey 
# Date created: 231204



library(pals)
library(dplyr)
library(ggpubr)
library(forcats)
library(gplots)
library(RColorBrewer)

##..Import data..##

# Read in necrotic and non-necrotic separately
setwd("/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2")
cell_table<-read.csv("cohort_cell_table.csv")

##..Make heatmap based on marker expression to validate..##

# go through all clusters and calculate mean for every channel
cluster_channels<-c("CD3", "CD4", "CD8a", "CD11c", "CD14", "CD20", "CD31", "CD45", "CD68", "Vimentin",
                    "CD163", "CD206", "Calprotectin", "Chy_Try", "SMA", "Fibronectin", "HLA.DR")
cluster_names <-unique(cell_table$pheno_corrected)

hm_allclusters <- matrix(, nrow = length(cluster_names), ncol = length(cluster_channels))
for(i in 1:length(cluster_names)) {
  temp_mat <- cell_table[cell_table[,"pheno_corrected"] == cluster_names[i], cluster_channels]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
rownames(hm_allclusters) <- cluster_names
colnames(hm_allclusters) <- cluster_channels
hm_allclusters

# get color key
order <- rownames(hm_allclusters)
colorkey<-read.csv('./keys/cell_color_key.csv')
colorkey<-droplevels(colorkey[-1,])
colorkey$Pheno<-factor(colorkey$Pheno, levels = order)
colorkey<-colorkey[order(colorkey$Pheno),]
color<-as.vector(colorkey$Hex)


# custom_pal<-c("#150151","#19165E","#1E2969","#243A73","#2A4B7D","#325C87",
#                        "#3A6C92","#437E9D","#4D8FA7","#5EA0AF","#82C3BF",
#                        "#97D3CA","#B4E1D7",'#D1F0E5','#EEFEF2')
custom_pal_div<-c("#0E1949","#1E356C","#31508C","#4272AE","#6A93C6","#98B1DA",
                           "#C8D0EF","#F8F0FE","#F0C5D8","#E19EB0","#D17486",
                           "#BD4B5C","#923346","#691D32","#43071E")
# custom_pal_simple<-c('#FFF',"#5EA0AF")
colfunc<-colorRampPalette(custom_pal_div)

# plot heatmap of all metaclusters
heatmap.2(hm_allclusters,
          scale = "column",
          Colv = T, Rowv = T,
          hclustfun = function(x) hclust(x, method="complete"),
          dendrogram = "both",
          trace = "none",
          col = colfunc(100),
          # # col = rev(as.vector((brewer.rdbu(100)))),
          # col = rev(as.vector((brewer.spectral(100)))),
          # col = (as.vector((ocean.delta(100)))),
          RowSideColors = color,
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),      
          sepcolor="white",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          breaks=seq(-3,3, length.out=101))
