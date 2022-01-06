# TBMIBIrunUMAP.R
# Author: Erin McCaffrey
# Date created: 190116
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it runs UMAP on chosen data and visualizes results colored by cell type (or any other selected variable)

library("umap")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library("ggrastr")
library(viridis)

##..Import data..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")

# TB samples
data_myco<-droplevels(data[data$Tissue %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])

##..Choose clustering channels..##

clusterChannels = c('CD14','CD16','CD68','CD163','CD11b','CD11c','CD209','CD206','CD45','MPO','MastChyTry')

##..Prepare data with only myeloid cells..##

myeloid<-c("CD206_Mac","CD209_DC","giant_cell","CD14_Mono","CD11c_DC/Mono","CD68_Mac",
           "CD16_CD14_Mono","CD163_Mac","CD11b/c_CD206_Mac/Mono","mast","neutrophil")
data_myeloid<-droplevels(data_myco[data_myco$cell_type %in% myeloid,])
data_myeloid[is.na(data_myeloid)]<-0

##..Run UMAP..##

set.seed(4)
gran_myeloid.umap<-umap(data_myeloid[,clusterChannels], min_dist = 0.5)
umap.output2<-gran_myeloid.umap$layout
data_myeloid$umap1<-umap.output2[,1]
data_myeloid$umap2<-umap.output2[,2]

##..Visualize..##

#Color by sample ID, tissue, or cell ID
data_myeloid$cell_type<-factor(data_myeloid$cell_type, levels=myeloid)
data_myeloid<-data_myeloid[order(data_myeloid$cell_type),]

##..Plot..##
colorkey<-read.csv("data/colorkey_R.csv")
myeloid_color<-droplevels(colorkey[colorkey$imm_order %in% myeloid,])
myeloid_color$imm_order<-factor(myeloid_color$imm_order, levels = myeloid)

myeloid_color<-myeloid_color[order(myeloid_color$imm_order),]
color<-as.vector(myeloid_color$code)

ggplot(data_myeloid, aes(umap1, umap2, color = as.factor(cell_type))) +
  geom_point_rast(size=2)  + 
  theme_bw() + 
  scale_color_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="UMAP 1") +
  labs(y="UMAP 2") +
  labs(color = "Immune Phenotype") +
  ggtitle("UMAP of Mycobacterial Myeloid Cells")


#overlay marker expression 

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())

data_myeloid.melt <- reshape2::melt(data_myeloid, id.vars = c('umap1', 'umap2'), measure.vars = c(clusterChannels,
                                                                                                  'IDO','PD.L1','HLA.DR.DQ.DP'))
data_myeloid.melt[data_myeloid.melt$value >1, ]$value = 1
# randomize row order so that no one patient is visualized more
rows <- sample(nrow(data_myeloid.melt))
data_myeloid.melt.rand <- data_myeloid.melt[rows, ]

ggplot(data_myeloid.melt.rand, aes(x = umap1, y = umap2, color = value)) +
  geom_point_rast(size=3) +
  scale_color_viridis(option = 'magma', limits=c(0, 1)) +
  facet_wrap(~ variable, ncol = 4) +
  theme +
  theme(legend.position = 'none')

