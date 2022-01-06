# TBMIBIexprsViolinPlot.R
# Author: Erin McCaffrey 
# Date created: 190415
# Overview: This script reads in the csv for normalized data and then visualizes violin plots of expression 
# values across the cell types for selected markers.

library(ggplot2)
library(dplyr)
library(forcats)

##..Read in data..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Subset just the TB samples and mono, mac, DC..##

data_myco<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])
myeloid<-c("CD14_Mono","CD11b/c_CD206_Mac/Mono","CD11c_DC/Mono","CD68_Mac","CD16_CD14_Mono",
           "CD206_Mac","CD163_Mac","CD209_DC","mast","neutrophil")
data_myeloid<-droplevels(data_myco[data_myco$cell_type %in% myeloid, ])


colorkey<-read.csv("data/colorkey_R.csv")
myeloid_color<-droplevels(colorkey[colorkey$imm_order %in% myeloid,])

##..Plot violin plot of expression across cell types IDO..##

# Get order for IDO
data_myeloid$cell_type= fct_reorder(.f = data_myeloid$cell_type,.x = data_myeloid$IDO,.fun = median,.desc = TRUE)
#reorder levels of data to reflect IDO order
IDO_order<-levels(data_myeloid$cell_type)
data_myeloid$cell_type = factor(data_myeloid$cell_type, level= IDO_order)

# get color in proper order
myeloid_color$imm_order<-factor(myeloid_color$imm_order, levels=IDO_order)
myeloid_color<-myeloid_color[order(myeloid_color$imm_order),]
IDO_color<-as.vector(myeloid_color$code)

IDO<-ggplot(data_myeloid, aes(x=fct_reorder(cell_type,IDO,.fun=median,.desc=TRUE), y=IDO, fill = cell_type)) +
  geom_violin(trim=FALSE, draw_quantiles = 0.5, width = 1) +
  geom_hline(yintercept=c(0.26), linetype="dashed") +
  scale_fill_manual(values = IDO_color) +
  theme_bw() +
  theme(axis.ticks.x=element_blank()) +
  theme(legend.position = "none") +
  labs(x="Cell Type") +
  labs(y="IDO1 Expression") 
IDO


##..Plot violin plot of expression across cell types IDO..##

# Get order for PDL1
data_myeloid$cell_type= fct_reorder(.f = data_myeloid$cell_type,.x = data_myeloid$PD.L1,.fun = median,.desc = TRUE)
PDL1_order<-levels(data_myeloid$cell_type)
data_myeloid$cell_type = factor(data_myeloid$cell_type, level= PDL1_order)

# get color in proper order
myeloid_color$imm_order<-factor(myeloid_color$imm_order, levels=PDL1_order)
myeloid_color<-myeloid_color[order(myeloid_color$imm_order),]
PDL1_color<-as.vector(myeloid_color$code)

PDL1<-ggplot(data_myeloid, aes(x=fct_reorder(cell_type,PD.L1,.fun=median,.desc=TRUE), y=PD.L1, fill = cell_type)) + 
  geom_violin(trim=FALSE, draw_quantiles = 0.5, width = 1) +
  geom_hline(yintercept=c(0.25), linetype="dashed") +
  scale_fill_manual(values = PDL1_color) +
  theme_bw() +
  theme(axis.ticks.x=element_blank()) +
  theme(legend.position = "none") +
  labs(x="Cell Type") +
  labs(y="PDL1 Expression") 
PDL1


