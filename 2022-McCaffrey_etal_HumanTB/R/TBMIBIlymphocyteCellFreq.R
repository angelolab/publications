# TBMIBIlymphocyteCellFreq.R
# Author: Erin McCaffrey 
# Date created: 190922
# Overview: This script reads in the csv for immune cell frequencies (out of total immune) across all 
# samples. Next it plots a piechart of the the lymphocyte totals across the mycobacterial samples.

library(dplyr)
library(viridis)
library(ggplot2)
library(tidyr)

##..Read in data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("allgran_cellpheno_noHIV.csv")

##..Subset just the mycobacterial samples and lymphs..##
data_myco<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])
lymph<-c("CD4_T","CD8_T","B_cell","Treg","gdT_cell")
data_lymph<-droplevels(data_myco[data_myco$cell_type %in% lymph, ])

##..Get a table of the major subsets and their total counts..##

lymph_freq<-as.data.frame(table(data_lymph$cell_type))
colnames(lymph_freq)<-c('cell_type','count')
lymph_freq$freq<-(lymph_freq$count / sum(lymph_freq$count)) *100
lymph_freq$cell_type<-factor(lymph_freq$cell_type, levels = lymph)
lymph_freq<-lymph_freq[order(lymph_freq$cell_type),]

##..Get appropriate color key..##

colorkey<-read.csv("colorkey_R.csv")
lymph_color<-droplevels(colorkey[colorkey$imm_order %in% lymph,])
lymph_color$imm_order<-factor(lymph_color$imm_order, levels = lymph)
lymph_color<-lymph_color[order(lymph_color$imm_order),]
color<-as.vector(lymph_color$code)

##..Plot..##

lymph_pie <- ggplot(lymph_freq, aes(x = "", y = freq, fill = cell_type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = color) +
  theme_void()
lymph_pie
