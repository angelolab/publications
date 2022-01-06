# TBMIBIquantifyPD1LymphsAcrossMEs.R
# Author: Erin McCaffrey
# Date created: 200225
# Overview: This script reads in the csv for the ME loadings for each cell with each cell assigned to a ME. 
# Next, it looks at the distribution of PD1+ lymphs across MEs in bulk.

library(ggplot2)
library(forcats)
library(dplyr)

##...Load in data..##

ME_data<-read.csv("data/allTB-ME_annotated.csv")

##...Indicate cell type(s) to look at...##
cell_types<-c('lymphocyte')

##...Create dataframe with just cell type of interest...##
ME_data_cell <- ME_data[ME_data$cell_lin %in% cell_types, ]

##..Counts of all PD1 positive cells across MEs..##

PD1_pos<-as.data.frame(table(ME_data_cell[ME_data_cell$PD.1>0.21, ]$maxME))
colnames(PD1_pos)<-c('ME','Count')

ggplot(data = PD1_pos, aes(x = as.factor(ME), y = Count, fill = as.factor(ME))) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = 'ME') + 
  labs(y = 'Count PD1+ Cells') +
  scale_fill_manual(name = 'ME', values = colors) 

