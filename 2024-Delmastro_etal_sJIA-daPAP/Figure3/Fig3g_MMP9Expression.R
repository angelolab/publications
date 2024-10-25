# Fig3g_MMP9Expression.R
# Author: Erin McCaffrey 
# Date updated: 07/15/2024

library(ggridges)
library(ggplot2)
library(reshape)

##..Import data..##
setwd("/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points")
data<-read.csv("celltable_05032024.csv")
labels<-read.csv('fov_labels.csv')
colorkey<-read.csv("cellpheno_num_colorkey.csv")

##..Include only PAP involved and uninvolved regions..##
keep<-labels[labels$FOV_Category %in% c("sJIA-daPAP", "non-sJIA-PAP", "uninvolved"),]$point
data <- data[data$point %in% keep,]

##..Drop giant cells..##
data<-droplevels(data[!data$name == 'Giant_cell',])

##..Melt..##
data.m <- melt(data, id.vars = 'name', measure.vars = 'MMP9')
colnames(data.m)<-c('name','MMP9','expression')

##..Get colors..##
color<-as.vector(colorkey_cells$Hex)
names(color) <- colorkey_cells$Pheno

##..Remove cells with zero MMP9 expression..##
data.m <- data.m[data.m$expression != 0,]

##..Get median..##
medians <- data.m %>%
  group_by(name) %>%
  summarize(median_value = median(expression)) %>%
  arrange(median_value)

data.m$name <- factor(data.m$name, levels = medians$name)

##..Ridge plot..##
ggplot(data.m, aes(x = expression, y = reorder(name, -expression), fill = name)) +
  geom_density_ridges2(alpha = 0.75) +
  scale_fill_manual(values = color) +
  geom_vline(xintercept = median(data.m[data.m$expression>0,]$expression), linetype = "dashed") +
  theme_ridges() + 
  theme(legend.position = "none")

ggsave('/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points/Sub-panels/Figure3/Fig3g_MMP9expression.pdf')

