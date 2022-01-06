# TBMIBIquantifyPDL1PosFreqAcrossMEsAndCells.R
# Author: Erin McCaffrey
# Date created: 200106
# Script reads in the ME-annotated data. It breaks down the proportion of PD-L1+ cells by FOV and 
# ME, plotting various examples/visualizations.

library(ggpubr)
library(ggplot2)
library(colorspace)
library(forcats)
library(dplyr)
library(tibble)
library(reshape2)
library(gplots)

##...Load in data..##

ME_data<-read.csv('data/allTB-ME_annotated.csv')

##...Filter non-myeloid cell types...##

cell_types<-c("CD206_Mac","CD209_DC","giant_cell","CD14_Mono","CD11c_DC/Mono","CD11b/c_CD206_Mac/Mono","CD16_CD14_Mono",
              "CD68_Mac","CD163_Mac")
ME_data <- droplevels(ME_data[ME_data$cell_type %in% cell_types, ])

##...Get percent of cells positive for PDL1 and ratio of pos:neg...##

marker_thresh <- 0.25
data_marker<-ME_data[ME_data$PD.L1>marker_thresh, ]

#Across samples
freq_sample<-as.data.frame(table(ME_data$SampleID, ME_data$maxME, ME_data$cell_type))
freq_marker_sample<-as.data.frame(table(data_marker$SampleID, data_marker$maxME, data_marker$cell_type))
freq_sample$PD.L1pos<-0
freq_sample[freq_sample$Var1 %in% freq_marker_sample$Var1 & 
              freq_sample$Var2 %in% freq_marker_sample$Var2 &
              freq_sample$Var3 %in% freq_marker_sample$Var3,]$PD.L1pos <- freq_marker_sample$Freq
names(freq_sample)<-c("SampleID","ME","cell_type","Total","Total_PDL1pos")
freq_sample$Total_PDL1neg <- freq_sample$Total - freq_sample$Total_PDL1pos
freq_sample$percentPDL1pos<-as.numeric(format((freq_sample$Total_PDL1pos / freq_sample$Total)*100),digits=3)
freq_sample$FC <- log2(freq_sample$Total_PDL1pos/freq_sample$Total_PDL1neg)

##...Visualize plot of FC across MEs (all cell types, boxplot)...##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

ggplot(data = freq_sample, aes(x = ME, y = FC, fill = ME)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'ME') + 
  labs(y = 'log2(PDL1+ / PDL1-)') +
  scale_fill_manual(values = colors)

##...Produce table of frequency of PDL1+ for all myeloid cells and then each ME...##

# across MEs and points
freq_ME<-as.data.frame(table(ME_data$maxME, ME_data$SampleID))
freq_marker_ME<-as.data.frame(table(data_marker$maxME, data_marker$SampleID))
freq_ME$PD.L1pos<-0
freq_ME[freq_ME$Var1 %in% freq_marker_ME$Var1 & 
              freq_ME$Var2 %in% freq_marker_ME$Var2,]$PD.L1pos <- freq_marker_ME$Freq
names(freq_ME)<-c("ME","SampleID","Total","Total_PDL1pos")

# total myeloid
freq_myeloid<-as.data.frame(table(ME_data$SampleID))
freq_myeloid<-add_column(freq_myeloid, ME = 'total_myeloid', .before = "Var1")
freq_myeloid_marker<-as.data.frame(table(data_marker$SampleID))
freq_myeloid$PD.L1pos<-0
freq_myeloid[freq_myeloid$Var1 %in% freq_myeloid_marker$Var1,]$PD.L1pos <- freq_myeloid_marker$Freq
names(freq_myeloid)<-c("ME","SampleID","Total","Total_PDL1pos")

# merge and get frequency
freq_ME<-rbind(freq_myeloid, freq_ME)
freq_ME$percentPDL1pos<-as.numeric((freq_ME$Total_PDL1pos / freq_ME$Total))
freq_ME$Total_PDL1neg <- freq_ME$Total - freq_ME$Total_PDL1pos
freq_ME$FC <- log2(freq_ME$Total_PDL1pos/freq_ME$Total_PDL1neg)


##...Plot the frequency and the FC as bar plots...##


data_summary <- freq_ME %>%
  group_by(ME) %>%
  summarize(mean = mean(percentPDL1pos), 
            n = n(), 
            sd = sd(percentPDL1pos), 
            se = sd/sqrt(n))

line<-data_summary[data_summary$ME=="total_myeloid",]$mean

# reoder to have total myeloid first
order<-c('total_myeloid',0,1,2,3,4,5,6,7)
freq_ME$ME<-factor(freq_ME$ME, levels = order)
freq_ME<-freq_ME[order(freq_ME$ME),]


ggplot(data = freq_ME, aes(x = ME, y = percentPDL1pos, fill = ME)) + 
  stat_summary(geom = "bar", fun = mean) +
  geom_point() +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3) +
  geom_hline(yintercept=line, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'ME') + 
  labs(y = '% PDL1+') +
  scale_fill_manual(values = c('#D3D3D3', colors))

##...Produce heatmap of the percent PDL1+ for a selected population across MEs...##

#by cell type
freq_ME_cell<-as.data.frame(table(ME_data$maxME, ME_data$cell_type))
freq_marker_cell<-as.data.frame(table(data_marker$maxME, data_marker$cell_type))
freq_ME_cell$PD.L1pos <- freq_marker_cell$Freq
names(freq_ME_cell)<-c("ME","cell_type","Total","Total_PDL1pos")
freq_ME_cell$percentPDL1pos<-as.numeric(format((freq_ME_cell$Total_PDL1pos / freq_ME_cell$Total)),digits=3)

#total myeloid
freq_ME_pooled<-as.data.frame(table(ME_data$maxME))
freq_marker_pooled<-as.data.frame(table(data_marker$maxME))
freq_ME_pooled$PD.L1pos <- freq_marker_pooled$Freq
names(freq_ME_pooled)<-c("ME","Total","Total_PDL1pos")
freq_ME_pooled$percentPDL1pos<-as.numeric(format((freq_ME_pooled$Total_PDL1pos / freq_ME_pooled$Total)),digits=3)

#add total myeloid data to the 
freq_ME_pooled<-add_column(freq_ME_pooled, d = 'myeloid', .after = "ME")
colnames(freq_ME_pooled)<-c("ME","cell_type","Total","Total_PDL1pos","percentPDL1pos")
freq_ME_cell<-rbind(freq_ME_pooled,freq_ME_cell)

# Plot the frequency of PDL1+ cell across all 
plot_data<-droplevels(freq_ME_cell[freq_ME_cell$cell_type=='CD163_Mac',])
plot_data<-droplevels(plot_data[plot_data$ME %in% c(1,2,3,4),])

ggplot(data = plot_data, aes(x = ME, y = percentPDL1pos, fill = ME)) + 
  geom_bar(stat='identity') +
  geom_hline(yintercept=0.2852326, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'ME') + 
  labs(y = '% PDL1+') +
  scale_fill_manual(values = c('#CA6627', '#7470AE', '#D53E88', '#74A439'))


#turn to heatmap format
cell_ME_PDL1_hmap<-dcast(freq_ME_cell, ME ~ cell_type, value.var = "percentPDL1pos")
cell_ME_PDL1_hmap<-as.matrix(cell_ME_PDL1_hmap[,-1])
cell_ME_PDL1_hmap[is.na(cell_ME_PDL1_hmap)] <- 0
rownames(cell_ME_PDL1_hmap)<-c(0,1,2,3,4,5,6,7)

heatmap.2(t(cell_ME_PDL1_hmap), 
          Colv = F, Rowv = F,
          dendrogram = 'none',
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          sepcolor="grey35",
          colsep=0:ncol(cell_ME_PDL1_hmap),
          rowsep=0:nrow(cell_ME_PDL1_hmap),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '',
          ColSideColors = colors,
          cexRow = 1, cexCol = 2, margins = c(8,14))

