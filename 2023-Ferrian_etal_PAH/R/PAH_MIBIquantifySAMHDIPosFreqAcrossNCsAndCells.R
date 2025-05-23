# PAH_MIBIquantifySAMHDIPosFreqAcrossNCsAndCells.R
# Author: Erin McCaffrey
# Date created: 200909

library(ggpubr)
library(ggplot2)
library(colorspace)
library(forcats)
library(dplyr)
library(tibble)
library(reshape2)
library(gplots)
library(RColorBrewer)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
data<-read.csv("celldata_region_annotated.csv")
data_neighborhood<-read.csv("PAH_AllCells_Neighborhood_K=10_PatientAnnotated.csv")

##..Drop healthy..##

data_pah<-droplevels(data[data$Subgroup %in% c('HPAH','IPAH'),])

##..Merge the neighborhood data with the expression data..##

NC_data<-merge(data_neighborhood, data_pah, by=c("Point_num","label",'PID'))

##...Filter to only include immune cells...##

data_immune<-droplevels(NC_data[!NC_data$cell_lineage %in% c('epithelial', 'mesenchymal', 
                                                             'fibroblast', 'endothelial'),])

##...Get percent of cells positive for TIM3 and ratio of pos:neg...##

marker_thresh <- 0.39
data_marker<-droplevels(data_immune[data_immune$SAMHD1>marker_thresh, ])

#Across samples
freq_sample<-as.data.frame(table(data_immune$PID, data_immune$cluster, data_immune$cell_lineage))
freq_marker_sample<-as.data.frame(table(data_marker$PID, data_marker$cluster, data_marker$cell_lineage))
freq_sample$SAMHD1pos <- freq_marker_sample$Freq
names(freq_sample)<-c("PID","cluster","cell_lineage","Total","Total_SAMHD1pos")
freq_sample$Total_SAMHD1neg <- freq_sample$Total - freq_sample$Total_SAMHD1pos
freq_sample$percentSAMHD1pos<-as.numeric(format((freq_sample$Total_SAMHD1pos / freq_sample$Total)*100),digits=3)
freq_sample$FC <- log2(freq_sample$Total_SAMHD1pos/freq_sample$Total_SAMHD1neg)

##...Visualize plot of FC across clusters (all cell types, boxplot)...##


ggplot(data = freq_sample, aes(x = cluster, y = FC, fill = cluster)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'cluster') + 
  labs(y = 'log2(SAMHD1+ / SAMHD1-)') +
  scale_fill_brewer(palette = "Paired")

##...Produce table of frequency of SAMHD1+ for all immune cells and then each cluster...##

# across clusters and points
freq_cluster<-as.data.frame(table(data_immune$cluster, data_immune$PID))
freq_marker_cluster<-as.data.frame(table(data_marker$cluster, data_marker$PID))
freq_cluster$SAMHD1pos <- freq_marker_cluster$Freq
names(freq_cluster)<-c("cluster","PID","Total","Total_SAMHD1pos")

# total immune
freq_immune<-as.data.frame(table(data_immune$PID))
freq_immune<-add_column(freq_immune, cluster = 'total_immune', .before = "Var1")
freq_immune_marker<-as.data.frame(table(data_marker$PID))
freq_immune$SAMHD1pos <- freq_immune_marker$Freq
names(freq_immune)<-c("cluster","PID","Total","Total_SAMHD1pos")

# merge and get frequency
freq_cluster<-rbind(freq_immune, freq_cluster)
freq_cluster$percentSAMHD1pos<-as.numeric((freq_cluster$Total_SAMHD1pos / freq_cluster$Total))
freq_cluster$Total_SAMHD1neg <- freq_cluster$Total - freq_cluster$Total_SAMHD1pos
freq_cluster$FC <- log2(freq_cluster$Total_SAMHD1pos/freq_cluster$Total_SAMHD1neg)


##...Plot the frequency and the FC as bar plots...##

data_summary <- freq_cluster %>%
  group_by(cluster) %>%
  summarize(mean = mean(percentSAMHD1pos), 
            n = n(), 
            sd = sd(percentSAMHD1pos), 
            se = sd/sqrt(n))

line<-data_summary[data_summary$cluster=="total_immune",]$mean

# reoder to have total immune first
order<-c('total_immune',1,2,3,4,5,6,7,8,9,10)
freq_cluster$cluster<-factor(freq_cluster$cluster, levels = order)
freq_cluster<-freq_cluster[order(freq_cluster$cluster),]

colors<-brewer.pal(10, 'Paired')

ggplot(data = freq_cluster, aes(x = cluster, y = percentSAMHD1pos, fill = cluster)) + 
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3) +
  geom_point(data = freq_cluster, aes(x = cluster, y = percentSAMHD1pos)) +
  geom_hline(yintercept=line, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'cluster') + 
  labs(y = '% SAMHD1+ Immune Cells (Of Total Immune per NC)') +
  scale_fill_manual(values = c('#D3D3D3', colors))


stats<-compare_means(percentSAMHD1pos~cluster,freq_cluster, method = 'wilcox.test')
# write.csv(stats, 'SAMHD1pos_total-immune_perNC_stats.csv',row.names = F)

##...Produce heatmap of the percent SAMHD1+ for a selected population across clusters...##

#by cell type
freq_cluster_cell<-as.data.frame(table(data_immune$cluster, data_immune$cell_lineage))
freq_marker_cell<-as.data.frame(table(data_marker$cluster, data_marker$cell_lineage))
freq_cluster_cell$SAMHD1pos <- freq_marker_cell$Freq
names(freq_cluster_cell)<-c("cluster","cell_lineage","Total","Total_SAMHD1pos")
freq_cluster_cell$percentSAMHD1pos<-as.numeric(format((freq_cluster_cell$Total_SAMHD1pos / freq_cluster_cell$Total)),digits=3)

#total immune
freq_cluster_pooled<-as.data.frame(table(data_immune$cluster))
freq_marker_pooled<-as.data.frame(table(data_marker$cluster))
freq_cluster_pooled$SAMHD1pos <- freq_marker_pooled$Freq
names(freq_cluster_pooled)<-c("cluster","Total","Total_SAMHD1pos")
freq_cluster_pooled$percentSAMHD1pos<-as.numeric(format((freq_cluster_pooled$Total_SAMHD1pos / freq_cluster_pooled$Total)),digits=3)

#add total immune data to the 
freq_cluster_pooled<-add_column(freq_cluster_pooled, d = 'immune', .after = "cluster")
colnames(freq_cluster_pooled)<-c("cluster","cell_lineage","Total","Total_SAMHD1pos","percentSAMHD1pos")
freq_cluster_cell<-rbind(freq_cluster_pooled,freq_cluster_cell)

#turn to heatmap format
cell_cluster_SAMHD1_hmap<-dcast(freq_cluster_cell, cluster ~ factor(cell_lineage, levels = unique(cell_lineage)), 
                               value.var = "percentSAMHD1pos")
cell_cluster_SAMHD1_hmap<-as.matrix(cell_cluster_SAMHD1_hmap[,-1])
cell_cluster_SAMHD1_hmap[is.na(cell_cluster_SAMHD1_hmap)] <- 0
rownames(cell_cluster_SAMHD1_hmap)<-c(1,2,3,4,5,6,7,8,9,10)

heatmap.2(t(cell_cluster_SAMHD1_hmap), 
          Colv = F, Rowv = F,
          dendrogram = 'none',
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Grays')),
          sepcolor="grey35",
          colsep=0:ncol(cell_cluster_SAMHD1_hmap),
          rowsep=0:nrow(cell_cluster_SAMHD1_hmap),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '',
          ColSideColors = colors,
          cexRow = 1, cexCol = 2, margins = c(8,14))
