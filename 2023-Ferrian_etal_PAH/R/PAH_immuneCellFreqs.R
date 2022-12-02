# PAH_immuneCellFreqs.R
# Author: Erin McCaffrey 
# Date created: 200307
# Overview: This script reads in the csv for the PAH dataset, determines the percent of 
# lineage of totl and the percentage of immune cell clusters of total immune. It produces 
# summary stats in PAH and healthy of bulk trends, pairwise PAH v Healthy, breaks down PAH type, 
# and produces stacked bars of frequency data.


library(dplyr)
library(viridis)
library(ggplot2)
library(forcats)
library(reshape2)
library(tidyr)
library(ggpubr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/PAH manuscript/Datasets")
data<-read.csv("ImmuneMetaclust_K=30.csv")

##..Create a table of cell frequencies..##

imm_freqs<-as.data.frame(table(data$Point_num,data$cell_lineage))
colnames(imm_freqs)<-c('PointNum','Cell_Type','Count')
totals<-as.numeric(table(data$Point_num))
imm_freqs$total<-rep(totals,8)
imm_freqs$freq<-as.numeric(imm_freqs$Count/imm_freqs$total)

##..Create a color key..##

# create numerical code

clusterID<-as.character(data$cell_lineage)

clusterID<-replace(clusterID,clusterID =='Bcell', 1)
clusterID<-replace(clusterID,clusterID =='DC', 2)
clusterID<-replace(clusterID,clusterID =='Macro', 3)
clusterID<-replace(clusterID,clusterID =='Mono', 4)
clusterID<-replace(clusterID,clusterID =='Neutro', 5)
clusterID<-replace(clusterID,clusterID =='NK', 6)
clusterID<-replace(clusterID,clusterID =='Tc', 7) 
clusterID<-replace(clusterID,clusterID =='Th', 8)

data$clusterID<-clusterID

# extract out cell type, code pairs, assign hexcode

cell_color_key<-unique(data[,c('cell_lineage','clusterID')])
cell_color_key$codes<-c('#F99157','#C594C5','#FAC863','#99C794',
         '#5FB3B3','#6699CC','#EC5f67','#AB7967')
order<-c('Tc','Th','Mono','DC','Neutro','Macro','NK','Bcell')

##..Order by decreasing median or plotting..##

imm_freqs$Cell_Type <- factor(imm_freqs$Cell_Type, levels=order)
imm_freqs<-imm_freqs[order(imm_freqs$Cell_Type),]

cell_color_key$cell_lineage<-factor(cell_color_key$cell_lineage, levels = order)
cell_color_key<-cell_color_key[order(cell_color_key$cell_lineage),]
color<-cell_color_key$codes

##..Bulk dataset stats..##

bulk_box<-ggplot(imm_freqs, aes(x=Cell_Type, y=freq, fill=Cell_Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Cell Types-All Samples") +
  guides(fill=guide_legend(title="Cell Type"))
bulk_box

#healthy

hlt<-unique(data[data$Subgroup == 'Healthy Controls',]$Point_num)
hlt_box<-ggplot(imm_freqs[imm_freqs$PointNum %in% hlt,], aes(x=Cell_Type, y=freq, fill=Cell_Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Cell Types-Healthy") +
  guides(fill=guide_legend(title="Cell Type"))
hlt_box

#hpah

hpah<-unique(data[data$Subgroup == 'HPAH',]$Point_num)
hpah_box<-ggplot(imm_freqs[imm_freqs$PointNum %in% hpah,], aes(x=Cell_Type, y=freq, fill=Cell_Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Cell Types-HPAH") +
  guides(fill=guide_legend(title="Cell Type"))
hpah_box


#ipah

ipah<-unique(data[data$Subgroup == 'IPAH',]$Point_num)
ipah_box<-ggplot(imm_freqs[imm_freqs$PointNum %in% ipah,], aes(x=Cell_Type, y=freq, fill=Cell_Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Cell Types-IPAH") +
  guides(fill=guide_legend(title="Cell Type"))
ipah_box

# total immune cells all groups

# add the hlt, ipah, and hpah 
imm_freqs$Group<-'hlt'
imm_freqs[imm_freqs$PointNum %in% hpah,]$Group<-'hpah'
imm_freqs[imm_freqs$PointNum %in% ipah,]$Group<-'ipah'

imm_freqs$PAH<-'pah'
imm_freqs[imm_freqs$Group =='hlt',]$PAH <-'hlt'

my_comparisons <- list(c("hlt","ipah"),
                       c("hlt","hpah"),
                       c("ipah","hpah"))

totals<-ggplot(unique(imm_freqs[,c('PointNum','Group','total')]), aes(x=Group, y=total, fill=Group)) + 
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  labs(x="Group") +
  labs(y="Total Immune Cells") +
  ggtitle("Total Immune Cell Counts") +
  guides(fill=guide_legend(title="Group"))
totals

##..Stacked bars all points..##

imm_freqs$Cell_Type <- factor(imm_freqs$Cell_Type, levels=rev(order))
imm_freqs<-imm_freqs[order(imm_freqs$Cell_Type),]

imm_bar_hlt<-ggplot(imm_freqs[imm_freqs$Group =='hlt',], aes(x=PointNum, y=freq, fill=Cell_Type)) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
imm_bar_hlt

imm_bar_hpah<-ggplot(imm_freqs[imm_freqs$Group =='hpah',], aes(x=PointNum, y=freq, fill=Cell_Type)) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
imm_bar_hpah

imm_bar_ipah<-ggplot(imm_freqs[imm_freqs$Group =='ipah',], aes(x=PointNum, y=freq, fill=Cell_Type)) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Point Number") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
imm_bar_ipah

##..Compare between all groups, all cell types..##

# individual plots hlt v ipah v hpah
compare_freqs<-ggplot(data = imm_freqs, aes(x = Group, y = freq, fill = Group)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  labs(x = 'Group') + 
  labs(y = 'Frequency') +
  facet_wrap(~Cell_Type, scale='free_y')
compare_freqs

compare_freqs_pahvhlt<-ggplot(data = imm_freqs, aes(x = PAH, y = freq, fill = PAH)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons=list(c('hlt','pah')), label = "p.signif", method= "wilcox.test") +
  labs(x = 'Group') + 
  labs(y = 'Frequency') +
  facet_wrap(~Cell_Type, scale='free_y')
compare_freqs_pahvhlt

# grouped bar chart

imm_freqs$Cell_Type <- factor(imm_freqs$Cell_Type, levels=order)
imm_freqs<-imm_freqs[order(imm_freqs$Cell_Type),]

grouped_bar<-ggplot(data = imm_freqs, aes(x = Cell_Type, y = freq, fill=PAH)) +
  geom_boxplot() +
  stat_compare_means(aes(group = PAH), method= "wilcox.test", label = "p.signif") +
  theme_bw() +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 
grouped_bar

grouped_bar_ipahvhpah<-ggplot(data = imm_freqs[!imm_freqs$Group=='hlt',], aes(x = Cell_Type, y = freq, fill=Group)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Group), method= "wilcox.test", label = "p.signif") +
  theme_bw() +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 
grouped_bar_ipahvhpah

