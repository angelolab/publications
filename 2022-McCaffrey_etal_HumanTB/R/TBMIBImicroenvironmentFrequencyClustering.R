# TBMIBImicroenvironmentFrequencyClustering.R
# Author: Erin McCaffrey 
# Date created: 200520
# Overview: This script creates a matrix of all patients and the frequency of each ME in that patient. 
# It next runs a clustering analysis of that data based on correlation and euclidean distance with 
# heirarchical clustering.

library(reshape2)
library(corrplot)
library(psych)
library(gplots)
library(ggplot2)
library(colorspace)
library("GMD")

##...Load in data..##

ME_data<-read.csv('data/allTB-ME_annotated.csv')

##...Enumerate frequency of each ME per patient..##

ME_freqs<-as.data.frame(table(ME_data$SampleID, ME_data$maxME))
colnames(ME_freqs)<-c('SampleID','ME','count')
patient_totals<-as.vector(table(ME_data$SampleID))
ME_freqs$total<-rep(patient_totals, n=8)
ME_freqs$freq<-ME_freqs$count / ME_freqs$total

##..Cast to a matrix that is all patients (columns) and all MEs (rows)..##

ME_mat<-reshape2::dcast(ME_freqs, ME ~ SampleID, value.var = "freq")
MEmat.m<-as.matrix(ME_mat[,-1])
row.names(MEmat.m)<-ME_mat$ME

##..Plot heatmap with hierarchical clustering along columns based on euclidean..##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

heatmap.2(MEmat.m, 
          Colv = T, Rowv = F,
          hclustfun = hclust,
          scale = "none",
          dendrogram = "column",
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(MEmat.m),
          rowsep=0:nrow(MEmat.m),
          RowSideColors = colors,
          breaks=seq(0, 1, length.out=101),
          cexRow = 2, cexCol = 1, margins = c(8,14))

##..Cluster analysis of hclust with euclidean..##

# get ideal number of clusters with wss
data<-MEmat.m
dist_mat<-dist(t(data), method = 'euclidean')
hc<-hclust(dist(t(data), method = 'euclidean'), method = 'complete')
plot(hc)
rect.hclust(hc, k = 5)

cluster_stats<-css.hclust(dist_mat, hclust.obj=hc, k=10)
plot(attr(cluster_stats,"meta")$hclust.obj)

# plot explained variance

ggplot(cluster_stats, aes(x=k, y=ev)) +
  geom_line() +
  geom_point()

##..Plot heatmap with hierarchical clustering along columns based on correlation..##

heatmap.2(MEmat.m, 
          Colv = T, Rowv = F,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="complete"),
          scale = "none",
          dendrogram = "column",
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(MEmat.m),
          rowsep=0:nrow(MEmat.m),
          RowSideColors = colors,
          breaks=seq(0, 1, length.out=101),
          cexRow = 2, cexCol = 1, margins = c(8,14))

##..Cluster analysis of Pearson correlation-based clustering..##

data<-MEmat.m
corr<-corr.test(data, y = NULL, use = "pairwise",method="pearson",adjust="fdr")

# get ideal number of clusters with wss
data_corr<-corr$r
dist_mat<-as.dist(1 - data_corr)
hc <- hclust(dist_mat, method = 'complete')
plot(hc)
rect.hclust(hc, k = 5)

cluster_stats<-css.hclust(dist_mat, hclust.obj=hc, k=10)
plot(attr(cluster_stats,"meta")$hclust.obj)

# plot explained variance

ggplot(cluster_stats, aes(x=k, y=ev)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,10,1)) 

