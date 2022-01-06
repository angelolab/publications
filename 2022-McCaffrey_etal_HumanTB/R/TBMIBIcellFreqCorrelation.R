# TBMIBIcellFreqCorrelation.R
# Author: Erin McCaffrey 
# Date created: 191103
# Overview: This script reads in the frequecy data for all TB granulomas and then runs a pearson correlation to produce
# a correlation matrix. 

library(factoextra)
library(data.table)
library(ggplot2)
library(PerformanceAnalytics)
library(psych)
library(Hmisc)
library(corrplot)
library(gplots)
library(colorspace)

#Function for flattening matrix from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# import data
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

# cell frequencies
imm_freqs_combo<-read.csv("immune_cell_freqs_noHIV.csv")

#reshape to a matrix where rows are samples and columns are counts of cell types
immfreq_data <- reshape2::dcast(imm_freqs_combo, SampleID ~ cell_type, value.var = "frequency")
immfreq_data[is.na(immfreq_data)] <- 0

#drop sarcoid
sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
immfreq_data <- immfreq_data[!immfreq_data$SampleID %in% sarcoid,]

# plot as a heatmap with cell types heirarchically clustered 

hm<-as.matrix(immfreq_data[,-1])
row.names(hm)<-immfreq_data$SampleID
              
heatmap.2(hm,
          scale = "none",
          Colv = T, Rowv = F,
          dendrogram = "column",
          trace = "none",
          col = colorRampPalette((brewer.pal(9,"Blues")))(100),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="black",
          colsep=0:ncol(hm),
          rowsep=0:nrow(hm),
          cexCol=0.4, cexRow =0.4)

# run correlation
immfreq_data.t <- as.data.frame(t(immfreq_data))
colnames(immfreq_data.t) <- immfreq_data$SampleID

corr<-corr.test(as.matrix(immfreq_data.t[-1,]), y = NULL, use = "pairwise",method="pearson",adjust="fdr")
corr_flat<-flattenCorrMatrix(corr$r,corr$p)
corr_flat

corrplot(corr$r, type = "upper", order = "hclust", tl.col = "black")

# plot as heatmaop with similar order to produce dendogram for freqeucny chart

heatmap.2(x = corr$r,
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          trace = "none",
          dendrogram = "col",
          distfun = function(x) as.dist(1-x),
          hclustfun = function(x) hclust(x, method="complete"),
          density.info = 'none',
          key.title = '',
          cexRow=0.6,
          cexCol=0.6,
          symm = TRUE)

# get ideal number of clusters with wss
data<-corr$r
dist_mat<-as.dist(1 - data)
hc <- hclust(dist_mat, method = 'complete')
plot(hc)
rect.hclust(hc, k = 5)


library("GMD")
cluster_stats<-css.hclust(dist_mat, hclust.obj=hc, k=10)
plot(attr(cluster_stats,"meta")$hclust.obj)

# plot explained variance

ggplot(cluster_stats, aes(x=k, y=ev)) +
  geom_line() +
  geom_point()

# get randomized data and randomized correlations

immfreq_rand_data.t<-immfreq_data.t[-1,]
immfreq_rand_data.t<- as.data.frame(lapply(immfreq_rand_data.t, sample))
colnames(immfreq_rand_data.t)<-colnames(immfreq_data.t)
row.names(immfreq_rand_data.t)<-row.names(immfreq_data.t[-1,])

corr_rand<-corr.test(as.matrix(immfreq_rand_data.t), y = NULL, use = "pairwise",method="pearson",adjust="fdr")

data_rand<-corr_rand$r
dist_mat_rand<-as.dist(1 - data_rand)
hc_rand<- hclust(dist_mat_rand, method = 'complete')
plot(hc_rand)

cluster_stats_rand<-css.hclust(dist_mat_rand, hclust.obj=hc_rand, k=10)

ggplot(cluster_stats_rand, aes(x=k, y=ev)) +
  geom_line() +
  geom_point()

# combine random and real

cluster_stats$data <- 'real'
cluster_stats_rand$data<-'randomized'
cluster_stats_comb<-rbind(cluster_stats,cluster_stats_rand)


ggplot(cluster_stats_comb, aes(x= k, y = ev, col = data, group=data))+
  geom_point()+
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,10,1)) 




