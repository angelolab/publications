#TBMIBIplotGeneESProgression.R
#Date created: 03/16/2020
#Author: Erin McCaffrey
#This script reads in the mean effect size relative to each stage of progression in TB infection
# and produces a heatmap clustered by column. 

library(gplots)
library(colorspace)

##...Load in data..##


ES_data<-read.csv("data/ESprogression_profiles_allgenes.csv")


##...Convert data to labeled matrix...##

hmap<-as.matrix(ES_data[-1,-1])
rownames(hmap)<-ES_data[-1,]$X

##..Visualize as heatmap..##

heatmap.2(hmap, 
          Colv = T, Rowv = F,
          hclustfun = hclust,
          scale = "none",
          dendrogram = c("column"),
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          sepcolor="grey35",
          colsep=0:ncol(hmap),
          rowsep=0:nrow(hmap),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '',
          cexRow = 2, cexCol = 0.5, margins = c(8,14))

