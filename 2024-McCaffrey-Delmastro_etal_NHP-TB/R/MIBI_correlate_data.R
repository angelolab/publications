# MIBI_correlate_data.R
# Created by: Erin McCaffrey 
# Date created: March 2, 2023
#
# Overview: Script takes a data sheet of feature (ie. cell frequencies) and 
# runs correlation between features and CFU in bulk or broken down by some 
# metadata (ie. time)

library(Hmisc)
library("PerformanceAnalytics")
library(corrplot)
library(ggplot2)
library(tidyr)
library(reshape2)
library(psych)
library(tibble)
library(ggstatsplot)
library(dplyr)
library(EnhancedVolcano)


# We will use this function from: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

##..Step 1: Read in data..##

data<-read.csv('cell_stats_all_samples_meta_data.csv')
data<-droplevels(data[data$category == 'pheno_of_total',])
data<-tibble::rowid_to_column(data, "ID")
color_key <- read.csv("./keys/cell_color_key.csv")

##..Step 2: Convert frequency data to format where there is 1 row per ID and columns are the freqs or counts..##

freq_data_bulk<-reshape2::dcast(data, sample + gran_CFU + necrosis_viable_ratio ~ variable, value.var = "freq_of_total", fun.aggregate = sum)
count_data_bulk<-reshape2::dcast(data, sample + gran_CFU + necrosis_viable_ratio ~ variable, value.var = "cell_density", fun.aggregate = sum)

##..Step 3: Convert NaN to 0..##

freq_data_bulk[is.na(freq_data_bulk)] <- 0
count_data_bulk[is.na(count_data_bulk)] <- 0

##..Step 4: Add in macrophage ratio..##

count_data_bulk$mac_ratio<-count_data_bulk$`CD11c+_Mac`/count_data_bulk$`CD14+CD11c+_Mac`
count_data_bulk$log_CFU<-log10(1+count_data_bulk$gran_CFU)
count_data_bulk$t_ratio<-count_data_bulk$`CD4+Tcell`/count_data_bulk$`CD8+Tcell`
freq_data_bulk$log_CFU<-log10(1+freq_data_bulk$gran_CFU)

##..Step 5: Remove unassigned cells due to inability to interpret..##

count_data_bulk<-count_data_bulk[ , -which(names(count_data_bulk) %in% c("Unassigned"))]
freq_data_bulk<-freq_data_bulk[ , -which(names(freq_data_bulk) %in% c("Unassigned"))]

##..Step 6: Perform correlation analysis of everything against everything..##

# option to subset
corrdata<-freq_data_bulk
corrmat<-rcorr(as.matrix(corrdata[,2:23]), type = 'spearman')
corr_results<-flattenCorrMatrix(corrmat$r, corrmat$P)

# fun bubbles
corrmatrix<- cor(as.matrix(corrdata[,2:23]), method = 'spearman')
# all combos
corrplot(corrmatrix, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
# cross out statistically insignificant
corrplot(corrmatrix, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, p.mat = corrmat$P, sig.level = 0.01, insig = "blank")


##..Step 7: Plot individual examples..##

ggscatterstats(
  data = corrdata, 
  x    = log_CFU,
  y    = "VIM+Stroma",
  type = "spearman")

ggplot(freq_data_bulk, aes(x=log_CFU, y = `CD4+Tcell`)) + 
  geom_point()

individual_corr <-corr.test(corrdata$gran_CFU, corrdata$`CD11c+_Mac`, method = 'spearman')
individual_corr

##..Step 8: Generate volcano of density relationships with CFU..##

CFU_corr_data<-corr_results[corr_results$row =='log_CFU' | corr_results$column == 'log_CFU',]
CFU_corr_data<-droplevels(CFU_corr_data[CFU_corr_data$row !='necrosis_viable_ratio' & 
                                          CFU_corr_data$row != 'mac_ratio' &
                                          CFU_corr_data$row != 'log_CFU' &
                                          CFU_corr_data$row != 'gran_CFU' &
                                          CFU_corr_data$column!= 't_ratio',])

#adjust p-values
CFU_corr_data$adj.p <- p.adjust(CFU_corr_data$p, method="fdr")


EnhancedVolcano(CFU_corr_data,
                lab = CFU_corr_data$row,
                title = 'Cell density versus CFU',
                x = 'cor',
                y = 'adj.p',
                pCutoff = 0.3,
                FCcutoff = 0,
                pointSize = 2.0,
                labSize = 3.0,
                legendLabels=c('NS','Rho > 0.5','Adj. p < 0.05',
                               'Adj. p < 0.05 & Rho > 0.5'),
                legendPosition = 'none',
                legendLabSize = 10,
                legendIconSize = 5.0,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                hline = c(0.05),
                ylab = '-log10(adj p)',
                xlab = 'Spearman Rho',
                xlim = c(-1,1),
                ylim = c(0, max(-log10(CFU_corr_data$adj.p), na.rm=T) + 0.3))

##..Step 9: Plot bi-axial scatters for all relationships..##

# melt data
corrdata.m <- melt(corrdata, id.vars = c('sample','log_CFU','gran_CFU'))

# plot
ggplot(corrdata.m, aes(x=log_CFU, y=as.numeric(value))) + 
  geom_point() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = "None") +
  facet_wrap(.~variable, scale='free_y', ncol = 4)


                      