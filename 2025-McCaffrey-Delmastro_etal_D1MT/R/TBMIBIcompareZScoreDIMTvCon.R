# TBMIBIcompareZScoreDIMTvCon
# Author: Erin McCaffrey 
# Date created: 201130

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(ggpmisc)

##..Import data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/")
data<-read.delim('spatialAll_dataset.txt', header = TRUE, sep = "\t", dec = ".")
rownames(data)<-data$Observations
data<-data[,-1]

##..Transpose..##

data.t<-as.data.frame(t(as.matrix(data)))
rownames(data.t) <- substring(rownames(data.t),2)
points<-rownames(data.t)
data.t$SampleID<-as.numeric(points)

##..Add group information..##

gran_D1MT<-c(8,9,13,14,15,16,18,25,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38) 
data.t$Group<-'gran_D1MT'
data.t[data.t$SampleID %in% gran_con,]$Group<-'gran_con'

##..Split into two dataframes..##

data_D1MT<-droplevels(data.t[data.t$Group == 'gran_D1MT',])
data_con<-droplevels(data.t[!data.t$Group == 'gran_D1MT',])

##..Run pairwise t test between all variables and store output..##

# T-TEST FCT
tfct <- function(v1, v2){
  t.test(v1, v2) 
}

# RUN T-TESTS BY COL, SAVE RESULTS TO LIST
p_values = as.data.frame(matrix(, nrow = 1, ncol = 900))
for(i in 1:900) {
  test <- t.test(data_D1MT[,i],data_con[,i])
  p_values[1, i] <- test$p.value
} 
colnames(p_values)<-colnames(data.t[,1:900])

p_values_adj<-p.adjust(p_values[1,], method = "BH")

##..Get significant interactions..##

sig_results<-p_values %>% select_if(~any(. <0.05))






