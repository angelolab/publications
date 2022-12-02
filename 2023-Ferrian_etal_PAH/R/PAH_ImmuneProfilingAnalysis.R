# PAH_MIBIFunctionalProfilingAnalysis
# Author: Selena Ferrian
# Date created: 04/21/2020
# Overview: Script imports the median freq of positive cell per group and display spider/radar graphs.


##########################################
##..Install packages and open libraries..##
###########################################

# install.packages("BiocManager")
# BiocManager::install("BiocUpgrade")  # upgrading Bioconductor
# BiocManager::install("ggplot2")      # for advanced data plotting
# BiocManager::install("gplots")       # for heatmaps
# BiocManager::install("RColorBrewer") # additional color schemes
# BiocManager::install("reshape2")     # reorganizing data
# BiocManager::install("plyr")         # for data handling
# install.packages("viridis")          # for badass color palettes
# install.packages("fmsb")             # for radar/spider plots


##..Open necessary libraries..##

          
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plyr)
library(viridis)
library(fmsb)


##..Subgroup analysis based on Median..##

options(stringsAsFactors = F)
rm(list =ls())
library(fmsb)
setwd("~/Desktop/FinalData/PAHanalyses/single_cell_data")
data<-read.csv("Figure3K_allpatients_freq-of-positivecells-persubset.csv")

###LINEAGE ANALYSIS###

#Subset Immune
dataImm <- data[c(36:37),]
head(dataImm)
dim(dataImm)
dataImm <- droplevels(dataImm[,3:13])
head(dataImm)
dataImm <- rbind(rep(1,11) , rep(0,11) , dataImm)
radarchart(dataImm, pcol=c("3","9"), pfcol=c(rgb(0.2,0.5,0.5,0.4), rgb(0.1,0.2,0.5,0.4)), cglcol = "grey80", plwd= 4, vlcex=1.5, plty=1, seg = 5)

#Subset Endo
dataEndo <- data[c(38:39),]
head(dataEndo)
dim(dataEndo)
dataEndo <- droplevels(dataEndo[,3:13])
head(dataEndo)
dataEndo <- rbind(rep(1,11) , rep(0,11) , dataEndo)
radarchart(dataEndo, pcol=c("3","9"), pfcol=c(rgb(0.2,0.5,0.5,0.4), rgb(0.1,0.2,0.5,0.4)), cglcol = "grey80", plwd= 4, vlcex=1.5, plty=1, seg = 5)

#Subset Epi
dataEpi <- data[c(40:41),]
head(dataEpi)
dim(dataEpi)
dataEpi <- droplevels(dataEpi[,3:13])
head(dataEpi)
dataEpi <- rbind(rep(1,11) , rep(0,11) , dataEpi)
radarchart(dataEpi, pcol=c("3","9"), pfcol=c(rgb(0.2,0.5,0.5,0.4), rgb(0.1,0.2,0.5,0.4)), cglcol = "grey80", plwd= 4, vlcex=1.5, plty=1, seg = 5)

#Subset Fibro
dataFi <- data[c(42:43),]
head(dataFi)
dim(dataFi)
dataFi <- droplevels(dataFi[,3:13])
head(dataFi)
dataFi <- rbind(rep(1,11) , rep(0,11) , dataFi)
radarchart(dataFi, pcol=c("3","9"), pfcol=c(rgb(0.2,0.5,0.5,0.4), rgb(0.1,0.2,0.5,0.4)), cglcol = "grey80", plwd= 4, vlcex=1.5, plty=1, seg = 5)

#Subset Mese
dataMese <- data[c(44:45),]
head(dataMese)
dim(dataMese)
dataMese <- droplevels(dataMese[,3:13])
head(dataMese)
dataMese <- rbind(rep(1,11) , rep(0,11) , dataMese)
radarchart(dataMese, pcol=c("3","9"), pfcol=c(rgb(0.2,0.5,0.5,0.4), rgb(0.1,0.2,0.5,0.4)), cglcol = "grey80", plwd= 4, vlcex=1.5, plty=1, seg = 5)



###IMMUNE SUBSETS ANALYSIS###


#Subset DCs
dataDC <- data[c(2:3),]
head(dataDC)
dim(dataDC)
dataDC <- droplevels(dataDC[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataDC)
dataDC <- rbind(rep(1,11) , rep(0,11) , dataDC)
radarchart(dataDC, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset Macro
dataMacro <- data[c(5:6),]
head(dataMacro)
dim(dataMacro)
dataMacro <- droplevels(dataMacro[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataMacro)
dataMacro <- rbind(rep(1,11) , rep(0,11) , dataMacro)
radarchart(dataMacro, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset Mono
dataMono <- data[c(8:9),]
head(dataMono)
dim(dataMono)
dataMono <- droplevels(dataMono[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataMono)
dataMono <- rbind(rep(1,11) , rep(0,11) , dataMono)
radarchart(dataMono, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset Neutro
dataNeutro <- data[c(11:12),]
head(dataNeutro)
dim(dataNeutro)
dataNeutro <- droplevels(dataNeutro[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataNeutro)
dataNeutro <- rbind(rep(1,11) , rep(0,11) , dataNeutro)
radarchart(dataNeutro, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset NK
dataNK <- data[c(14:15),]
head(dataNK)
dim(dataNK)
dataNK<- droplevels(dataNK[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataNK)
dataNK <- rbind(rep(1,11) , rep(0,11) , dataNK)
radarchart(dataNK, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset Bcell
dataB <- data[c(22:23),]
head(dataB)
dim(dataB)
dataB <- droplevels(dataB[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataB)
dataB <- rbind(rep(1,11) , rep(0,11) , dataB)
radarchart(dataB, pcol=c("green","purple","blue"), cglcol = "grey80", seg = 10)
radarchart(dataB, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset Tc
dataTc <- data[c(17:18),]
head(dataTc)
dim(dataTc)
dataTc <- droplevels(dataTc[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataTc)
dataTc <- rbind(rep(1,11) , rep(0,11) , dataTc)
radarchart(dataTc, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

#Subset Th
dataTh <- data[c(20:21),]
head(dataTh)
dim(dataTh)
dataTh <- droplevels(dataTh[,2:12])
#data <- data[,-c(3:22,31,32,33,36,38:41)] #leave PID and lineage
head(dataTh)
dataTh <- rbind(rep(1,11) , rep(0,11) , dataTh)
radarchart(dataTh, pcol=c("6","4"), pfcol=c(rgb(0.8,0.2,0.5,0.4) , rgb(0.6,0.5,0.9,0.4)), cglcol = "grey80", axislabcol= "grey", caxislabels=seq(0,1,0.2), plwd= 4, vlcex=1.5, plty=1, seg = 5) 

