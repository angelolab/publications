# TBMIBIdropHIVSamples.R
# Author: Erin McCaffrey
# Date created: 210621
# Overview: Script reads in the ME annotated data and the single cell data and removes the rows for 
# SampleIDs 64, 65, 85, 13, 35, and 36 which belong to patients 1, 5, and 6 which have confirmed HIV-coninfection.

library(dplyr)

##..Read in annotated data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
cell_data<-read.csv("allgran_cellpheno_CS-asinh-norm_combo_revised.csv")
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
me_data<-read.csv("all_TB_topic_annotation_combo.csv")

##..Drop HIV+ samples..##

drop<-c(64, 65, 85, 13, 35, 36)
cell_data_noHIV<-droplevels(cell_data[!cell_data$SampleID %in% drop,])

me_data_noHIV<-droplevels(me_data[!me_data$SampleID %in% drop,])

##..Export new csv..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
write.csv(me_data_noHIV, "all_TB_topic_annotation_noHIV.csv", row.names = F)

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
write.csv(cell_data_noHIV, "allgran_cellpheno_noHIV.csv", row.names = F)
