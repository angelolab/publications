# TBMIBIconcatenatecsv.R
# Author: Erin McCaffrey
# Date created: 190116
# Overview: This script reads in the csv for raw, cell size normalized, and asinh-transformed 
# cell size normalized expression data for all points in the study, appends a column with the 
# sample ID number and tissue type, then concatenates dataframes based on sample groupings of 
# data type. Exports as csv files for later use.


##..Set wd and load packages/libraries..##

require(gtools)
library(dplyr)

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")


##..Assign tissue type to sample IDs..##

LN<-c(3,4) 
spleen<-c(5,6)
gran_D1MT<-c(7,8,9,13,14,15,16,17,18,25,26,27,29,30,31,39,40,41)
gran_con<-c(1,2,10,11,12,19,20,21,22,23,24,28,32,33,34,35,36,37,38) 

##..Read in all raw, size-normed, and size-normed+asinh-transformed csv..##

#raw
files_raw<-list.files(pattern='dataFCS_.*\\.csv') 
files_raw<-mixedsort(files_raw) 
data_raw<-lapply(files_raw,read.csv)

#cs normed
files_sizenorm<-list.files(pattern='dataScaleSizeFCS_.*\\.csv')
files_sizenorm<-mixedsort(files_sizenorm) 
data_sizenorm<-lapply(files_sizenorm,read.csv)

##..Create concetenated matrices..##

#create dataframes of all with sample ID
all_raw<-bind_rows(data_raw, .id="SampleID") #raw

all_sizenorm<-bind_rows(data_sizenorm, .id="SampleID") #normalized by cell size


#add tissue type as label

all_raw<-all_raw %>% mutate(Tissue=case_when(all_raw$SampleID %in% LN ~ "tonsil",
                                          all_raw$SampleID %in% spleen ~ "spleen",
                                          all_raw$SampleID %in% gran_con ~ "gran_control",
                                          all_raw$SampleID %in% gran_D1MT ~ "gran_D1MT"))

all_sizenorm$Tissue<-all_raw$Tissue

#add internal Sample ID to the mycobacterial granulomas

all_raw<-all_raw %>% mutate(PatientID=case_when(all_raw$SampleID %in% c(3,4,22,23,24,33,34,35) ~ '15A092',
                                                all_raw$SampleID %in% c(19,20,21,10,11,12) ~ '15A093',
                                                all_raw$SampleID %in% c(29,30,31,39,40,41) ~ '15A209',
                                                all_raw$SampleID %in% c(16,17,18,25,26,27) ~ '16A101',
                                                all_raw$SampleID %in% c(7,8,9,13,14,16,5,6) ~ '16A105',
                                                all_raw$SampleID %in% c(1,2,36,37,38) ~ '17A235',
                                                all_raw$SampleID %in% c(28,32) ~ '17A236'))

all_sizenorm$PatientID<-all_raw$PatientID

##..Save concatenated dataframes..##
write.csv(all_raw, file="allsamples_data_annotated.csv",row.names = FALSE)
write.csv(all_sizenorm, file="allsamples_dataCS_annotated.csv",row.names = FALSE)
