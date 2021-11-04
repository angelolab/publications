###########################################
##..Install packages and open libraries..##
###########################################

source("~/Risom2022_R_scripts/Risom2022_SC_FUNCTIONS")

##..Open necessary libraries..##

library(flowCore)           
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(FlowSOM)
library(plyr)
library(dplyr)
library(viridis)
library(data.table) # to create data table objects
library(wesanderson) # color pallate
library(scales) # for non-scientific notation
library(umap) # umap library
library(forcats) # reverse ordering of factors
library(gdata)
################################# Sam NORMALIZE #######################################################

#setwd("/Volumes/MIBI Files/remote/dataPerCell")
setwd("~/RDATA/DCIS/200428_FRANKENSTEIN5/")

##### Made a file where I added three columns as ERstatus ARstatus HER2status all filled with zeros, will populate with "1" if positive
data<-read.csv("200706_CellTable_Fullinfo.csv")
# tumor <- c("tumor")
# tumordata<-droplevels(tumordata[tumordata$celllineage %in% tumor,])

info.dir = "~/RDATA/DCIS/200428_FRANKENSTEIN5/metadata/"
info.csv <- read.csv(file = paste0(info.dir, "200707_COHORT_METADATA.csv"))


data$ERstatus<-0
data$HER2status<-0
data$ARstatus<-0
data$Ki67status<-0
data$pS6status<-0
data$GLUT1status<-0
data$HIF1astatus<-0
data$COX2status<-0
data$CD36status<-0
data$CD44status<-0
data$PD1status<-0
data$PDL1status<-0
data$IDO1status<-0
data$GZMBstatus<-0
data$MMP9status<-0
# data$COLIstatus<-0
# data$MPOstatus<-0
# data$VIMstatus<-0
# data$SMAstatus<-0
# data$CK5status<-0
data$ECADstatus<-0
# data$PanKRTstatus<-0
# data$FOXP3status<-0
# data$CD56status<-0


HER2pos_row_idx<-which(data$HER2 > 0.2) #use positivity threshold of 0.5
HER2_rows<-rownames(data[HER2pos_row_idx,])
# levels(data$phenotype)<-c(levels(concatenatedcohort2$phenotype), "NEUT1" ) 
data[HER2_rows,]$HER2status<-"1"

HER2amp_row_idx<-which(data$HER2 > 0.65) #use positivity threshold of 0.5
HER2amp_rows<-rownames(data[HER2amp_row_idx,])
data[HER2amp_rows,]$HER2status<-"2"

ER_row_idx<-which(data$ER > 0.25) #use positivity threshold of 0.5
ER_rows<-rownames(data[ER_row_idx,])
data[ER_rows,]$ERstatus<-"1"

AR_row_idx<-which(data$AR > 0.25) #use positivity threshold of 0.5
AR_rows<-rownames(data[AR_row_idx,])
data[AR_rows,]$ARstatus<-"1"

Ki67_row_idx<-which(data$Ki67 > 0.25) #use positivity threshold of 0.5
Ki67_rows<-rownames(data[Ki67_row_idx,])
data[Ki67_rows,]$Ki67status<-"1"

GLUT1_row_idx<-which(data$GLUT1 > 0.25) #use positivity threshold of 0.5
GLUT1_rows<-rownames(data[GLUT1_row_idx,])
data[GLUT1_rows,]$GLUT1status<-"1"

pS6_row_idx<-which(data$pS6 > 0.25) #use positivity threshold of 0.5
pS6_rows<-rownames(data[pS6_row_idx,])
data[pS6_rows,]$pS6status<-"1"

HIF1a_row_idx<-which(data$HIF1a > 0.25) #use positivity threshold of 0.5
HIF1a_rows<-rownames(data[HIF1a_row_idx,])
data[HIF1a_rows,]$HIF1astatus<-"1"

COX2_row_idx<-which(data$COX2 > 0.25) #use positivity threshold of 0.5
COX2_rows<-rownames(data[COX2_row_idx,])
data[COX2_rows,]$COX2status<-"1"

CD36_row_idx<-which(data$CD36 > 0.25) #use positivity threshold of 0.5
CD36_rows<-rownames(data[CD36_row_idx,])
data[CD36_rows,]$CD36status<-"1"

CD44_row_idx<-which(data$CD44 > 0.25) #use positivity threshold of 0.5
CD44_rows<-rownames(data[CD44_row_idx,])
data[CD44_rows,]$CD44status<-"1"

PD1_row_idx<-which(data$PD1 > 0.25) #use positivity threshold of 0.5
PD1_rows<-rownames(data[PD1_row_idx,])
data[PD1_rows,]$PD1status<-"1"

PDL1_row_idx<-which(data$PDL1 > 0.25) #use positivity threshold of 0.5
PDL1_rows<-rownames(data[PDL1_row_idx,])
data[PDL1_rows,]$PDL1status<-"1"

IDO1_row_idx<-which(data$IDO1 > 0.25) #use positivity threshold of 0.5
IDO1_rows<-rownames(data[IDO1_row_idx,])
data[IDO1_rows,]$IDO1status<-"1"

MMP9_row_idx<-which(data$MMP9 > 0.25) #use positivity threshold of 0.5
MMP9_rows<-rownames(data[MMP9_row_idx,])
data[MMP9_rows,]$MMP9status<-"1"

GZMB_row_idx<-which(data$GZMB > 0.25) #use positivity threshold of 0.5
GZMB_rows<-rownames(data[GZMB_row_idx,])
data[GZMB_rows,]$GZMBstatus<-"1"

CK7_row_idx<-which(data$CK7 > 0.25) #use positivity threshold of 0.5
CK7_rows<-rownames(data[CK7_row_idx,])
data[CK7_rows,]$CK7status<-"1"

ECAD_row_idx<-which(data$ECAD > 0.25) #use positivity threshold of 0.5
ECAD_rows<-rownames(data[ECAD_row_idx,])
data[ECAD_rows,]$ECADstatus<-"1"




write.csv(data, file="200909_CellTable_Fullinfo_Status.csv",row.names = FALSE)



# LOAD DATA

data <- read.csv("200707_CellTable_Fullinfo_Status.csv")

#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

# Independent Marker Reports (later merged report will need this)

#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

#Select denominator


population <- c("tumor")
# cell_data<-data
cell_data <- droplevels(data[data$phenotype %in% population, ])




###########
### ER ####
###########

ER_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$ERstatus))
names(ER_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

cell_totals_tumor_ER<-aggregate(ER_Freqs_tumor$count, by=list(Category=ER_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(ER_Freqs_tumor$SampleID)) {
  frequencies <- ER_Freqs_tumor[ER_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_ER[cell_totals_tumor_ER$Category==i,2]
  ER_Freqs_tumor[ER_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
#..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_tumor_ER<-info.csv
# get list of pointnNums is frequency data
pointNumsER<-unique(cell_data$Point_Num)

# filter annotation data by PointNum
annotation_data_tumor_ER <- droplevels(annotation_data_tumor_ER[annotation_data_tumor_ER$PointNumber %in% pointNumsER, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_ER$PointNumber <- factor(annotation_data_tumor_ER$PointNumber, levels=pointNumsER)
annotation_data_tumor_ER<-annotation_data_tumor_ER[order(annotation_data_tumor_ER$CohortNumber),]

# cast the annotations to the frequency data
anno_data_tumor_ER<-rep(annotation_data_tumor_ER,1)

# meantumor<-tumordata %>% select(Point_Num, HER2, ER, AR, Ki67, ECAD)
# ER_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(ER), mean = mean(ER), se = sqrt(var(ER)/length(ER)))

annotated_Freqs_tumor_ER<-cbind(ER_Freqs_tumor, anno_data_tumor_ER)


# write.csv(annotated_Freqs_tumor_ER, file="200707_ER_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ AR ###################
########################################

##..Create a dataframe with the counts of each cell cluster across point number..##
AR_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$ARstatus))
names(AR_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

cell_totals_tumor_AR<-aggregate(AR_Freqs_tumor$count, by=list(Category=AR_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(AR_Freqs_tumor$SampleID)) {
  frequencies <- AR_Freqs_tumor[AR_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_AR[cell_totals_tumor_AR$Category==i,2]
  AR_Freqs_tumor[AR_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_tumor_AR<-info.csv
# get list of pointnNums is frequency data
pointNumsAR<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_AR <- droplevels(annotation_data_tumor_AR[annotation_data_tumor_AR$PointNumber %in% pointNumsAR, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_AR$PointNumber <- factor(annotation_data_tumor_AR$PointNumber, levels=pointNumsAR)
annotation_data_tumor_AR<-annotation_data_tumor_AR[order(annotation_data_tumor_AR$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_AR<-rep(annotation_data_tumor_AR,1)

# meantumor<-tumordata %>% select(Point_Num, HER2, ER, AR, Ki67)
# AR_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(AR), mean = mean(AR), se = sqrt(var(AR)/length(AR)))
annotated_Freqs_tumor_AR<-cbind(AR_Freqs_tumor, anno_data_tumor_AR)

# write.csv(annotated_Freqs_tumor_AR, file="200707_AR_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ HER2 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
HER2_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$HER2status))
names(HER2_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_HER2<-aggregate(HER2_Freqs_tumor$count, by=list(Category=HER2_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(HER2_Freqs_tumor$SampleID)) {
  frequencies <- HER2_Freqs_tumor[HER2_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_HER2[cell_totals_tumor_HER2$Category==i,2]
  HER2_Freqs_tumor[HER2_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_HER2<-info.csv
# get list of pointnNums is frequency data
pointNumsHER2<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_HER2 <- droplevels(annotation_data_tumor_HER2[annotation_data_tumor_HER2$PointNumber %in% pointNumsHER2, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_HER2$PointNumber <- factor(annotation_data_tumor_HER2$PointNumber, levels=pointNumsHER2)
annotation_data_tumor_HER2<-annotation_data_tumor_HER2[order(annotation_data_tumor_HER2$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_HER2<-rep(annotation_data_tumor_HER2,1)
# meantumor<-tumordata %>% select(Point_Num, HER2, ER, AR, Ki67)
# HER2_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(HER2), mean = mean(HER2), se = sqrt(var(HER2)/length(HER2)))
annotated_Freqs_tumor_HER2<-cbind(HER2_Freqs_tumor, anno_data_tumor_HER2)
# write.csv(annotated_Freqs_tumor_HER2, file="200707_HER2_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ Ki67 ###################
########################################

##..Create a dataframe with the counts of each cell cluster across point number..##
Ki67_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$Ki67status))
names(Ki67_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

cell_totals_tumor_Ki67<-aggregate(Ki67_Freqs_tumor$count, by=list(Category=Ki67_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Ki67_Freqs_tumor$SampleID)) {
  frequencies <- Ki67_Freqs_tumor[Ki67_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_Ki67[cell_totals_tumor_Ki67$Category==i,2]
  Ki67_Freqs_tumor[Ki67_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_tumor_Ki67<-info.csv
# get list of pointnNums is frequency data
pointNumsKi67<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_Ki67 <- droplevels(annotation_data_tumor_Ki67[annotation_data_tumor_Ki67$PointNumber %in% pointNumsKi67, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_Ki67$PointNumber <- factor(annotation_data_tumor_Ki67$PointNumber, levels=pointNumsKi67)
annotation_data_tumor_Ki67<-annotation_data_tumor_Ki67[order(annotation_data_tumor_Ki67$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_Ki67<-rep(annotation_data_tumor_Ki67,1)

# meantumor<-tumordata %>% select(Point_Num, HER2, ER, AR, Ki67, ECAD)
# Ki67_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(Ki67), mean = mean(Ki67), se = sqrt(var(Ki67)/length(Ki67)))
annotated_Freqs_tumor_Ki67<-cbind(Ki67_Freqs_tumor, anno_data_tumor_Ki67)

# write.csv(annotated_Freqs_tumor_Ki67, file="200402_Ki67_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ pS6 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
pS6_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$pS6status))
names(pS6_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_pS6<-aggregate(pS6_Freqs_tumor$count, by=list(Category=pS6_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(pS6_Freqs_tumor$SampleID)) {
  frequencies <- pS6_Freqs_tumor[pS6_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_pS6[cell_totals_tumor_pS6$Category==i,2]
  pS6_Freqs_tumor[pS6_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_pS6<-info.csv
# get list of pointnNums is frequency data
pointNumspS6<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_pS6 <- droplevels(annotation_data_tumor_pS6[annotation_data_tumor_pS6$PointNumber %in% pointNumspS6, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_pS6$PointNumber <- factor(annotation_data_tumor_pS6$PointNumber, levels=pointNumspS6)
annotation_data_tumor_pS6<-annotation_data_tumor_pS6[order(annotation_data_tumor_pS6$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_pS6<-rep(annotation_data_tumor_pS6,1)
# meantumor<-tumordata %>% select(Point_Num, pS6, ER, AR, Ki67)
# pS6_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(pS6), mean = mean(pS6), se = sqrt(var(pS6)/length(pS6)))
annotated_Freqs_tumor_pS6<-cbind(pS6_Freqs_tumor, anno_data_tumor_pS6)
# write.csv(annotated_Freqs_tumor_pS6, file="200707_pS6_freqs_per_sample_annotated.csv",row.names = FALSE)



########################################
################ GLUT1 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
GLUT1_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$GLUT1status))
names(GLUT1_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_GLUT1<-aggregate(GLUT1_Freqs_tumor$count, by=list(Category=GLUT1_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(GLUT1_Freqs_tumor$SampleID)) {
  frequencies <- GLUT1_Freqs_tumor[GLUT1_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_GLUT1[cell_totals_tumor_GLUT1$Category==i,2]
  GLUT1_Freqs_tumor[GLUT1_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_GLUT1<-info.csv
# get list of pointnNums is frequency data
pointNumsGLUT1<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_GLUT1 <- droplevels(annotation_data_tumor_GLUT1[annotation_data_tumor_GLUT1$PointNumber %in% pointNumsGLUT1, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_GLUT1$PointNumber <- factor(annotation_data_tumor_GLUT1$PointNumber, levels=pointNumsGLUT1)
annotation_data_tumor_GLUT1<-annotation_data_tumor_GLUT1[order(annotation_data_tumor_GLUT1$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_GLUT1<-rep(annotation_data_tumor_GLUT1,1)
# meantumor<-tumordata %>% select(Point_Num, GLUT1, ER, AR, Ki67)
# GLUT1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(GLUT1), mean = mean(GLUT1), se = sqrt(var(GLUT1)/length(GLUT1)))
annotated_Freqs_tumor_GLUT1<-cbind(GLUT1_Freqs_tumor, anno_data_tumor_GLUT1)
# write.csv(annotated_Freqs_tumor_GLUT1, file="200707_GLUT1_freqs_per_sample_annotated.csv",row.names = FALSE)



########################################
################ HIF1a ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
HIF1a_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$HIF1astatus))
names(HIF1a_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_HIF1a<-aggregate(HIF1a_Freqs_tumor$count, by=list(Category=HIF1a_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(HIF1a_Freqs_tumor$SampleID)) {
  frequencies <- HIF1a_Freqs_tumor[HIF1a_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_HIF1a[cell_totals_tumor_HIF1a$Category==i,2]
  HIF1a_Freqs_tumor[HIF1a_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_HIF1a<-info.csv
# get list of pointnNums is frequency data
pointNumsHIF1a<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_HIF1a <- droplevels(annotation_data_tumor_HIF1a[annotation_data_tumor_HIF1a$PointNumber %in% pointNumsHIF1a, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_HIF1a$PointNumber <- factor(annotation_data_tumor_HIF1a$PointNumber, levels=pointNumsHIF1a)
annotation_data_tumor_HIF1a<-annotation_data_tumor_HIF1a[order(annotation_data_tumor_HIF1a$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_HIF1a<-rep(annotation_data_tumor_HIF1a,1)
# meantumor<-tumordata %>% select(Point_Num, HIF1a, ER, AR, Ki67)
# HIF1a_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(HIF1a), mean = mean(HIF1a), se = sqrt(var(HIF1a)/length(HIF1a)))
annotated_Freqs_tumor_HIF1a<-cbind(HIF1a_Freqs_tumor, anno_data_tumor_HIF1a)
# write.csv(annotated_Freqs_tumor_HIF1a, file="200707_HIF1a_freqs_per_sample_annotated.csv",row.names = FALSE)



########################################
################ COX2 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
COX2_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$COX2status))
names(COX2_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_COX2<-aggregate(COX2_Freqs_tumor$count, by=list(Category=COX2_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(COX2_Freqs_tumor$SampleID)) {
  frequencies <- COX2_Freqs_tumor[COX2_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_COX2[cell_totals_tumor_COX2$Category==i,2]
  COX2_Freqs_tumor[COX2_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_COX2<-info.csv
# get list of pointnNums is frequency data
pointNumsCOX2<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_COX2 <- droplevels(annotation_data_tumor_COX2[annotation_data_tumor_COX2$PointNumber %in% pointNumsCOX2, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_COX2$PointNumber <- factor(annotation_data_tumor_COX2$PointNumber, levels=pointNumsCOX2)
annotation_data_tumor_COX2<-annotation_data_tumor_COX2[order(annotation_data_tumor_COX2$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_COX2<-rep(annotation_data_tumor_COX2,1)
# meantumor<-tumordata %>% select(Point_Num, COX2, ER, AR, Ki67)
# COX2_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(COX2), mean = mean(COX2), se = sqrt(var(COX2)/length(COX2)))
annotated_Freqs_tumor_COX2<-cbind(COX2_Freqs_tumor, anno_data_tumor_COX2)
# write.csv(annotated_Freqs_tumor_COX2, file="200707_COX2_freqs_per_sample_annotated.csv",row.names = FALSE)



########################################
################ CD36 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
CD36_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$CD36status))
names(CD36_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_CD36<-aggregate(CD36_Freqs_tumor$count, by=list(Category=CD36_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(CD36_Freqs_tumor$SampleID)) {
  frequencies <- CD36_Freqs_tumor[CD36_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_CD36[cell_totals_tumor_CD36$Category==i,2]
  CD36_Freqs_tumor[CD36_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_CD36<-info.csv
# get list of pointnNums is frequency data
pointNumsCD36<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_CD36 <- droplevels(annotation_data_tumor_CD36[annotation_data_tumor_CD36$PointNumber %in% pointNumsCD36, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_CD36$PointNumber <- factor(annotation_data_tumor_CD36$PointNumber, levels=pointNumsCD36)
annotation_data_tumor_CD36<-annotation_data_tumor_CD36[order(annotation_data_tumor_CD36$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_CD36<-rep(annotation_data_tumor_CD36,1)
# meantumor<-tumordata %>% select(Point_Num, CD36, ER, AR, Ki67)
# CD36_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(CD36), mean = mean(CD36), se = sqrt(var(CD36)/length(CD36)))
annotated_Freqs_tumor_CD36<-cbind(CD36_Freqs_tumor, anno_data_tumor_CD36)
# write.csv(annotated_Freqs_tumor_CD36, file="200707_CD36_freqs_per_sample_annotated.csv",row.names = FALSE)



########################################
################ CD44 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
CD44_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$CD44status))
names(CD44_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_CD44<-aggregate(CD44_Freqs_tumor$count, by=list(Category=CD44_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(CD44_Freqs_tumor$SampleID)) {
  frequencies <- CD44_Freqs_tumor[CD44_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_CD44[cell_totals_tumor_CD44$Category==i,2]
  CD44_Freqs_tumor[CD44_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_CD44<-info.csv
# get list of pointnNums is frequency data
pointNumsCD44<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_CD44 <- droplevels(annotation_data_tumor_CD44[annotation_data_tumor_CD44$PointNumber %in% pointNumsCD44, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_CD44$PointNumber <- factor(annotation_data_tumor_CD44$PointNumber, levels=pointNumsCD44)
annotation_data_tumor_CD44<-annotation_data_tumor_CD44[order(annotation_data_tumor_CD44$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_CD44<-rep(annotation_data_tumor_CD44,1)
# meantumor<-tumordata %>% select(Point_Num, CD44, ER, AR, Ki67)
# CD44_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(CD44), mean = mean(CD44), se = sqrt(var(CD44)/length(CD44)))
annotated_Freqs_tumor_CD44<-cbind(CD44_Freqs_tumor, anno_data_tumor_CD44)
# write.csv(annotated_Freqs_tumor_CD44, file="200707_CD44_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ PD1 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
PD1_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$PD1status))
names(PD1_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_PD1<-aggregate(PD1_Freqs_tumor$count, by=list(Category=PD1_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(PD1_Freqs_tumor$SampleID)) {
  frequencies <- PD1_Freqs_tumor[PD1_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_PD1[cell_totals_tumor_PD1$Category==i,2]
  PD1_Freqs_tumor[PD1_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_PD1<-info.csv
# get list of pointnNums is frequency data
pointNumsPD1<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_PD1 <- droplevels(annotation_data_tumor_PD1[annotation_data_tumor_PD1$PointNumber %in% pointNumsPD1, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_PD1$PointNumber <- factor(annotation_data_tumor_PD1$PointNumber, levels=pointNumsPD1)
annotation_data_tumor_PD1<-annotation_data_tumor_PD1[order(annotation_data_tumor_PD1$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_PD1<-rep(annotation_data_tumor_PD1,1)
# meantumor<-tumordata %>% select(Point_Num, PD1, ER, AR, Ki67)
# PD1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(PD1), mean = mean(PD1), se = sqrt(var(PD1)/length(PD1)))
annotated_Freqs_tumor_PD1<-cbind(PD1_Freqs_tumor, anno_data_tumor_PD1)
# write.csv(annotated_Freqs_tumor_PD1, file="200707_PD1_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ PDL1 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
PDL1_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$PDL1status))
names(PDL1_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_PDL1<-aggregate(PDL1_Freqs_tumor$count, by=list(Category=PDL1_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(PDL1_Freqs_tumor$SampleID)) {
  frequencies <- PDL1_Freqs_tumor[PDL1_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_PDL1[cell_totals_tumor_PDL1$Category==i,2]
  PDL1_Freqs_tumor[PDL1_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_PDL1<-info.csv
# get list of pointnNums is frequency data
pointNumsPDL1<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_PDL1 <- droplevels(annotation_data_tumor_PDL1[annotation_data_tumor_PDL1$PointNumber %in% pointNumsPDL1, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_PDL1$PointNumber <- factor(annotation_data_tumor_PDL1$PointNumber, levels=pointNumsPDL1)
annotation_data_tumor_PDL1<-annotation_data_tumor_PDL1[order(annotation_data_tumor_PDL1$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_PDL1<-rep(annotation_data_tumor_PDL1,1)
# meantumor<-tumordata %>% select(Point_Num, PDL1, ER, AR, Ki67)
# PDL1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(PDL1), mean = mean(PDL1), se = sqrt(var(PDL1)/length(PDL1)))
annotated_Freqs_tumor_PDL1<-cbind(PDL1_Freqs_tumor, anno_data_tumor_PDL1)
# write.csv(annotated_Freqs_tumor_PDL1, file="200707_PDL1_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ IDO1 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
IDO1_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$IDO1status))
names(IDO1_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_IDO1<-aggregate(IDO1_Freqs_tumor$count, by=list(Category=IDO1_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(IDO1_Freqs_tumor$SampleID)) {
  frequencies <- IDO1_Freqs_tumor[IDO1_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_IDO1[cell_totals_tumor_IDO1$Category==i,2]
  IDO1_Freqs_tumor[IDO1_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_IDO1<-info.csv
# get list of pointnNums is frequency data
pointNumsIDO1<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_IDO1 <- droplevels(annotation_data_tumor_IDO1[annotation_data_tumor_IDO1$PointNumber %in% pointNumsIDO1, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_IDO1$PointNumber <- factor(annotation_data_tumor_IDO1$PointNumber, levels=pointNumsIDO1)
annotation_data_tumor_IDO1<-annotation_data_tumor_IDO1[order(annotation_data_tumor_IDO1$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_IDO1<-rep(annotation_data_tumor_IDO1,1)
# meantumor<-tumordata %>% select(Point_Num, IDO1, ER, AR, Ki67)
# IDO1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(IDO1), mean = mean(IDO1), se = sqrt(var(IDO1)/length(IDO1)))
annotated_Freqs_tumor_IDO1<-cbind(IDO1_Freqs_tumor, anno_data_tumor_IDO1)
# write.csv(annotated_Freqs_tumor_IDO1, file="200707_IDO1_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ GZMB ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
GZMB_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$GZMBstatus))
names(GZMB_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_GZMB<-aggregate(GZMB_Freqs_tumor$count, by=list(Category=GZMB_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(GZMB_Freqs_tumor$SampleID)) {
  frequencies <- GZMB_Freqs_tumor[GZMB_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_GZMB[cell_totals_tumor_GZMB$Category==i,2]
  GZMB_Freqs_tumor[GZMB_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_GZMB<-info.csv
# get list of pointnNums is frequency data
pointNumsGZMB<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_GZMB <- droplevels(annotation_data_tumor_GZMB[annotation_data_tumor_GZMB$PointNumber %in% pointNumsGZMB, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_GZMB$PointNumber <- factor(annotation_data_tumor_GZMB$PointNumber, levels=pointNumsGZMB)
annotation_data_tumor_GZMB<-annotation_data_tumor_GZMB[order(annotation_data_tumor_GZMB$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_GZMB<-rep(annotation_data_tumor_GZMB,1)
# meantumor<-tumordata %>% select(Point_Num, GZMB, ER, AR, Ki67)
# GZMB_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(GZMB), mean = mean(GZMB), se = sqrt(var(GZMB)/length(GZMB)))
annotated_Freqs_tumor_GZMB<-cbind(GZMB_Freqs_tumor, anno_data_tumor_GZMB)
# write.csv(annotated_Freqs_tumor_GZMB, file="200707_GZMB_freqs_per_sample_annotated.csv",row.names = FALSE)


########################################
################ MMP9 ###################
########################################
##..Create a dataframe with the counts of each cell cluster across point number..##
MMP9_Freqs_tumor <- as.data.frame(table(cell_data$Point_Num, cell_data$MMP9status))
names(MMP9_Freqs_tumor) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##
cell_totals_tumor_MMP9<-aggregate(MMP9_Freqs_tumor$count, by=list(Category=MMP9_Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##
for(i in unique(MMP9_Freqs_tumor$SampleID)) {
  frequencies <- MMP9_Freqs_tumor[MMP9_Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor_MMP9[cell_totals_tumor_MMP9$Category==i,2]
  MMP9_Freqs_tumor[MMP9_Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}
##..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
annotation_data_tumor_MMP9<-info.csv
# get list of pointnNums is frequency data
pointNumsMMP9<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_MMP9 <- droplevels(annotation_data_tumor_MMP9[annotation_data_tumor_MMP9$PointNumber %in% pointNumsMMP9, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_MMP9$PointNumber <- factor(annotation_data_tumor_MMP9$PointNumber, levels=pointNumsMMP9)
annotation_data_tumor_MMP9<-annotation_data_tumor_MMP9[order(annotation_data_tumor_MMP9$CohortNumber),]
# cast the annotations to the frequency data
anno_data_tumor_MMP9<-rep(annotation_data_tumor_MMP9,1)
# meantumor<-tumordata %>% select(Point_Num, MMP9, ER, AR, Ki67)
# MMP9_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(MMP9), mean = mean(MMP9), se = sqrt(var(MMP9)/length(MMP9)))
annotated_Freqs_tumor_MMP9<-cbind(MMP9_Freqs_tumor, anno_data_tumor_MMP9)
# write.csv(annotated_Freqs_tumor_MMP9, file="200707_MMP9_freqs_per_sample_annotated.csv",row.names = FALSE)







#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

########## MERGE IMPORTANT FREQUENCIES WITH MEAN TABLE ##########

cell_totals<-aggregate(ER_Freqs_tumor$count, by=list(Category=ER_Freqs_tumor$SampleID), FUN=sum)
names(cell_totals) <- c("SampleID","cell_total")
keep_total<-c("SampleID","cell_total")
info <- merge(info.csv, cell_totals[ , keep_total, ], by.x="PointNumber", by.y="SampleID", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv


write.csv(info, file="200709_Info_celltotal.csv",row.names = FALSE)

Freqs_ER_merge<-drop.levels(annotated_Freqs_tumor_ER[annotated_Freqs_tumor_ER$cell_type %in% c(1),])
names(Freqs_ER_merge) <- c("SampleID","ERtype","ER_pos_count","ER_pos_freq")
keep_ER<-c("SampleID","ER_pos_count","ER_pos_freq")
info_add <- merge(info, Freqs_ER_merge[ , keep_ER, ], by.x="PointNumber", by.y="SampleID", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv

Freqs_HER2_merge<-drop.levels(annotated_Freqs_tumor_HER2[annotated_Freqs_tumor_HER2$cell_type %in% c(1),])
names(Freqs_HER2_merge) <- c("SampleID","HER2type","HER2_pos_count","HER2_pos_freq")
keep_HER2<-c("SampleID","HER2_pos_count","HER2_pos_freq")
info_add <- merge(info_add, Freqs_HER2_merge[ , keep_HER2, ], by.x="PointNumber", by.y="SampleID", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv

Freqs_HER2intense_merge<-drop.levels(annotated_Freqs_tumor_HER2[annotated_Freqs_tumor_HER2$cell_type %in% c(2),])
names(Freqs_HER2intense_merge) <- c("SampleID","HER2intensetype","HER2intense_pos_count","HER2intense_pos_freq")
keep_HER2intense<-c("SampleID","HER2intense_pos_count","HER2intense_pos_freq")
info_add <- merge(info_add, Freqs_HER2intense_merge[ , keep_HER2intense, ], by.x="PointNumber", by.y="SampleID", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv

Freqs_AR_merge<-drop.levels(annotated_Freqs_tumor_AR[annotated_Freqs_tumor_AR$cell_type %in% c(1),])
names(Freqs_AR_merge) <- c("SampleID","ARtype","AR_pos_count","AR_pos_freq")
keep_AR<-c("SampleID","AR_pos_count","AR_pos_freq")
info_add <- merge(info_add, Freqs_AR_merge[ , keep_AR, ], by.x="PointNumber", by.y="SampleID", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv

Freqs_Ki67_merge<-drop.levels(annotated_Freqs_tumor_Ki67[annotated_Freqs_tumor_Ki67$cell_type %in% c(1),])
names(Freqs_Ki67_merge) <- c("SampleID","Ki67type","Ki67_pos_count","Ki67_pos_freq")
keep_Ki67<-c("SampleID","Ki67_pos_count","Ki67_pos_freq")
info_add <- merge(info_add, Freqs_Ki67_merge[ , keep_Ki67, ], by.x="PointNumber", by.y="SampleID", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv


#### ADD "status" of ER HER2 etc

info_add_status <- info_add
info_add_status$ERstatus <- 0
info_add_status$HER2status <- 0
info_add_status$ARstatus <- 0
info_add_status$Ki67status <- 0

ER_pos_row_idx<-which(info_add_status$ER_pos_freq > 0.01) #use positivity threshold of 0.5
ER_rows<-rownames(data[ER_pos_row_idx,])
info_add_status[ER_rows,]$ERstatus<-"1"

HER2_pos_row_idx<-which(info_add_status$HER2intense_pos_freq > 0.10) #use positivity threshold of 0.5
HER2_rows<-rownames(data[HER2_pos_row_idx,])
info_add_status[HER2_rows,]$HER2status<-"1"

AR_pos_row_idx<-which(info_add_status$AR_pos_freq > 0.01) #use positivity threshold of 0.5
AR_rows<-rownames(data[AR_pos_row_idx,])
info_add_status[AR_rows,]$ARstatus<-"1"

Ki67_pos_row_idx<-which(info_add_status$Ki67_pos_freq > 0.14) #use positivity threshold of 0.5
Ki67_rows<-rownames(data[Ki67_pos_row_idx,])
info_add_status[Ki67_rows,]$Ki67status<-"1"

########### add mean median summaries ########


#Calculate means
meantumor<-cell_data %>% select(Point_Num, HER2, ER, AR, ECAD, Ki67, pS6, GLUT1, HIF1a, COX2, CD36, CD44, PD1, PDL1, IDO1, GZMB, MMP9)
#COLI, MPO, VIM, SMA, CK5, ECAD, PanKRT, FOXP3, CD56
ER_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(ER), mean = mean(ER), se = sqrt(var(ER)/length(ER)))
names(ER_summary_data) <- c("Point_Num","ER_median","ER_mean","ER_se")

HER2_summary_data <- meantumor %>% group_by(Point_Num) %>% summarize(median = median(HER2), mean = mean(HER2), se = sqrt(var(HER2)/length(HER2)))
names(HER2_summary_data) <- c("Point_Num","HER2_median","HER2_mean","HER2_se")

AR_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(AR), mean = mean(AR), se = sqrt(var(AR)/length(AR)))
names(AR_summary_data) <- c("Point_Num","AR_median","AR_mean","AR_se")

ECAD_summary_data <- meantumor %>% group_by(Point_Num) %>% summarize(median = median(ECAD), mean = mean(ECAD), se = sqrt(var(ECAD)/length(ECAD)))
names(ECAD_summary_data) <- c("Point_Num","ECAD_median","ECAD_mean","ECAD_se")

Ki67_summary_data <- meantumor %>% group_by(Point_Num) %>% summarize(median = median(Ki67), mean = mean(Ki67), se = sqrt(var(Ki67)/length(Ki67)))
names(Ki67_summary_data) <- c("Point_Num","Ki67_median","Ki67_mean","Ki67_se")

pS6_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(pS6), mean = mean(pS6), se = sqrt(var(pS6)/length(pS6)))
names(pS6_summary_data) <- c("Point_Num","pS6_median","pS6_mean","pS6_se")

GLUT1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(GLUT1), mean = mean(GLUT1), se = sqrt(var(GLUT1)/length(GLUT1)))
names(GLUT1_summary_data) <- c("Point_Num","GLUT1_median","GLUT1_mean","GLUT1_se")

HIF1a_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(HIF1a), mean = mean(HIF1a), se = sqrt(var(HIF1a)/length(HIF1a)))
names(HIF1a_summary_data) <- c("Point_Num","HIF1a_median","HIF1a_mean","HIF1a_se")

COX2_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(COX2), mean = mean(COX2), se = sqrt(var(COX2)/length(COX2)))
names(COX2_summary_data) <- c("Point_Num","COX2_median","COX2_mean","COX2_se")

CD36_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(CD36), mean = mean(CD36), se = sqrt(var(CD36)/length(CD36)))
names(CD36_summary_data) <- c("Point_Num","CD36_median","CD36_mean","CD36_se")

CD44_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(CD44), mean = mean(CD44), se = sqrt(var(CD44)/length(CD44)))
names(CD44_summary_data) <- c("Point_Num","CD44_median","CD44_mean","CD44_se")

PD1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(PD1), mean = mean(PD1), se = sqrt(var(PD1)/length(PD1)))
names(PD1_summary_data) <- c("Point_Num","PD1_median","PD1_mean","PD1_se")

PDL1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(PDL1), mean = mean(PDL1), se = sqrt(var(PDL1)/length(PDL1)))
names(PDL1_summary_data) <- c("Point_Num","PDL1_median","PDL1_mean","PDL1_se")

IDO1_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(IDO1), mean = mean(IDO1), se = sqrt(var(IDO1)/length(IDO1)))
names(IDO1_summary_data) <- c("Point_Num","IDO1_median","IDO1_mean","IDO1_se")

MMP9_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(MMP9), mean = mean(MMP9), se = sqrt(var(MMP9)/length(MMP9)))
names(MMP9_summary_data) <- c("Point_Num","MMP9_median","MMP9_mean","MMP9_se")

GZMB_summary_data<-meantumor %>% group_by(Point_Num) %>% summarize(median = median(GZMB), mean = mean(GZMB), se = sqrt(var(GZMB)/length(GZMB)))
names(GZMB_summary_data) <- c("Point_Num","GZMB_median","GZMB_mean","GZMB_se")







#### ADD mean median data onto info.csv for later plotting

InfoAdd_summary_aggregate <- merge(info_add_status, ER_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv


Subtype_summary_aggregate <- merge(InfoAdd_summary_aggregate, HER2_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, AR_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, ECAD_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, Ki67_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, pS6_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, GLUT1_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, HIF1a_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, COX2_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, CD36_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, CD44_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, PD1_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, PDL1_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, IDO1_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, MMP9_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv
Subtype_summary_aggregate <- merge(Subtype_summary_aggregate, GZMB_summary_data[, , ], by.x="PointNumber", by.y="Point_Num", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv





write.csv(Subtype_summary_aggregate, file="200707_Info_ALL.csv",row.names = FALSE)






























#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
##############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
setwd("~/RDATA/DCIS/200428_FRANKENSTEIN5/")

data_plot<-read.csv("200707_Info_All_Lineage.csv")


library(ggplot2)
library(forcats)
library(dplyr)
library(rstatix)
library(ggpubr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(5)




#plot_data<-droplevels(data_immune[data_immune$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-drop.levels(data_plot[data_plot$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
#plot_data<-drop.levels(plot_data[plot_data$ERstatus %in% c(1),])
#plot_data<-drop.levels(plot_data[plot_data$HER2status %in% c(0),])


#plot_data<-droplevels(data_immune[data_immune$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$Lineage),plot_data$cell_total,.fun=median,.desc=FALSE))
plot_data$Lineage <- factor(plot_data$Lineage, levels=cluster_order)
plot_data<-plot_data[order(plot_data$Lineage),]



theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 40, hjust=1))

#  Ki67    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_Ki67_mean.pdf", height =5, width = 10) 
Tissue_Ki67_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=Ki67_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_Ki67_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_Ki67_Mean.pdf", height =4, width = 10) 
Ki67_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=Ki67_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Ki67_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_Ki67_Mean_INVASIVEONLY.pdf", height =4, width = 10) 
Ki67_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=Ki67_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Ki67_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################



#  pS6    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_pS6_mean.pdf", height =5, width = 10) 
Tissue_pS6_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=pS6_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_pS6_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_pS6_Mean", height =4, width = 10) 
pS6_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=pS6_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
pS6_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_pS6_Mean_INVASIVEONLY", height =4, width = 10) 
pS6_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=pS6_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
pS6_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  GLUT1    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_GLUT1_mean", height =5, width = 10) 
Tissue_GLUT1_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=GLUT1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_GLUT1_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_GLUT1_Mean", height =4, width = 10) 
GLUT1_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=GLUT1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
GLUT1_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_GLUT1_Mean_INVASIVEONLY", height =4, width = 10) 
GLUT1_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=GLUT1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
GLUT1_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  COX2    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_COX2_mean", height =5, width = 10) 
Tissue_COX2_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=COX2_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_COX2_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_COX2_Mean", height =4, width = 10) 
COX2_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=COX2_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
COX2_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_COX2_Mean_INVASIVEONLY", height =4, width = 10) 
COX2_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=COX2_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
COX2_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  CD36    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_CD36_mean", height =5, width = 10) 
Tissue_CD36_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=CD36_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_CD36_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_CD36_Mean", height =4, width = 10) 
CD36_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=CD36_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
CD36_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_CD36_Mean_INVASIVEONLY", height =4, width = 10) 
CD36_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=CD36_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
CD36_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  HIF1a    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_HIF1a_mean", height =5, width = 10) 
Tissue_HIF1a_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=HIF1a_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_HIF1a_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_HIF1a_Mean", height =4, width = 10) 
HIF1a_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=HIF1a_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
HIF1a_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_HIF1a_Mean_INVASIVEONLY", height =4, width = 10) 
HIF1a_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=HIF1a_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
HIF1a_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  IDO1    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_IDO1_mean", height =5, width = 10) 
Tissue_IDO1_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=IDO1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_IDO1_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_IDO1_Mean", height =4, width = 10) 
IDO1_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=IDO1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
IDO1_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_IDO1_Mean_INVASIVEONLY", height =4, width = 10) 
IDO1_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=IDO1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
IDO1_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  PDL1    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_PDL1_mean", height =5, width = 10) 
Tissue_PDL1_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=PDL1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_PDL1_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_PDL1_Mean", height =4, width = 10) 
PDL1_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=PDL1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
PDL1_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_PDL1_Mean_INVASIVEONLY", height =4, width = 10) 
PDL1_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=PDL1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
PDL1_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  PD1    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_PD1_mean", height =5, width = 10) 
Tissue_PD1_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=PD1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_PD1_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_PD1_Mean", height =4, width = 10) 
PD1_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=PD1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
PD1_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_PD1_Mean_INVASIVEONLY", height =4, width = 10) 
PD1_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=PD1_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
PD1_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  MMP9    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_MMP9_mean", height =5, width = 10) 
Tissue_MMP9_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=MMP9_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_MMP9_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_MMP9_Mean", height =4, width = 10) 
MMP9_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=MMP9_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
MMP9_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_MMP9_Mean_INVASIVEONLY", height =4, width = 10) 
MMP9_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=MMP9_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
MMP9_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  GZMB    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_GZMB_mean", height =5, width = 10) 
Tissue_GZMB_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=GZMB_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_GZMB_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_GZMB_Mean", height =4, width = 10) 
GZMB_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=GZMB_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
GZMB_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_GZMB_Mean_INVASIVEONLY", height =4, width = 10) 
GZMB_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=GZMB_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
GZMB_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  CD44    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_CD44_mean", height =5, width = 10) 
Tissue_CD44_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=CD44_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_CD44_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_CD44_Mean", height =4, width = 10) 
CD44_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=CD44_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
CD44_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_CD44_Mean_INVASIVEONLY", height =4, width = 10) 
CD44_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=CD44_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
CD44_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  ER    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_ER_mean", height =5, width = 10) 
Tissue_ER_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=ER_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_ER_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_ER_Mean", height =4, width = 10) 
ER_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=ER_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
ER_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_ER_Mean_INVASIVEONLY", height =4, width = 10) 
ER_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=ER_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
ER_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  HER2    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_HER2_mean", height =5, width = 10) 
Tissue_HER2_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=HER2_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_HER2_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_HER2_Mean", height =4, width = 10) 
HER2_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=HER2_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
HER2_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_HER2_Mean_INVASIVEONLY", height =4, width = 10) 
HER2_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=HER2_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
HER2_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################


#  AR    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_AR_mean", height =5, width = 10) 
Tissue_AR_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=AR_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_AR_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_AR_Mean", height =4, width = 10) 
AR_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=AR_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
AR_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_AR_Mean_INVASIVEONLY", height =4, width = 10) 
AR_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=AR_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
AR_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################

#  ECAD    ######################################################################################################################################################################################
#############################################################################################################################################################################################
my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'concurrent'), c('normal', 'NewInv'), c('normal', 'NewCIS'), c('DCIS', 'concurrent'), c('DCIS', 'NewInv'), c('DCIS', 'NewCIS'), c('concurrent', 'NewInv'), c('NewInv', 'NewCIS'))
#
pdf(file = "200708_Tissue_ECAD_mean", height =5, width = 10) 
Tissue_ECAD_mean<-ggplot(plot_data,aes(x=factor(Tissue, level = c('normal', 'DCIS','concurrent','NewInv','NewCIS')), y=ECAD_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Tissue') + 
  labs(y = 'Mean') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
Tissue_ECAD_mean
dev.off()
#### CASE CTRL - all Immunes ####
my_comparisons <- list( c('normal', 'case'), c('normal', 'ctrl'), c('ctrl', 'case'))
#
pdf(file = "200708_Status_ECAD_Mean.pdf", height =4, width = 10) 
ECAD_status_means<-ggplot(plot_data,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=ECAD_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
ECAD_status_means
dev.off()
#### INVASIVE ONLY CASE CTRL #####
marker_inv <- plot_data[plot_data$Reccurence %in% c( "continv","ipsinv","na"),]
marker_inv2 <- marker_inv[!marker_inv$Tissue_Type %in% c( "Ips ILC","Cont ILC"),]
#
pdf(file = "200708_Status_ECAD_Mean_INVASIVEONLY", height =4, width = 10) 
ECAD_inv_status_means<-ggplot(marker_inv2,aes(x=factor(Status, level = c('normal', 'ctrl', 'case')), y=ECAD_mean)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method= "wilcox.test") +
  geom_boxplot(aes(fill = Lineage), outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_fill_manual(values=color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Immune Phenotypes') + 
  labs(y = 'Frequency') +
  facet_wrap(.~Lineage, scales='free_y', ncol=5) +
  theme(legend.position = 'none')
ECAD_inv_status_means
dev.off()
##########################################################################################################################################################################################################################################################################################################################################################################################
#############################################################################################################################################################################################



















