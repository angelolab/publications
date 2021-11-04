## DCIS_CombineAndPlotMaskData ###

# Author: Tyler Risom 
# Date created: 200407
# Overview: Combine an annotated cell table with the a cell table showing the occupancy of cells in different masks and generate summary statistics and plots.
# Masks produced from MatLab "Generate MasksTR "MIBIcreateRegionMasksTR" and MIBIannotateCellsinMaskRegions"scripts

#If you need to make any custom combinations of mask location values, edit the .csv first e.g if you want to count all in mask 1 + 2, create a sum column of mask 1 and 2. All nonzeros occupy the mask

###########################################
##..Install packages and open libraries..##
###########################################

source("~/Risom2022_R_scripts/Risom2022_SC_FUNCTIONS")

library(dplyr)
library(viridis)
library(forcats)
library(reshape2)
library(tidyr)
library(mefa)

#setwd("/Volumes/MIBI Files/remote/dataPerCell")
setwd("~/Desktop/DCIS/Segmentation/200303_CLEANED_RESEGMENT/")
info.csv <- read.csv(file = paste0(info.dir, "200317_COHORT_METADATA.csv"))

###Combine Mask Data####################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #

data_cohort<-read.csv("200406_CellTable_TumorImmuneStromaMyoepOtherNeut_ClusterCoded.csv") #import annotated cell table csv

data_mask<-read.csv("200407_cell_mask_annotations.csv") #import cell mask locations file

data_cohort_masklocations<-cbind(data_cohort, data_mask) #merge files


write.csv(data_cohort_masklocations, file="200407_CellTable_TumorImmuneStromaMyoepOtherNeut_ClusterCoded_MaskLocations.csv",row.names = FALSE) #write merged

# CHECK THE CSV AND THAT ALL CELL LABELS ALLIGN !!!!

###### Drop cells not occupying specific masks #####

selectmask<-c(1) 
deselectmask<-c(0) 

data_cohort_masklocations_Emask<-droplevels(data_cohort_masklocations[data_cohort_masklocations$duct_mask %in% selectmask,]) #drop the 
write.csv(data_cohort_masklocations_Emask, file="200407_CellTable_Emask_Only.csv",row.names = FALSE) #write merged

data_cohort_masklocations_Smask<-droplevels(data_cohort_masklocations[!data_cohort_masklocations$S_mask %in% deselectmask,]) #drop the 
write.csv(data_cohort_masklocations_Smask, file="200407_CellTable_Smask_Only.csv",row.names = FALSE) #write merged

data_cohort_masklocations_MYOEPmask<-droplevels(data_cohort_masklocations[data_cohort_masklocations$myoep_mask %in% selectmask,]) #drop the 
write.csv(data_cohort_masklocations_MYOEPmask, file="200407_CellTable_MYOEPmask_Only.csv",row.names = FALSE) #write merged

data_cohort_masklocations_PERIPHmask<-droplevels(data_cohort_masklocations[data_cohort_masklocations$periph_mask %in% selectmask,]) #drop the 
write.csv(data_cohort_masklocations_PERIPHmask, file="200407_CellTable_PERIPHmask_Only.csv",row.names = FALSE) #write merged

data_cohort_masklocations_DISTALHmask<-droplevels(data_cohort_masklocations[data_cohort_masklocations$distal_mask %in% selectmask,]) #drop the 
write.csv(data_cohort_masklocations_DISTALHmask, file="200407_CellTable_DISTALmask_Only.csv",row.names = FALSE) #write merged




### Summary Stats  #####################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
##########################################################################################################################################################################################


datamask<-read.csv("200407_CellTable_Emask_Only.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_mask<-datamask %>% select(Point_Num, sublineage)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_Emask <- as.data.frame(table(cell_data_mask$Point_Num, cell_data_mask$sublineage))
# Freqs_Emask <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$fs_clusters_stroma))
names(Freqs_Emask) <- c("SampleID","sublineage","count")

##..Get overall totals on a per sample basis..##

cell_totals_Emask<-aggregate(Freqs_Emask$count, by=list(Category=Freqs_Emask$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_Emask$SampleID)) {
  frequencies <- Freqs_Emask[Freqs_Emask$SampleID==i,"count"] / cell_totals_Emask[cell_totals_Emask$Category==i,2]
  Freqs_Emask[Freqs_Emask$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data_Emask<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_mask$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort_Emask <- droplevels(annotation_data_Emask[annotation_data_Emask$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort_Emask$PointNumber <- factor(annotation_data_cohort_Emask$PointNumber, levels=pointNums)
annotation_data_cohort_Emask<-annotation_data_cohort_Emask[order(annotation_data_cohort_Emask$PointNumber),]
# cast the annotations to the frequency data
anno_data_Emask<-rep(annotation_data_cohort_Emask,1)
annotated_Freqs_Emask<-cbind(Freqs_Emask, anno_data_Emask)

##..Save as a csv later use and plotting..##

# cell totals per sample
write.csv(cell_totals_Emask, file="200407_Emask_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_Emask, file="200407_Emask_Freqs_Emask_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_Emask, file="200407_Emask_Freqs_Emask_per_sample_annotated.csv",row.names = FALSE)

### DROP EPITHELIAL CELLS ####


datamask_nonep <- droplevels(datamask[!datamask$sublineage %in% c('TUMOR','MYOEP','OTHER'),])

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_mask_nonep<-datamask_nonep %>% select(Point_Num, sublineage)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_Emask_Nonep <- as.data.frame(table(cell_data_mask_nonep$Point_Num, cell_data_mask_nonep$sublineage))

names(Freqs_Emask_Nonep) <- c("SampleID","sublineage","count")
cell_totals_Emask_nonep<-aggregate(Freqs_Emask_Nonep$count, by=list(Category=Freqs_Emask_Nonep$SampleID), FUN=sum)
# annotated_Freqs_Emask_Nonepithelial <- droplevels(annotated_Freqs_Emask[!annotated_Freqs_Emask$sublineage %in% c('TUMOR','MYOEP'),])
##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_Emask_Nonep$SampleID)) {
  frequencies <- Freqs_Emask_Nonep[Freqs_Emask_Nonep$SampleID==i,"count"] / cell_totals_Emask_nonep [cell_totals_Emask_nonep $Category==i,2]
  Freqs_Emask_Nonep[Freqs_Emask_Nonep$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data_Emask_nonep <-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_mask_nonep$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort_Emask_nonep <- droplevels(annotation_data_Emask_nonep [annotation_data_Emask_nonep $PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort_Emask_nonep $PointNumber <- factor(annotation_data_cohort_Emask_nonep $PointNumber, levels=pointNums)
annotation_data_cohort_Emask_nonep <-annotation_data_cohort_Emask_nonep [order(annotation_data_cohort_Emask_nonep $PointNumber),]
# cast the annotations to the frequency data
anno_data_Emask_nonep <-rep(annotation_data_cohort_Emask_nonep ,1)
annotated_Freqs_Emask_nonep <-cbind(Freqs_Emask_Nonep , anno_data_Emask_nonep )

write.csv(cell_totals_Emask_nonep, file="200407_Emask_nonep_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_Emask_Nonep, file="200407_Emask_Freqs_Emask_nonep_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations

write.csv(annotated_Freqs_Emask_nonep , file="200407_Emask_Nonep_Freqs_per_sample_annotated.csv",row.names = FALSE)

################################ PLOT IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(16)

##..Import data..##
# setwd("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/")
# data_meta_clusters30<-read.csv("200303_DCIScohortFrankenstein4_ScaleTransNormtrimbad_FLOWSOMMetacluster20_ANNOTATED.csv")
data<-read.csv("200407_Emask_Nonep_Freqs_per_sample_annotated.csv")
plot_data<-data
plot_data<-droplevels(data[data$Tissue %in% c('DCIS','concurrent'),])

#plot_data_nonep<-droplevels(data[!data$sublineage %in% c('TUMOR','MYOEP'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$sublineage),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$sublineage <- factor(plot_data$sublineage, levels=cluster_order)
plot_data<-plot_data[order(plot_data$sublineage),]
# reorder by descending median frequency (nonep clusters)
# cluster_order_nonep<-levels(fct_reorder(as.factor(plot_data_nonep$sublineage),plot_data_nonep$frequency,.fun=median,.desc=FALSE))
# plot_data_nonep$sublineage <- factor(plot_data_nonep$sublineage, levels=cluster_order_nonep)
# plot_data_nonep<-plot_data_nonep[order(plot_data_nonep$sublineage),]
# reorder by descending median frequency (DCIS clusters)
# cluster_order_DCIS<-levels(fct_reorder(as.factor(plot_data_DCIS$sublineage),plot_data_DCIS$frequency,.fun=median,.desc=FALSE))
# plot_data_DCIS$sublineage <- factor(plot_data_DCIS$sublineage, levels=cluster_order_nonep)
# plot_data_DCIS<-plot_data_DCIS[order(plot_data_DCIS$sublineage),]

##..Plot all clusters across all samples..##
pdf(file = "200407_Emask__nonep_DCIS_Cellbox", height = 12, width = 16) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

# pdf(file = "200407_Emask_DCIS_NonEp_Cellbox", height = 12, width = 16) 
# cell_box<-ggplot(plot_data_nonep, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
#   geom_boxplot() +
#   scale_fill_manual(values = rev(color)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
#   labs(x="Cluster") +
#   labs(y="Frequency of Total") +
#   ggtitle("Frequency of Cell Types") +
#   guides(fill=guide_legend(title="Cell Type"))
# cell_box
# dev.off()

##..Plot frequency broken down by Point Number..##
pdf(file = "200407_Emask_NONEP_nonep_Cellbar", height = 12, width = 40)
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(sublineage))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust=1)) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
cell_bar
dev.off()

##..Plot each cluster frequency as box across conditions..##

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 20, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = sublineage)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~sublineage, scales = "free")
# session
pdf(file = "200407_Emask_CaseCtrl", height = 15, width = 15)
status<-ggplot(plot_data[!plot_data$Status %in% c('continv','ipsinv','contcis','ipscis','tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
status
dev.off()


pdf(file = "200407_Emask_Freq_Progression", height = 6, width = 7)
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS', 'NewCIS', 'concurrent','NewInv')), y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
Tissue
dev.off()


plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])


pdf(file = "200407_Emask_CaseCtrl_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

pdf(file = "200407_Emask_CaseCtrlConcurrent_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case','concurrent')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()


## S MASK ############################################################################################################################################################################################################
##########################################################################################################################################################################################


datamask<-read.csv("200407_CellTable_Smask_Only.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_mask<-datamask %>% select(Point_Num, sublineage)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_Smask <- as.data.frame(table(cell_data_mask$Point_Num, cell_data_mask$sublineage))
# Freqs_Smask  <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$fs_clusters_stroma))
names(Freqs_Smask) <- c("SampleID","sublineage","count")

##..Get overall totals on a per sample basis..##

cell_totals_Smask<-aggregate(Freqs_Smask $count, by=list(Category=Freqs_Smask$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_Smask$SampleID)) {
  frequencies <- Freqs_Smask[Freqs_Smask$SampleID==i,"count"] / cell_totals_Smask[cell_totals_Smask$Category==i,2]
  Freqs_Smask[Freqs_Smask$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data_Smask<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_mask$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort_Smask <- droplevels(annotation_data_Smask[annotation_data_Smask$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort_Smask$PointNumber <- factor(annotation_data_cohort_Smask$PointNumber, levels=pointNums)
annotation_data_cohort_Smask<-annotation_data_cohort_Smask[order(annotation_data_cohort_Smask$PointNumber),]
# cast the annotations to the frequency data
anno_data_Smask<-rep(annotation_data_cohort_Smask,1)
annotated_Freqs_Smask<-cbind(Freqs_Smask, anno_data_Smask)

##..Save as a csv later use and plotting..##

# cell totals per sample
write.csv(cell_totals_Smask, file="200407_Smask_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_Smask, file="200407_Smask_Freqs_Smask_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_Smask, file="200407_Smask_Freqs_Smask_per_sample_annotated.csv",row.names = FALSE)


####### PLOT SMASK DATA #######


data<-read.csv("200407_Smask_Freqs_Smask_per_sample_annotated.csv")
plot_data<-data
plot_data<-droplevels(data[data$Tissue %in% c('DCIS','concurrent'),])

#plot_data_nonep<-droplevels(data[!data$sublineage %in% c('TUMOR','MYOEP'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$sublineage),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$sublineage <- factor(plot_data$sublineage, levels=cluster_order)
plot_data<-plot_data[order(plot_data$sublineage),]
# reorder by descending median frequency (nonep clusters)
# cluster_order_nonep<-levels(fct_reorder(as.factor(plot_data_nonep$sublineage),plot_data_nonep$frequency,.fun=median,.desc=FALSE))
# plot_data_nonep$sublineage <- factor(plot_data_nonep$sublineage, levels=cluster_order_nonep)
# plot_data_nonep<-plot_data_nonep[order(plot_data_nonep$sublineage),]
# reorder by descending median frequency (DCIS clusters)
# cluster_order_DCIS<-levels(fct_reorder(as.factor(plot_data_DCIS$sublineage),plot_data_DCIS$frequency,.fun=median,.desc=FALSE))
# plot_data_DCIS$sublineage <- factor(plot_data_DCIS$sublineage, levels=cluster_order_nonep)
# plot_data_DCIS<-plot_data_DCIS[order(plot_data_DCIS$sublineage),]

##..Plot all clusters across all samples..##
pdf(file = "200407_Smask_DCIS_Cellbox", height = 12, width = 16) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

# pdf(file = "200407_Smask_DCIS_NonEp_Cellbox", height = 12, width = 16) 
# cell_box<-ggplot(plot_data_nonep, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
#   geom_boxplot() +
#   scale_fill_manual(values = rev(color)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
#   labs(x="Cluster") +
#   labs(y="Frequency of Total") +
#   ggtitle("Frequency of Cell Types") +
#   guides(fill=guide_legend(title="Cell Type"))
# cell_box
# dev.off()

##..Plot frequency broken down by Point Number..##
pdf(file = "200407_Smask_NONEP_nonep_Cellbar", height = 12, width = 40)
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(sublineage))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust=1)) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
cell_bar
dev.off()

##..Plot each cluster frequency as box across conditions..##

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 20, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = sublineage)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~sublineage, scales = "free")
# session
pdf(file = "200407_Smask_CaseCtrl", height = 15, width = 15)
status<-ggplot(plot_data[!plot_data$Status %in% c('continv','ipsinv','contcis','ipscis','tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
status
dev.off()

pdf(file = "200407_Smask_FreqNEWEVENT", height = 15, width = 15)
reccur<-ggplot(plot_data[!plot_data$Status %in% c('ctrl', 'case', 'tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
reccur
dev.off()



pdf(file = "200407_Smask_Freq_Progression", height = 6, width = 7)
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS', 'NewCIS', 'concurrent','NewInv')), y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
Tissue
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])



pdf(file = "200407_Smask_CaseCtrl_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

pdf(file = "200407_Smask_CaseCtrlConcurrent_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case','concurrent')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

sublin1<-c('TUMOR')
sublin2<-c('MYOEP')

table(datamask[datamask$sublineage==sublin1,]$Point_Num)
table(datamask[datamask$sublineage==sublin2,]$Point_Num)







## MYOEP MASK ############################################################################################################################################################################################################
##########################################################################################################################################################################################

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(16)


datamask<-read.csv("200407_CellTable_MYOEPmask_Only.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_mask<-datamask %>% select(Point_Num, sublineage)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_MYOEP <- as.data.frame(table(cell_data_mask$Point_Num, cell_data_mask$sublineage))
# Freqs_MYOEP  <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$fs_clusters_stroma))
names(Freqs_MYOEP) <- c("SampleID","sublineage","count")

##..Get overall totals on a per sample basis..##

cell_totals_MYOEP<-aggregate(Freqs_MYOEP $count, by=list(Category=Freqs_MYOEP$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_MYOEP$SampleID)) {
  frequencies <- Freqs_MYOEP[Freqs_MYOEP$SampleID==i,"count"] / cell_totals_MYOEP[cell_totals_MYOEP$Category==i,2]
  Freqs_MYOEP[Freqs_MYOEP$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data_MYOEP<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_mask$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort_MYOEP <- droplevels(annotation_data_MYOEP[annotation_data_MYOEP$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort_MYOEP$PointNumber <- factor(annotation_data_cohort_MYOEP$PointNumber, levels=pointNums)
annotation_data_cohort_MYOEP<-annotation_data_cohort_MYOEP[order(annotation_data_cohort_MYOEP$PointNumber),]
# cast the annotations to the frequency data
anno_data_MYOEP<-rep(annotation_data_cohort_MYOEP,1)
annotated_Freqs_MYOEP<-cbind(Freqs_MYOEP, anno_data_MYOEP)

##..Save as a csv later use and plotting..##

# cell totals per sample
write.csv(cell_totals_MYOEP, file="200407_MYOEP_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_MYOEP, file="200407_MYOEP_Freqs_MYOEP_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_MYOEP, file="200407_MYOEP_Freqs_MYOEP_per_sample_annotated.csv",row.names = FALSE)


####### PLOT MYOEP DATA #######


data<-read.csv("200407_MYOEP_Freqs_MYOEP_per_sample_annotated.csv")
plot_data<-data
plot_data<-droplevels(data[data$Tissue %in% c('DCIS','concurrent'),])

#plot_data_nonep<-droplevels(data[!data$sublineage %in% c('TUMOR','MYOEP'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$sublineage),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$sublineage <- factor(plot_data$sublineage, levels=cluster_order)
plot_data<-plot_data[order(plot_data$sublineage),]
# reorder by descending median frequency (nonep clusters)
# cluster_order_nonep<-levels(fct_reorder(as.factor(plot_data_nonep$sublineage),plot_data_nonep$frequency,.fun=median,.desc=FALSE))
# plot_data_nonep$sublineage <- factor(plot_data_nonep$sublineage, levels=cluster_order_nonep)
# plot_data_nonep<-plot_data_nonep[order(plot_data_nonep$sublineage),]
# reorder by descending median frequency (DCIS clusters)
# cluster_order_DCIS<-levels(fct_reorder(as.factor(plot_data_DCIS$sublineage),plot_data_DCIS$frequency,.fun=median,.desc=FALSE))
# plot_data_DCIS$sublineage <- factor(plot_data_DCIS$sublineage, levels=cluster_order_nonep)
# plot_data_DCIS<-plot_data_DCIS[order(plot_data_DCIS$sublineage),]

##..Plot all clusters across all samples..##
pdf(file = "200407_MYOEP_DCIS_Cellbox", height = 12, width = 16) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

# pdf(file = "200407_MYOEP_DCIS_NonEp_Cellbox", height = 12, width = 16) 
# cell_box<-ggplot(plot_data_nonep, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
#   geom_boxplot() +
#   scale_fill_manual(values = rev(color)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
#   labs(x="Cluster") +
#   labs(y="Frequency of Total") +
#   ggtitle("Frequency of Cell Types") +
#   guides(fill=guide_legend(title="Cell Type"))
# cell_box
# dev.off()

##..Plot frequency broken down by Point Number..##
pdf(file = "200407_MYOEP__Cellbar", height = 12, width = 40)
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(sublineage))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust=1)) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
cell_bar
dev.off()

##..Plot each cluster frequency as box across conditions..##

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 20, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = sublineage)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~sublineage, scales = "free")
# session
pdf(file = "200407_MYOEP_CaseCtrl", height = 15, width = 15)
status<-ggplot(plot_data[!plot_data$Status %in% c('continv','ipsinv','contcis','ipscis','tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
status
dev.off()

pdf(file = "200407_MYOEP_FreqNEWEVENT", height = 15, width = 15)
reccur<-ggplot(plot_data[!plot_data$Status %in% c('ctrl', 'case', 'tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
reccur
dev.off()



pdf(file = "200407_MYOEP_Freq_Progression", height = 6, width = 7)
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS', 'NewCIS', 'concurrent','NewInv')), y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
Tissue
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])



pdf(file = "200407_MYOEP_CaseCtrl_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

pdf(file = "200407_MYOEP_CaseCtrlConcurrent_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case','concurrent')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

sublin1<-c('TUMOR')
sublin2<-c('MYOEP')



## PERIPH MASK ############################################################################################################################################################################################################
##########################################################################################################################################################################################


datamask<-read.csv("200407_CellTable_PERIPHmask_Only.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_mask<-datamask %>% select(Point_Num, sublineage)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_PERIPH <- as.data.frame(table(cell_data_mask$Point_Num, cell_data_mask$sublineage))
# Freqs_PERIPH  <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$fs_clusters_stroma))
names(Freqs_PERIPH) <- c("SampleID","sublineage","count")

##..Get overall totals on a per sample basis..##

cell_totals_PERIPH<-aggregate(Freqs_PERIPH $count, by=list(Category=Freqs_PERIPH$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_PERIPH$SampleID)) {
  frequencies <- Freqs_PERIPH[Freqs_PERIPH$SampleID==i,"count"] / cell_totals_PERIPH[cell_totals_PERIPH$Category==i,2]
  Freqs_PERIPH[Freqs_PERIPH$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data_PERIPH<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_mask$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort_PERIPH <- droplevels(annotation_data_PERIPH[annotation_data_PERIPH$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort_PERIPH$PointNumber <- factor(annotation_data_cohort_PERIPH$PointNumber, levels=pointNums)
annotation_data_cohort_PERIPH<-annotation_data_cohort_PERIPH[order(annotation_data_cohort_PERIPH$PointNumber),]
# cast the annotations to the frequency data
anno_data_PERIPH<-rep(annotation_data_cohort_PERIPH,1)
annotated_Freqs_PERIPH<-cbind(Freqs_PERIPH, anno_data_PERIPH)

##..Save as a csv later use and plotting..##

# cell totals per sample
write.csv(cell_totals_PERIPH, file="200407_PERIPH_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_PERIPH, file="200407_PERIPH_Freqs_PERIPH_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_PERIPH, file="200407_PERIPH_Freqs_PERIPH_per_sample_annotated.csv",row.names = FALSE)


####### PLOT PERIPH DATA #######


data<-read.csv("200407_PERIPH_Freqs_PERIPH_per_sample_annotated.csv")
plot_data<-data
plot_data<-droplevels(data[data$Tissue %in% c('DCIS','concurrent', 'NewInv', 'NewCIS', 'normal'),])

#plot_data_nonep<-droplevels(data[!data$sublineage %in% c('TUMOR','PERIPH'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$sublineage),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$sublineage <- factor(plot_data$sublineage, levels=cluster_order)
plot_data<-plot_data[order(plot_data$sublineage),]
# reorder by descending median frequency (nonep clusters)
# cluster_order_nonep<-levels(fct_reorder(as.factor(plot_data_nonep$sublineage),plot_data_nonep$frequency,.fun=median,.desc=FALSE))
# plot_data_nonep$sublineage <- factor(plot_data_nonep$sublineage, levels=cluster_order_nonep)
# plot_data_nonep<-plot_data_nonep[order(plot_data_nonep$sublineage),]
# reorder by descending median frequency (DCIS clusters)
# cluster_order_DCIS<-levels(fct_reorder(as.factor(plot_data_DCIS$sublineage),plot_data_DCIS$frequency,.fun=median,.desc=FALSE))
# plot_data_DCIS$sublineage <- factor(plot_data_DCIS$sublineage, levels=cluster_order_nonep)
# plot_data_DCIS<-plot_data_DCIS[order(plot_data_DCIS$sublineage),]

##..Plot all clusters across all samples..##
pdf(file = "200407_PERIPH_DCIS_Cellbox", height = 12, width = 16) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

# pdf(file = "200407_PERIPH_DCIS_NonEp_Cellbox", height = 12, width = 16) 
# cell_box<-ggplot(plot_data_nonep, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
#   geom_boxplot() +
#   scale_fill_manual(values = rev(color)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
#   labs(x="Cluster") +
#   labs(y="Frequency of Total") +
#   ggtitle("Frequency of Cell Types") +
#   guides(fill=guide_legend(title="Cell Type"))
# cell_box
# dev.off()

##..Plot frequency broken down by Point Number..##
pdf(file = "200407_PERIPH_NONEP_nonep_Cellbar", height = 12, width = 40)
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(sublineage))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust=1)) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
cell_bar
dev.off()

##..Plot each cluster frequency as box across conditions..##

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 20, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = sublineage)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~sublineage, scales = "free")
# session
pdf(file = "200407_PERIPH_CaseCtrl", height = 15, width = 15)
status<-ggplot(plot_data[!plot_data$Status %in% c('continv','ipsinv','contcis','ipscis','tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
status
dev.off()

pdf(file = "200407_PERIPH_FreqNEWEVENT", height = 15, width = 15)
reccur<-ggplot(plot_data[!plot_data$Status %in% c('ctrl', 'case', 'tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
reccur
dev.off()



pdf(file = "200407_PERIPH_Freq_Progression", height = 6, width = 7)
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS', 'NewCIS', 'concurrent','NewInv')), y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
Tissue
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])



pdf(file = "200407_PERIPH_CaseCtrl_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

pdf(file = "200407_PERIPH_Progression", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('normal','ctrl', 'case','concurrent','continv','ipsinv')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

sublin1<-c('TUMOR')
sublin2<-c('PERIPH')






## DISTAL MASK ############################################################################################################################################################################################################
##########################################################################################################################################################################################


datamask<-read.csv("200407_CellTable_DISTALmask_Only.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_mask<-datamask %>% select(Point_Num, sublineage)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_DISTAL <- as.data.frame(table(cell_data_mask$Point_Num, cell_data_mask$sublineage))
# Freqs_DISTAL  <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$fs_clusters_stroma))
names(Freqs_DISTAL) <- c("SampleID","sublineage","count")

##..Get overall totals on a per sample basis..##

cell_totals_DISTAL<-aggregate(Freqs_DISTAL $count, by=list(Category=Freqs_DISTAL$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_DISTAL$SampleID)) {
  frequencies <- Freqs_DISTAL[Freqs_DISTAL$SampleID==i,"count"] / cell_totals_DISTAL[cell_totals_DISTAL$Category==i,2]
  Freqs_DISTAL[Freqs_DISTAL$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data_DISTAL<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_mask$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort_DISTAL <- droplevels(annotation_data_DISTAL[annotation_data_DISTAL$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort_DISTAL$PointNumber <- factor(annotation_data_cohort_DISTAL$PointNumber, levels=pointNums)
annotation_data_cohort_DISTAL<-annotation_data_cohort_DISTAL[order(annotation_data_cohort_DISTAL$PointNumber),]
# cast the annotations to the frequency data
anno_data_DISTAL<-rep(annotation_data_cohort_DISTAL,1)
annotated_Freqs_DISTAL<-cbind(Freqs_DISTAL, anno_data_DISTAL)

##..Save as a csv later use and plotting..##

# cell totals per sample
write.csv(cell_totals_DISTAL, file="200407_DISTAL_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_DISTAL, file="200407_DISTAL_Freqs_DISTAL_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_DISTAL, file="200407_DISTAL_Freqs_DISTAL_per_sample_annotated.csv",row.names = FALSE)


####### PLOT DISTAL DATA #######


data<-read.csv("200407_DISTAL_Freqs_DISTAL_per_sample_annotated.csv")
plot_data<-data
plot_data<-droplevels(data[data$Tissue %in% c('DCIS','concurrent'),])

#plot_data_nonep<-droplevels(data[!data$sublineage %in% c('TUMOR','DISTAL'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$sublineage),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$sublineage <- factor(plot_data$sublineage, levels=cluster_order)
plot_data<-plot_data[order(plot_data$sublineage),]
# reorder by descending median frequency (nonep clusters)
# cluster_order_nonep<-levels(fct_reorder(as.factor(plot_data_nonep$sublineage),plot_data_nonep$frequency,.fun=median,.desc=FALSE))
# plot_data_nonep$sublineage <- factor(plot_data_nonep$sublineage, levels=cluster_order_nonep)
# plot_data_nonep<-plot_data_nonep[order(plot_data_nonep$sublineage),]
# reorder by descending median frequency (DCIS clusters)
# cluster_order_DCIS<-levels(fct_reorder(as.factor(plot_data_DCIS$sublineage),plot_data_DCIS$frequency,.fun=median,.desc=FALSE))
# plot_data_DCIS$sublineage <- factor(plot_data_DCIS$sublineage, levels=cluster_order_nonep)
# plot_data_DCIS<-plot_data_DCIS[order(plot_data_DCIS$sublineage),]

##..Plot all clusters across all samples..##
pdf(file = "200407_DISTAL_DCIS_Cellbox", height = 12, width = 16) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

# pdf(file = "200407_DISTAL_DCIS_NonEp_Cellbox", height = 12, width = 16) 
# cell_box<-ggplot(plot_data_nonep, aes(x=fct_reorder(as.factor(sublineage),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(sublineage))) + 
#   geom_boxplot() +
#   scale_fill_manual(values = rev(color)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
#   labs(x="Cluster") +
#   labs(y="Frequency of Total") +
#   ggtitle("Frequency of Cell Types") +
#   guides(fill=guide_legend(title="Cell Type"))
# cell_box
# dev.off()

##..Plot frequency broken down by Point Number..##
pdf(file = "200407_DISTAL_NONEP_nonep_Cellbar", height = 12, width = 40)
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(sublineage))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust=1)) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
cell_bar
dev.off()

##..Plot each cluster frequency as box across conditions..##

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 20, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = sublineage)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~sublineage, scales = "free")
# session
pdf(file = "200407_DISTAL_CaseCtrl", height = 15, width = 15)
status<-ggplot(plot_data[!plot_data$Status %in% c('continv','ipsinv','contcis','ipscis','tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
status
dev.off()

pdf(file = "200407_DISTAL_FreqNEWEVENT", height = 15, width = 15)
reccur<-ggplot(plot_data[!plot_data$Status %in% c('ctrl', 'case', 'tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
reccur
dev.off()



pdf(file = "200407_DISTAL_Freq_Progression", height = 6, width = 7)
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS', 'NewCIS', 'concurrent','NewInv')), y=frequency, fill = as.factor(sublineage))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
Tissue
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])



pdf(file = "200407_DISTAL_CaseCtrl_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

pdf(file = "200407_DISTAL_CaseCtrlConcurrent_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case','concurrent')), y=frequency, fill = as.factor(sublineage))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~sublineage, scales = "free")
statusinvasive
dev.off()

sublin1<-c('TUMOR')
sublin2<-c('DISTAL')

