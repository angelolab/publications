# TBMIBIcomparePDL1PD1FreqTNBC.R
# Author: Erin McCaffrey 
# Date created: 190822
# Overview: This script reads in the normalized expression matrix for the mycobacteria and sarcoid samples.
# It also reads in the cell-size normalized data for a tonsil, placenta, and spleen sample from the cohort. 
# The TNBC data (cell-size normalized, asinh transformed, z-scored) is also imported. First the immune controls
# are linearly scaled (x100), asinh-transformed, and scaled to the 99th percentile. The TNBC cohort is scaled 
# 0-1. Then the frequency of PD1, Lag3, and PDL1+ cells is determined for each sample. There are two variations
# of this: 1. Using the cutoff of 0.5 for the TNBC cohort as published and 2. Using the custom thresholds for each
# marker decided for my cohort. 

library(dplyr)
library(reshape2)
library(ggplot2)

#####..Load in mycobacteria cohort and TNBC data..#####

dataAll_CSscaled<-read.csv("data/allsamples_dataCS3px_annotated.csv") #cell size normalized and untransformed for immune controls
dataGran_norm<-read.csv("data/allTB-sarcoid-scdata.csv") #cell size normalized, linearly + asinh scaled, 99th percentile scaled
dataTNBC_zscore<-read.csv("data/TNBC_cellData.csv") #TNBC data, cell size normalized, asinh, z-scored
dataTNBCimm_zscore<-droplevels(dataTNBC_zscore[dataTNBC_zscore$Group==2,])

# specify points for each category
imm_tons <- c(2,10,17,24,31,38,45,51,56,78,83,86) 
imm_spl <-c(1,9,16,23,30,37,44,62)
imm_plac <- c(3,11,18,25,32,46,63)
granA_lung <- c(21,84,42,88,28,89,14,15,98,99,90,91,94,95,96,97) 
granA_pleura <- c(33,34,26,27,40,61)
granA_endo <- c(47,48)
granA_LN <- c(54,55,92,93)
granA_vert <- c(6,7)
gran_sarc <- c(67,68,69,70,71,72,73,74,75,76)
TNBC <- unique(dataTNBC_zscore$SampleID)

# for easy indexing define larger sample groups where needed
all_myco <- c(granA_lung,granA_pleura,granA_LN,granA_endo,granA_vert) #all A granulomas
all_imm <- c(imm_tons,imm_spl,imm_plac) #all immune controls

#####..Transform the immune control data..#####

#get just immune data
data_control <- droplevels(dataAll_CSscaled[dataAll_CSscaled$SampleID %in% all_imm,])

# linear transform the data by 100 (only the expression data and not sample ID, cell label, cell size, or tissue)

data_control[,4:50]<- data_control[,4:50]*100

# arcsinh transform data (only the expression data and not sample ID, cell label, cell size, or tissue)

asinh_scale <- 5
data_trans<-data_control
data_trans[,4:50]<- asinh(data_control[,4:50]/ asinh_scale)

#####..Percent normalization and scale data 0-99th percentile..#####

v <- 1:1000
v
quantile_value <- 0.999
quantile(v, quantile_value)

percentile.vector <- apply(data_trans[,4:50], 2, function(x) quantile(x, quantile_value, names = F))
percentile.vector
data_imm_norm<-data_trans
data_imm_norm[,4:50] <- data.frame(t(t(data_trans[,4:50]) / as.numeric(percentile.vector)))
data_imm_norm<-data_imm_norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)


#####..For each sample get the % positive for PD1, PDL1, Lag3, and IDO1..#####

##..granulomas..##

PD1_thresh = 0.21
PDL1_thresh = 0.25
Lag3_thresh = 0.09

#PD1
data_PD1<-data.frame(table(dataGran_norm$SampleID))
PD1pos<-dataGran_norm[dataGran_norm$PD.1 >= PD1_thresh,]
data_PD1pos<-data.frame(table(PD1pos$SampleID))

data_PD1$PD1pos<-0
data_PD1[data_PD1$Var1 %in% data_PD1pos$Var1,]$PD1pos<-data_PD1pos$Freq
names(data_PD1)<-c("SampleID","Total","Total_PD1pos")
data_PD1$percentPD1<-as.numeric(format((data_PD1$Total_PD1pos / data_PD1$Total)*100),digits=3)

#PDL1
data_PDL1<-data.frame(table(dataGran_norm$SampleID))
PDL1pos<-dataGran_norm[dataGran_norm$PD.L1 >= PDL1_thresh,]
data_PDL1pos<-data.frame(table(PDL1pos$SampleID))

data_PDL1$PDL1pos<-0
data_PDL1[data_PDL1$Var1 %in% data_PDL1pos$Var1,]$PDL1pos<-data_PDL1pos$Freq
names(data_PDL1)<-c("SampleID","Total","Total_PDL1pos")
data_PDL1$percentPDL1<-as.numeric(format((data_PDL1$Total_PDL1pos / data_PDL1$Total)*100),digits=3)

#Lag3
data_Lag3<-data.frame(table(dataGran_norm$SampleID))
Lag3pos<-dataGran_norm[dataGran_norm$Lag3 >= Lag3_thresh,]
data_Lag3pos<-data.frame(table(Lag3pos$SampleID))

data_Lag3$Lag3pos<-0
data_Lag3[data_Lag3$Var1 %in% data_Lag3pos$Var1,]$Lag3pos<-data_Lag3pos$Freq
names(data_Lag3)<-c("SampleID","Total","Total_Lag3pos")
data_Lag3$percentLag3<-as.numeric(format((data_Lag3$Total_Lag3pos / data_Lag3$Total)*100),digits=3)

#concatenate
freq_gran<-cbind(data_PD1,data_PDL1[,-c(1,2)],data_Lag3[,-c(1,2)])


##..immune controls..##

#PD1
data_PD1<-data.frame(table(data_imm_norm$SampleID))
PD1pos<-data_imm_norm[data_imm_norm$PD.1 >= PD1_thresh,]
data_PD1pos<-data.frame(table(PD1pos$SampleID))

data_PD1$PD1pos<-0
data_PD1[data_PD1$Var1 %in% data_PD1pos$Var1,]$PD1pos<-data_PD1pos$Freq
names(data_PD1)<-c("SampleID","Total","Total_PD1pos")
data_PD1$percentPD1<-as.numeric(format((data_PD1$Total_PD1pos / data_PD1$Total)*100),digits=3)

#PDL1
data_PDL1<-data.frame(table(data_imm_norm$SampleID))
PDL1pos<-data_imm_norm[data_imm_norm$PD.L1 >= PDL1_thresh,]
data_PDL1pos<-data.frame(table(PDL1pos$SampleID))

data_PDL1$PDL1pos<-0
data_PDL1[data_PDL1$Var1 %in% data_PDL1pos$Var1,]$PDL1pos<-data_PDL1pos$Freq
names(data_PDL1)<-c("SampleID","Total","Total_PDL1pos")
data_PDL1$percentPDL1<-as.numeric(format((data_PDL1$Total_PDL1pos / data_PDL1$Total)*100),digits=3)

#Lag3
data_Lag3<-data.frame(table(data_imm_norm$SampleID))
Lag3pos<-data_imm_norm[data_imm_norm$Lag3 >= Lag3_thresh,]
data_Lag3pos<-data.frame(table(Lag3pos$SampleID))

data_Lag3$Lag3pos<-0
data_Lag3[data_Lag3$Var1 %in% data_Lag3pos$Var1,]$Lag3pos<-data_Lag3pos$Freq
names(data_Lag3)<-c("SampleID","Total","Total_Lag3pos")
data_Lag3$percentLag3<-as.numeric(format((data_Lag3$Total_Lag3pos / data_Lag3$Total)*100),digits=3)

#concatenate
freq_imm<-cbind(data_PD1,data_PDL1[,-c(1,2)],data_Lag3[,-c(1,2)])

##..TNBC..##

thresh = 0.50

#PD1
data_PD1<-data.frame(table(dataTNBCimm_zscore$SampleID))
PD1pos<-dataTNBCimm_zscore[dataTNBCimm_zscore$PD1 >= thresh,]
data_PD1pos<-data.frame(table(PD1pos$SampleID))

data_PD1$PD1pos<-0
data_PD1[data_PD1$Var1 %in% data_PD1pos$Var1,]$PD1pos<-data_PD1pos$Freq
names(data_PD1)<-c("SampleID","Total","Total_PD1pos")
data_PD1$percentPD1<-as.numeric(format((data_PD1$Total_PD1pos / data_PD1$Total)*100),digits=3)

#PDL1
data_PDL1<-data.frame(table(dataTNBCimm_zscore$SampleID))
PDL1pos<-dataTNBCimm_zscore[dataTNBCimm_zscore$PD.L1 >= thresh,]
data_PDL1pos<-data.frame(table(PDL1pos$SampleID))

data_PDL1$PDL1pos<-0
data_PDL1[data_PDL1$Var1 %in% data_PDL1pos$Var1,]$PDL1pos<-data_PDL1pos$Freq
names(data_PDL1)<-c("SampleID","Total","Total_PDL1pos")
data_PDL1$percentPDL1<-as.numeric(format((data_PDL1$Total_PDL1pos / data_PDL1$Total)*100),digits=3)

#Lag3
data_Lag3<-data.frame(table(dataTNBCimm_zscore$SampleID))
Lag3pos<-dataTNBCimm_zscore[dataTNBCimm_zscore$Lag3 >= thresh,]
data_Lag3pos<-data.frame(table(Lag3pos$SampleID))

data_Lag3$Lag3pos<-0
data_Lag3[data_Lag3$Var1 %in% data_Lag3pos$Var1,]$Lag3pos<-data_Lag3pos$Freq
names(data_Lag3)<-c("SampleID","Total","Total_Lag3pos")
data_Lag3$percentLag3<-as.numeric(format((data_Lag3$Total_Lag3pos / data_Lag3$Total)*100),digits=3)


freq_TNBC<-cbind(data_PD1,data_PDL1[,-c(1,2)],data_Lag3[,-c(1,2)])

#####..Add tissue types to each dataframe before merging into one..#####

freq_TNBC$tissue<-'TNBC'
freq_imm<-freq_imm %>% mutate(tissue=case_when(freq_imm$SampleID %in% imm_tons ~ "tonsil",
                                               freq_imm$SampleID %in% imm_spl ~ "spleen",
                                               freq_imm$SampleID %in% imm_plac ~ "placenta"))
freq_gran<-freq_gran %>% mutate(tissue=case_when(freq_gran$SampleID %in% all_myco ~ "mycobacteria",
                                                 freq_gran$SampleID %in% gran_sarc ~ "sarcoid"))

all_samples_freq<-rbind(freq_TNBC,freq_imm,freq_gran)

####..Plot paired boxplots of frequency of lag3+, pd1+, and pdl1+ cells per condition..####

dfmelt.freq <- melt(all_samples_freq, id.vars = "tissue", measure.vars = c("percentPD1", "percentPDL1","percentLag3"))
dfmelt.count <- melt(all_samples_freq, id.vars = "tissue", measure.vars = c("Total_PD1pos", "Total_PDL1pos","Total_Lag3pos"))

boxplot.freq <- ggplot(dfmelt.freq, aes(x=tissue, y=value, fill=variable))+
  geom_boxplot()+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, vjust=0.4)) +
  labs(y="% Positive for Marker") + 
  labs(x="Tissue Type") 
boxplot.freq

scatter.freq <- ggplot(all_samples_freq, aes(x=percentPD1, y=percentPDL1, color=tissue))+
  geom_point(size=3)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, vjust=0.4)) +
  labs(y="% PDL1+") + 
  labs(x="% PD1+") 
scatter.freq

multiboxplot.freq <- ggplot(dfmelt.freq, aes(x=tissue, y=value, fill=variable))+
  geom_boxplot()+
  facet_wrap(.~variable,scales=c("free_y")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, vjust=0.4)) +
  labs(y="% Positive for Marker") + 
  labs(x="Tissue Type")
multiboxplot.freq

####..Plot the percent lymphs positivie for PD1 and Lag3 in sarcoid, TNBC, and TB...####

# Sarcoid
sarc<-c(67,68,69,70,71,72,73,74,75,76)
data_sarc<-droplevels(dataGran_norm[dataGran_norm$SampleID %in% sarc & dataGran_norm$cell_lin=="lymphocyte",])
sarcoid_summary<-as.data.frame(table(data_sarc$SampleID))
colnames(sarcoid_summary)<-c('SampleID','Total')
sarcoid_summary$PD1pos<-(table(data_sarc$SampleID,data_sarc$PD.1 > PD1_thresh )[,2])
sarcoid_summary$Lag3pos<-(table(data_sarc$SampleID,data_sarc$Lag3 > Lag3_thresh)[,2])
sarcoid_summary$Freq_PD1<-as.numeric(format((sarcoid_summary$PD1pos / sarcoid_summary$Total)*100),digits=3)
sarcoid_summary$Freq_Lag3<-as.numeric(format((sarcoid_summary$Lag3pos / sarcoid_summary$Total)*100),digits=3)
sarcoid_summary$Tissue<-'sarcoid'


# TB
data_TB<-droplevels(dataGran_norm[dataGran_norm$SampleID %in% all_myco & dataGran_norm$cell_lin=="lymphocyte",])
TB_summary<-as.data.frame(table(data_TB$SampleID))
colnames(TB_summary)<-c('SampleID','Total')
TB_summary$PD1pos<-(table(data_TB$SampleID,data_TB$PD.1 > PD1_thresh )[,2])
TB_summary$Lag3pos<-(table(data_TB$SampleID,data_TB$Lag3 > Lag3_thresh)[,2])
TB_summary$Freq_PD1<-as.numeric(format((TB_summary$PD1pos / TB_summary$Total)*100),digits=3)
TB_summary$Freq_Lag3<-as.numeric(format((TB_summary$Lag3pos / TB_summary$Total)*100),digits=3)
TB_summary$Tissue<-'TB'

# TNBC
lymphs<-c(1,2,3,4,6)
data_TNBC<-droplevels(dataTNBCimm_zscore[dataTNBCimm_zscore$immuneGroup %in% lymphs,])
TNBC_summary<-as.data.frame(table(data_TNBC$SampleID))
colnames(TNBC_summary)<-c('SampleID','Total')
TNBC_summary$PD1pos<-(table(data_TNBC$SampleID,data_TNBC$PD1 > 0.5 )[,2])
TNBC_summary$Lag3pos<-(table(data_TNBC$SampleID,data_TNBC$Lag3 > 0.5)[,2])
TNBC_summary$Freq_PD1<-as.numeric(format((TNBC_summary$PD1pos / TNBC_summary$Total)*100),digits=3)
TNBC_summary$Freq_Lag3<-as.numeric(format((TNBC_summary$Lag3pos / TNBC_summary$Total)*100),digits=3)
TNBC_summary$Tissue<-'TNBC'

# PD1  and Lag3 summary all tissue
all_tissues_summary<-rbind(sarcoid_summary, TB_summary, TNBC_summary)

# Get PD1 summmary stats


# Plot

# plot_data <- all_tissues_summary[!all_tissues_summary$Tissue=='sarcoid',]
plot_data <- all_tissues_summary

PD1<-ggplot(plot_data, aes(x=Tissue, y=Freq_PD1))+
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  geom_point() +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", aes(width=0.25)) +
  # stat_compare_means(method = "wilcox.test", label = "p.signif") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons=list(c(1,3), c(1,2), c(2,3))) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, vjust=0.4)) +
  labs(y="Percent Pos") + 
  labs(x="Tissue Type")
PD1

compare_means(Freq_PD1 ~ Tissue, plot_data, method= "wilcox.test")

Lag3<-ggplot(plot_data, aes(x=Tissue, y=Freq_Lag3))+
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  geom_point() +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",aes(width=0.25)) +
  # stat_compare_means(method = "wilcox.test", label = "p.signif") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons=list(c(1,3), c(1,2), c(2,3))) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, vjust=0.4)) +
  labs(y="Percent Pos") + 
  labs(x="Tissue Type")
Lag3

compare_means(Freq_Lag3 ~ Tissue, plot_data, method= "wilcox.test")

####..For each sample get ratio of Lag3/PDL1 and PD1/PDL1..####

all_samples_freq$Lag3.PDL1.Ratio <- log2(all_samples_freq$Total_Lag3pos/all_samples_freq$Total_PDL1pos)
all_samples_freq$PD1.PDL1.Ratio <- log2(all_samples_freq$Total_PD1pos/all_samples_freq$Total_PDL1pos)
all_samples_freq$Lag3.PD1.Ratio <- log(all_samples_freq$Total_Lag3pos/all_samples_freq$Total_PD1pos)

# Annotate new v old data
all_samples_freq$cohort<-'old'
all_samples_freq[all_samples_freq$SampleID %in% c(90,91,92,93,94,95,96,97,98,99),]$cohort<-'new'

dfmelt.ratio <- melt(all_samples_freq, id.vars = c("tissue","cohort"), measure.vars = c("Lag3.PDL1.Ratio", "PD1.PDL1.Ratio","Lag3.PD1.Ratio"))


boxplot.ratio <- ggplot(dfmelt.ratio, aes(x=tissue, y=value, fill=variable))+
  geom_boxplot(outlier.shape=NA)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, vjust=0.4)) +
  labs(y="Ratio") + 
  labs(x="Tissue Type") +
  facet_wrap(.~variable,scales=c("free_y"))
boxplot.ratio


plot_data<-all_samples_freq[all_samples_freq$tissue %in% c('mycobacteria','TNBC'),]
plot_data<-plot_data %>% filter((!is.infinite(PD1.PDL1.Ratio)))
boxplot.individualratio <- ggplot(plot_data, aes(x=tissue, y=PD1.PDL1.Ratio, fill = tissue))+
  geom_boxplot(fill='grey', outlier.shape = NA)+
  geom_point(size = 3, position = position_jitterdodge()) +
  geom_hline(yintercept=c(0), linetype="dashed") +
  theme_bw() + 
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme(legend.position = 'none') +
  labs(y="Ratio") + 
  labs(x="Tissue Type") 
boxplot.individualratio






