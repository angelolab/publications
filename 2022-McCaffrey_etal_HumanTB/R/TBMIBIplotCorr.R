# TBMIBIplotCorr.R
# Author: Erin McCaffrey 
# Date created: 190319
# Overview: This script reads in the csv for cell-size normalized data and then plots two markers, 
# determines their Pearson correlation and fits with a linear regression. 

library(ggplot2)

##..Import data..##

data<-read.csv("data/allTB-sarcoid-scdata.csv")

##..Keep TB samples and myeloid cells only..##

data_myco<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN','gran_vert'),])
myeloid<-c("CD14_Mono","CD11b/c_CD206_Mac/Mono","CD11c_DC/Mono","CD68_Mac","CD16_CD14_Mono",
           "CD206_Mac","CD163_Mac","CD209_DC","giant_cell","mast","neutrophil")
data_myeloid<-droplevels(data_myco[data_myco$cell_type %in% myeloid, ])

##..Run pearson correlation and get coefficient..##

corr<-cor.test(data_myco$IDO,data_myco$PD.L1,method="pearson")
corrCoeff<-corr$estimate

corr_myeloid<-cor.test(data_myeloid$IDO,data_myeloid$PD.L1,method="pearson")
corrCoeff_myeloid<-corr_myeloid$estimate

##..Break down correlation by cohort and tissue..##

resection<-c(21,84,42,88,28,89)
corr_resection<-cor.test(data_myco[data_myco$SampleID %in% resection,]$IDO,
                         data_myco[data_myco$SampleID %in% resection,]$PD.L1,method="pearson")
corrCoeff_resection<-corr_resection$estimate

diagnostic_pulm<-c(14,15,98,99)
corr_diagnostic_pulm<-cor.test(data_myco[data_myco$SampleID %in% diagnostic_pulm,]$IDO,
                         data_myco[data_myco$SampleID %in% diagnostic_pulm,]$PD.L1,method="pearson")
corrCoeff_diagnostic_pulm<-corr_diagnostic_pulm$estimate

autopsy<-c(90,91,94,95,96,97)
corr_autopsy<-cor.test(data_myco[data_myco$SampleID %in% autopsy,]$IDO,
                         data_myco[data_myco$SampleID %in% autopsy,]$PD.L1,method="pearson")
corrCoeff_autopsy<-corr_autopsy$estimate

pulm<-c(resection, diagnostic_pulm, autopsy)
corr_pulm<-cor.test(data_myco[data_myco$SampleID %in% pulm,]$IDO,
                         data_myco[data_myco$SampleID %in% pulm,]$PD.L1,method="pearson")
corrCoeff_pulm<-corr_pulm$estimate

diagnostic_expulm<-c(6,7,33,34,26,27,40,61,47,48,54,55,92,93)
corr_diagnostic_expulm<-cor.test(data_myco[data_myco$SampleID %in% diagnostic_expulm,]$IDO,
                         data_myco[data_myco$SampleID %in% diagnostic_expulm,]$PD.L1,method="pearson")
corrCoeff_diagnostic_expulm<-corr_diagnostic_expulm$estimate

diagnostic<-c(diagnostic_pulm, diagnostic_expulm)
corr_diagnostic<-cor.test(data_myco[data_myco$SampleID %in% diagnostic,]$IDO,
                                 data_myco[data_myco$SampleID %in% diagnostic,]$PD.L1,method="pearson")
corrCoeff_diagnostic<-corr_diagnostic$estimate

##..Plot IDO and PD.L1 with correlation and linear regression..##

corrCoeff_txt<-as.numeric(round(corrCoeff_myeloid,digits=2))
corrString=toString(corrCoeff_txt)


ggplot(data_myeloid,aes(x=IDO,y=PD.L1)) +
  geom_point() + 
  labs(x="IDO Expression") + 
  labs(y="PD-L1 Expression") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


