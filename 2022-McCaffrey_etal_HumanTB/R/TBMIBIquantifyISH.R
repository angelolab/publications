# TBMIBIquantifyISH.R
# Author: Erin McCaffrey 
# Date created: 201006
# Overview: This script reads in the summarized ISH data for all granulomas and produces
# several visualizations quantifying the results.

library(factoextra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)


# import data for whole granuloma regions
ish_data<-read.csv("data/ISH_all-cells.csv")

# drop necrosis
ish_data_nonnec<-droplevels(ish_data[!ish_data$Parent=='Necrosis',])

# Summarise total dots and total cells per region

ish_summary<-as.data.frame(table(ish_data_nonnec$Image))
colnames(ish_summary)<-c('image','total_cells')
ish_summary$Probe<-str_sub(ish_summary$image,-8,-5)

ish_total_dots <- ish_data_nonnec %>%
  group_by(Image) %>%
  summarize(total_dots = sum(Subcellular..DAB..Num.spots.estimated))

ish_summary$total_dots<-0
ish_summary[ish_summary$image %in% ish_total_dots$Image, ]$total_dots<-ish_total_dots$total_dots

# import and add region area data
ish_region_data<-read.csv('data/ISH_dot_summary_region.csv')
ish_region_imm<-droplevels(ish_region_data[ish_region_data$Name=="Immune cells",]) #only keep immune cell regions

ish_summary$Area.µm.2<-0
ish_summary[ish_summary$image %in% ish_region_imm$Image,]$Area.µm.2<-ish_region_imm$Area.µm.2
ish_summary$Area.mm.2<-ish_summary$Area.µm.2 / 1000

##### repeat all of the above for mac and lymph zones #####

ish_data_zones<-read.csv("data/ISH_all-cells_zones.csv")

# drop necrosis
ish_data_nonnec_zone<-droplevels(ish_data_zones[!ish_data_zones$Parent=='Necrosis',])

# get Mac and Lymph separated
ish_data_mac <- droplevels(ish_data_nonnec_zone[ish_data_nonnec_zone$Parent == 'Mac',])
ish_data_lymph <- droplevels(ish_data_nonnec_zone[ish_data_nonnec_zone$Parent == 'Lymph',])

# Summarise total dots and total cells per region

# Mac
ish_summary_mac<-as.data.frame(table(ish_data_mac$Image))
colnames(ish_summary_mac)<-c('image','total_cells')
ish_summary_mac$Probe<-str_sub(ish_summary_mac$image,-8,-5)

ish_total_dots_mac <- ish_data_mac %>%
  group_by(Image) %>%
  summarize(total_dots = sum(Subcellular..DAB..Num.spots.estimated))

ish_summary_mac$total_dots<-0
ish_summary_mac[ish_summary_mac$image %in% ish_total_dots_mac$Image, ]$total_dots<-ish_total_dots_mac$total_dots

# lymph
ish_summary_lymph<-as.data.frame(table(ish_data_lymph$Image))
colnames(ish_summary_lymph)<-c('image','total_cells')
ish_summary_lymph$Probe<-str_sub(ish_summary_lymph$image,-8,-5)

ish_total_dots_lymph <- ish_data_lymph %>%
  group_by(Image) %>%
  summarize(total_dots = sum(Subcellular..DAB..Num.spots.estimated))

ish_summary_lymph$total_dots<-0
ish_summary_lymph[ish_summary_lymph$image %in% ish_total_dots_lymph$Image, ]$total_dots<-ish_total_dots_lymph$total_dots

# add region area data

#Mac
ish_region_mac<-droplevels(ish_region_data[ish_region_data$Name=="Mac",]) #only keep immune cell regions
ish_summary_mac$Area.µm.2<-0
ish_summary_mac[ish_summary_mac$image %in% ish_region_mac$Image,]$Area.µm.2<-ish_region_mac$Area.µm.2
ish_summary_mac$Area.mm.2<-ish_summary_mac$Area.µm.2 / 1000

#Lymph
ish_region_lymph<-droplevels(ish_region_data[ish_region_data$Name=="Lymph",]) #only keep immune cell regions
ish_summary_lymph$Area.µm.2<-0
ish_summary_lymph[ish_summary_lymph$image %in% ish_region_lymph$Image,]$Area.µm.2<-ish_region_lymph$Area.µm.2
ish_summary_lymph$Area.mm.2<-ish_summary_lymph$Area.µm.2 / 1000

# Merge into single summary with region name
ish_summary$Class<-'Immune'
ish_summary_mac$Class<-'Mac'
ish_summary_lymph$Class<-'Lymph'


ish_summary_all<-rbind(ish_summary,ish_summary_mac)
ish_summary_all<-rbind(ish_summary_all,ish_summary_lymph)

# get dots / area and dots / cell
ish_summary_all$dots_per_mm.2<-ish_summary_all$total_dots / ish_summary_all$Area.mm.2
ish_summary_all$dots_per_cell<-ish_summary_all$total_dots / ish_summary_all$total_cells

# add unifying patient-roi indication to each row (take string before second '_')
ish_summary_all$Patient.ROI<-gsub("(.+?_.+?)_.*" ,"\\1",ish_summary_all$image)

# get total region only
ish_summary_imm<-droplevels(ish_summary_all[ish_summary_all$Class == 'Immune',])

# restructure the dataframe to have one ROI per row
ish_data_perROI_count <- reshape2::dcast(ish_summary_imm, Patient.ROI ~ Probe, value.var = "total_dots")
ish_data_perROI_perarea <- reshape2::dcast(ish_summary_imm, Patient.ROI ~ Probe, value.var = "dots_per_mm.2")
ish_data_perROI_percell <- reshape2::dcast(ish_summary_imm, Patient.ROI ~ Probe, value.var = "dots_per_cell")

# determine log2 FC TGFb: IFNg for dots per mm and dots per cell
ish_data_perROI_perarea$TGFB_IFNG_ratio <- log2(ish_data_perROI_perarea$TGFB/ish_data_perROI_perarea$IFNG)
ish_data_perROI_percell$TGFB_IFNg_ratio <- log2(ish_data_perROI_percell$TGFB/ish_data_perROI_percell$IFNG)
ish_data_perROI_count$TGFB_IFNg_ratio <- log2(ish_data_perROI_count$TGFB/ish_data_perROI_count$IFNG)

# Plot a paired barchart of all RNA counts per ROI including cell an area normalized #

plot_data<-ish_data_perROI_count[,-6]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI'))

ggplot(plot_data_melt, aes(fill=variable, y=value, x=Patient.ROI)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_grey() +
  theme_bw() +
  labs(x = 'ROI') + 
  labs(y = 'Total Dot Count') +
  theme(axis.text.x = element_text(angle=35,hjust=1)) 


plot_data<-ish_data_perROI_percell[,-6]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI'))

ggplot(plot_data_melt, aes(fill=variable, y=value, x=Patient.ROI)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_grey() +
  theme_bw() +
  labs(x = 'ROI') + 
  labs(y = 'Dots/Cell') +
  theme(axis.text.x = element_text(angle=35,hjust=1)) 

plot_data<-ish_data_perROI_perarea[,-6]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI'))

ggplot(plot_data_melt, aes(fill=variable, y=value, x=Patient.ROI)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_grey() +
  theme_bw() +
  labs(x = 'ROI') + 
  labs(y = 'Dots/mm^2') +
  theme(axis.text.x = element_text(angle=35,hjust=1)) 

# Append patient data back to summary dataframes #

ish_data_perROI_count$Patient<-c(2384,2384,2384,51210,51210,51210,72215,72215,72215,72215,72215)
ish_data_perROI_percell$Patient<-c(2384,2384,2384,51210,51210,51210,72215,72215,72215,72215,72215)
ish_data_perROI_perarea$Patient<-c(2384,2384,2384,51210,51210,51210,72215,72215,72215,72215,72215)

# Plot dots/mm^2 and dots/cell TGFb v IFNg w/ dots shape by patient #

# per area
plot_data<-ish_data_perROI_perarea[,c(1,4,5,7)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG','TGFB'))

data_summary <- plot_data_melt[,c(3,4)] %>%
  group_by(variable) %>%
  summarize_all(list(median=median))

neg_line<-median(ish_data_perROI_perarea$DapB, na.rm = T)

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  geom_hline(yintercept=neg_line, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Dots/mm^2")

# version without 72215 #

outlier_removed<-plot_data_melt[!plot_data_melt$Patient==72215, ]

data_summary <- outlier_removed[,c(3,4)] %>%
  group_by(variable) %>%
  summarize_all(list(median=median))

ggplot(data=outlier_removed, aes(x=variable, y=value)) +
  geom_point(data=outlier_removed, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  geom_hline(yintercept=neg_line, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Dots/mm^2")

#per cell
plot_data<-ish_data_perROI_percell[,c(1,4,5,7)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG','TGFB'))

data_summary <- plot_data_melt[,c(3,4)] %>%
  group_by(variable) %>%
  summarize_all(list(median=median))

neg_line<-median(ish_data_perROI_percell$DapB, na.rm = T)

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  geom_hline(yintercept=neg_line, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Dots/cell")

_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Dots/cell")


# Plot FC per cell, per count, and per area for all ROI #

ggplot(data=ish_data_perROI_count, aes(x=reorder(Patient.ROI, -TGFB_IFNg_ratio), y=TGFB_IFNg_ratio)) +
  geom_bar(stat="Identity") +
  scale_fill_grey() +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(angle=35,hjust=1)) +
  labs(y="log2(TGFB dots / IFNG dots)")


ggplot(data=ish_data_perROI_perarea, aes(x=reorder(Patient.ROI, -TGFB_IFNG_ratio), y=TGFB_IFNG_ratio)) +
  geom_bar(stat="Identity") +
  scale_fill_grey() +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(angle=35,hjust=1)) +
  labs(y="log2(TGFB dots per mm^2 / IFNG dots per mm^2")


ggplot(data=ish_data_perROI_percell, aes(x=reorder(Patient.ROI, -TGFB_IFNg_ratio), y=TGFB_IFNg_ratio)) +
  geom_bar(stat="Identity") +
  scale_fill_grey() +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(angle=35,hjust=1)) +
  labs(y="log2(TGFB dots per cell / IFNG dots per cell")

# Compare DapB to IFNG and TGFB #

# per area
plot_data<-ish_data_perROI_perarea[,c(1,3,4,5,7)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG','TGFB','DapB'))

data_summary <- plot_data_melt %>%
  group_by(variable) %>%
  summarize(median = median(value, na.rm=T))

my_comparisons<-list(c("DapB","IFNG"),
                  c("DapB","TGFB"),
                  c("IFNG","TGFB"))

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons = my_comparisons) +
  labs(y="Dots/mm^2")

# per cell
plot_data<-ish_data_perROI_percell[,c(1,3,4,5,7)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG','TGFB','DapB'))

data_summary <- plot_data_melt %>%
  group_by(variable) %>%
  summarize(median = median(value, na.rm=T))

my_comparisons<-list(c("DapB","IFNG"),
                     c("DapB","TGFB"),
                     c("IFNG","TGFB"))

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons = my_comparisons) +
  labs(y="Dots/cell")
  
# Look at correlation between IFNG and TGFB (area and cell normalized)
  
corr_cell <-cor.test(ish_data_perROI_percell$IFNG,ish_data_perROI_percell$TGFB ,method="pearson")
corrCoeff_cell<-corr_cell$estimate

ggplot(ish_data_perROI_percell,aes(x=IFNG,y=TGFB)) +
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  labs(x="IFNG dots per cell") + 
  labs(y="IFNG dots per cell") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

corr_area <-cor.test(ish_data_perROI_perarea$IFNG,ish_data_perROI_perarea$TGFB ,method="pearson")
corrCoeff_area<-corr_area$estimate

ggplot(ish_data_perROI_perarea,aes(x=IFNG,y=TGFB)) +
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  labs(x="IFNG dots per mm^2") + 
  labs(y="IFNG dots per mm^2") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


# Analyze myeloid core and lymphocytic cuff differently

# get just the mac and lymph data
ish_summary_zones<-droplevels(ish_summary_all[ish_summary_all$Class %in% c('Mac','Lymph'),])

# add zone specific probe names
ish_summary_zones[ish_summary_zones$Class=='Mac',]$Probe<-paste(ish_summary_zones[ish_summary_zones$Class=='Mac',]$Probe,'_Mac', sep='')
ish_summary_zones[ish_summary_zones$Class=='Lymph',]$Probe<-paste(ish_summary_zones[ish_summary_zones$Class=='Lymph',]$Probe,'_Lymph', sep='')

# restructure the dataframe to have one ROI per row
ish_data_perROI_zone_count <- reshape2::dcast(ish_summary_zones, Patient.ROI ~ Probe, value.var = "total_dots")
ish_data_perROI_zone_perarea <- reshape2::dcast(ish_summary_zones, Patient.ROI ~ Probe, value.var = "dots_per_mm.2")
ish_data_perROI_zone_percell <- reshape2::dcast(ish_summary_zones, Patient.ROI ~ Probe, value.var = "dots_per_cell")

# Append patient data back to summary dataframes #

ish_data_perROI_zone_count$Patient<-c(2384,2384,2384,51210,51210,51210,72215,72215,72215,72215,72215)
ish_data_perROI_zone_perarea$Patient<-c(2384,2384,2384,51210,51210,51210,72215,72215,72215,72215,72215)
ish_data_perROI_zone_percell$Patient<-c(2384,2384,2384,51210,51210,51210,72215,72215,72215,72215,72215)

# Plot data

#per area
plot_data<-ish_data_perROI_zone_perarea[,c(1,6,7,8,9,10)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG_Lymph',
                                                                                        'IFNG_Mac',
                                                                                        'TGFB_Lymph',
                                                                                        'TGFB_Mac'))

data_summary <- plot_data_melt[,c(3,4)] %>%
  group_by(variable) %>%
  summarize_all(list(median=median))

neg_line1<-median(ish_data_perROI_zone_perarea$DapB_Lymph, na.rm = T)
neg_line2<-median(ish_data_perROI_zone_perarea$DapB_Mac, na.rm = T)

my_comparisons = list(c('IFNG_Lymph',
                        'IFNG_Mac'),
                      c('TGFB_Lymph',
                        'TGFB_Mac'))

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  geom_hline(yintercept=neg_line1, linetype="dashed") +
  geom_hline(yintercept=neg_line2, linetype="dashed") +
  theme_bw() + 
  stat_compare_means(label = "p.format", method= "wilcox.test",comparisons = my_comparisons) +
  theme(axis.ticks.x=element_blank()) + 
  labs(y="Dots/area")

# Plot data

#per cell
plot_data<-ish_data_perROI_zone_percell[,c(1,6,7,8,9,10)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG_Lymph',
                                                                                        'IFNG_Mac',
                                                                                        'TGFB_Lymph',
                                                                                        'TGFB_Mac'))

data_summary <- plot_data_melt[,c(3,4)] %>%
  group_by(variable) %>%
  summarize_all(list(median=median))

neg_line1<-median(ish_data_perROI_zone_percell$DapB_Lymph, na.rm = T)
neg_line2<-median(ish_data_perROI_zone_percell$DapB_Mac, na.rm = T)

my_comparisons = list(c('IFNG_Lymph',
                        'IFNG_Mac'),
                      c('TGFB_Lymph',
                        'TGFB_Mac'))

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  geom_hline(yintercept=neg_line1, linetype="dashed") +
  geom_hline(yintercept=neg_line2, linetype="dashed") +
  theme_bw() + 
  theme(legend.position='none') +
  stat_compare_means(label = "p.format", method= "wilcox.test",comparisons = my_comparisons) +
  theme(axis.ticks.x=element_blank()) + 
  labs(y="Dots/cell")

# version without 72215 #

#per cell
plot_data<-ish_data_perROI_zone_percell[!ish_data_perROI_zone_percell$Patient == 72215, c(1,6,7,8,9,10)]
plot_data_melt <- melt(plot_data, id.vars = c('Patient.ROI','Patient'),measure.vars = c('IFNG_Lymph',
                                                                                        'IFNG_Mac',
                                                                                        'TGFB_Lymph',
                                                                                        'TGFB_Mac'))

data_summary <- plot_data_melt[,c(3,4)] %>%
  group_by(variable) %>%
  summarize_all(list(median=median))

neg_line1<-median(ish_data_perROI_zone_percell$DapB_Lymph, na.rm = T)
neg_line2<-median(ish_data_perROI_zone_percell$DapB_Mac, na.rm = T)

my_comparisons = list(c('IFNG_Lymph',
                        'IFNG_Mac'),
                      c('TGFB_Lymph',
                        'TGFB_Mac'))

ggplot(data=plot_data_melt, aes(x=variable, y=value)) +
  geom_point(data=plot_data_melt, aes(x=variable, y=value, shape=as.factor(Patient)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  geom_errorbar(aes(y = median, ymin = median, ymax = median),
                color = "black", width = 0.2, data = data_summary) +
  geom_hline(yintercept=neg_line1, linetype="dashed") +
  geom_hline(yintercept=neg_line2, linetype="dashed") +
  theme_bw() + 
  theme(legend.position='none') +
  stat_compare_means(label = "p.format", method= "wilcox.test",comparisons = my_comparisons) +
  theme(axis.ticks.x=element_blank()) + 
  labs(y="Dots/cell")

# Get count of TGFB and IFNG expressing cells per ROI

# Summarise total dots and total cells per region

ish_data_pos_cells<-droplevels(ish_data_nonnec[ish_data_nonnec$Subcellular..DAB..Num.spots.estimated > 1, ])
ish_summary_pos_cells<-as.data.frame(table(ish_data_pos_cells$Image))
colnames(ish_summary_pos_cells)<-c('image','total_pos_cells')
ish_summary_pos_cells$Probe<-str_sub(ish_summary_pos_cells$image,-8,-5)

# add region area data

ish_summary_pos_cells$Area.µm.2<-0
ish_summary_pos_cells[ish_summary_pos_cells$image %in% ish_region_imm$Image,]$Area.µm.2<-ish_region_imm$Area.µm.2
ish_summary_pos_cells$Area.mm.2<-ish_summary_pos_cells$Area.µm.2 / 1000

# For TGFB and IFNG plot the # of positive cells per region

plot_data<-ish_summary_pos_cells[ish_summary_pos_cells$Probe %in% c("TGFB",'IFNG'),]
plot_data$Patient<-c('2384-1','2384-1','2384-2','2384-2','2384-3','2384-3','51210-1','51210-1',
                     '51210-2','51210-2','51210-3','51210-3','72215-1','72215-1','72215-2','72215-2',
                     '72215-3','72215-3','72215-4','72215-4','72215-5','72215-5')
plot_data_melt <- melt(plot_data, id.vars = c('Patient','Probe'), measure.vars = c('total_pos_cells'))

ggplot(plot_data_melt, aes(fill=Probe, y=value, x=Patient)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_grey() +
  theme_bw() +
  labs(x = 'ROI') + 
  labs(y = 'Count of Pos Cells')

# For TGFB and IFNG plot the # of positive cells per region per 500 um^2

ish_summary_pos_cells$scale_factor<-ish_summary_pos_cells$Area.µm.2 / (500*500)
ish_summary_pos_cells$scaled_cells<-ish_summary_pos_cells$total_pos_cells / ish_summary_pos_cells$scale_factor

plot_data<-ish_summary_pos_cells[ish_summary_pos_cells$Probe %in% c("TGFB",'IFNG'),]
plot_data$Patient<-c('2384-1','2384-1','2384-2','2384-2','2384-3','2384-3','51210-1','51210-1',
                     '51210-2','51210-2','51210-3','51210-3','72215-1','72215-1','72215-2','72215-2',
                     '72215-3','72215-3','72215-4','72215-4','72215-5','72215-5')
plot_data_melt <- melt(plot_data, id.vars = c('Patient','Probe'), measure.vars = c('scaled_cells'))

ggplot(plot_data_melt, aes(fill=Probe, y=value, x=Patient)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_grey() +
  theme_bw() +
  labs(x = 'ROI') + 
  labs(y = 'Count of Pos Cells')





