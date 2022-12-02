# PAH_FunctionalMarkerVesselScore.R
# Author: Erin McCaffrey 
# Date created: 200628
# Overview: This script reads in vessel scoring data and percent positive immune cells per 
# marker per point. Next it produces plots comparing these frequencies across vessel scores
# for all functional markers. It also breaks this down by certain cell subsets and markers.

require(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(psych)
library(devtools)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Collaborations/PAH manuscript/Datasets")
vessel_data<-read.csv("VesselsSummaryImmunePerPoint.csv")
marker_freqs<-read.csv("allpoints_freq-of-positivecells-persubset.csv")

##..Append the vessel scoring to the frequency data..##

marker_freqs$vessel_score<-rep(vessel_data$Vessel_Scoring, 13)

##..Drop regions without a vessel score..##
marker_freqs_vessel<- marker_freqs %>% filter(!is.na(vessel_score)) 

##..Melt..##

marker_freqs_vessel.m<-melt(marker_freqs_vessel, id=c('Point_num','cell_lineage','vessel_score', 'Subgroup'))

##..Define functional markers..##

functional_markers<-c('CD45RO','GranzB','CD163','HLA.Class.I','HLA.DR','IDO1','IFNg','iNOS','KI67',
                          'MPO','PDL1b','SAMHD1','TIM.3','bCatenin','CD141','NaK.ATPase')

##..Plot markers for immune cells across stage..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='immune',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)
            

all<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
all

total_immune<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
total_immune$cell_type<-'total_immune'

##.NK cells..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='NK',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)


NK<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
NK

NK_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
NK_stats$cell_type<-'NK'

##.CD4 T..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Th',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)


Th<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Th

Th_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
Th_stats$cell_type<-'Th'

##.CD8 T..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Tc',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)


Tc<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Tc

Tc_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
Tc_stats$cell_type<-"Tc"

##.Macro..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Macro',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)


Macro<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Macro

Macro_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
Macro_stats$cell_type<-'Macro'


##.Mono..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Mono',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)


Mono<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Mono

Mono_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
Mono_stats$cell_type<-'Mono'

##..DC..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='DC',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)

DC<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
DC

DC_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
DC_stats$cell_type<-'DC'

##..Neutro..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Neutro',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)

Neutro<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Neutro

Neutro_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
Neutro_stats$cell_type<-'Neutro'

##..Bcell..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Bcell',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.factor(plot_data$vessel_score)

Bcell<-ggplot(data = plot_data, aes(x = as.factor(vessel_score), y = as.numeric(value))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Bcell

Bcell_stats<-compare_means(value ~ vessel_score, group.by = "variable", data = plot_data)
Bcell_stats$cell_type<-'Bcell'

##..Create combined stats file..##

all_stats<-rbind(total_immune,NK_stats,Th_stats,Tc_stats,Macro_stats,Mono_stats,DC_stats,Neutro_stats,Bcell_stats)
write.csv(all_stats, "markerfreqs-across-vessel-stage.csv",row.names = F)

##..Correlation between markers and score..##

##..Total immune..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='immune',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)

immune_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
immune_lin

##..NK..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='NK',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


NK_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
NK_lin

##..Tc..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Tc',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


Tc_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Tc_lin


##..Th..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Th',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


Th_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Th_lin

##..DC..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='DC',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


DC_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
DC_lin

##..Macro..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Macro',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


Macro_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Macro_lin


##..Mono..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Mono',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


Mono_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Mono_lin

##..Neutro..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Neutro',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


Neutro_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Neutro_lin

##..Bcell..##

plot_data<-droplevels(marker_freqs_vessel.m[marker_freqs_vessel.m$cell_lineage=='Bcell',])
plot_data<-droplevels(plot_data[plot_data$variable %in% functional_markers,])
plot_data$value<-as.numeric(plot_data$value)
plot_data$vessel_score<-as.numeric(plot_data$vessel_score)


Bcell_lin<-ggplot(data = plot_data, aes(x = vessel_score, y = value)) + 
  geom_smooth(method='lm') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 3) +
  geom_point(shape=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Score') + 
  labs(y = 'Frequency of Total Immune') +
  facet_wrap(.~variable, scales='free_y')
Bcell_lin

##..Produce correlation analysis of immune cell percent pos with vessel stage..##

immune_vessel_freqs<-droplevels(marker_freqs_vessel[marker_freqs_vessel$cell_lineage=='DC',])
immune_vessel_freqs<-droplevels(immune_vessel_freqs[,c(functional_markers,'vessel_score')])
immune_vessel_freqs$vessel_score<-as.numeric(immune_vessel_freqs$vessel_score)

#Function for flattening matrix from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

freq_corr<-corr.test(as.matrix(immune_vessel_freqs), y = NULL, use = "pairwise",method="pearson",adjust="fdr")
freq_corr_flat<-flattenCorrMatrix(freq_corr$r,freq_corr$p)
freq_corr_flat$sig<-freq_corr_flat$p<0.05

corrplot(freq_corr$r, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45, method = 'square')
corrplot(freq_corr$r, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45, 
         p.mat = freq_corr$p, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",
         method='square',  col = diverging_hcl(100, palette = 'Blue-Red 3'))
                                