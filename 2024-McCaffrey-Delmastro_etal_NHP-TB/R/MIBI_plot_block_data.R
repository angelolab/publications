# MIBI_plot_block_data.R
# Author: Erin McCaffrey 
# Date created: 210609
#
# Overview: This script plots the meta data for the cohort

require(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
require(ggbeeswarm)
library(psych)
library(gtools)
library(pals)
library(ggpubr)
library(ggstatsplot)

##..Step 1: Import data..##

data_gran <- read.csv("./cohort_metadata/study_cohort_metadata.csv")
data_animal <- read.csv("./cohort_metadata/study_cohort_animal_metadata.csv")
animal_color_key <- read.csv("./keys/animal_color_key.csv")

##..Step 2: Append granuloma age data..##

data_animal$Weeks_Necropsy <- data_animal$days_to_necropsy/7
data_gran$Week_Necropsy <- data_animal$Weeks_Necropsy[match(data_gran$Animal,data_animal$Animal)]
data_gran <- data_gran %>% relocate(Week_Necropsy, .before=Granuloma)
data_gran$gran_age <- data_gran$Week_Necropsy - data_gran$Week_Detected
data_gran <- data_gran %>% relocate(gran_age, .before=Granuloma)

##..Step 3: Add a column to reflect high v low burden from median cutoff..##

meta_data <- data_gran
thresh <- median(meta_data$log_CFU)            
meta_data$burden <- 'high'
meta_data$burden[meta_data$log_CFU < thresh] <- 'low'
meta_data <- meta_data %>% relocate(burden, .before=CT_size)

##..Step 4: Melt..##

measure_vars <- c("gran_CFU", "CT_size", "FDG_SUV", "Week_Detected", "gran_age", 
                 "Week_Necropsy", "Granuloma", "Multifocal", "log_CFU",
                 "Necrotic", "Non.necrotic", "Fibrinoid_debri", "Fibrosis", 
                 "Mineralization", "Collagenization", "Neutrophillic", "Early_Evolving")
data_gran.m<-reshape2::melt(data_gran, id.vars = c('Animal_Code','sample'), measure.vars = measure_vars )

##..Step 5: Plot data..##

# plot number of granulomas per animal in descending order..##
n_gran_data<-as.data.frame(table(data_gran$Animal))
colnames(n_gran_data)<-c('Animal','n_gran')
median_CFU<-data_gran %>%
  group_by(Animal) %>%
  summarise(med_CFU = median(gran_CFU))
median_CFU$Animal<-as.factor(median_CFU$Animal)
n_gran_data_CFU<-left_join(n_gran_data, median_CFU)
n_gran_data_CFU$log_med_CFU<-log10(as.numeric(n_gran_data_CFU$med_CFU)+1)

ggplot(n_gran_data_CFU, aes(x=reorder(Animal, -n_gran), y=n_gran, fill=log_med_CFU)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_continuous(low="blue", high="red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  labs(y="# Specimens", x='Animal') 


# plot animal versus CFU with one point per gran in ascending mean order
plot_data<-data_gran
plot_data$gran_CFU<-as.numeric(plot_data$gran_CFU) + 1
CFU_plot<-ggplot(plot_data, aes(x=reorder(as.factor(Animal_Code), -gran_CFU), 
                      y=as.numeric(gran_CFU), color = as.factor(Animal))) + 
  geom_jitter(aes(shape=as.factor(Necrotic)), width=0.25, size = 5) +
  scale_color_manual(values = as.vector(glasbey(16))) +
  theme_bw() +
  coord_flip() +
  scale_y_continuous(trans='log10') +
  labs(y="CFU/granuloma", x="Animal") +
  theme(legend.position = "None") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CFU_plot

# plot other granuloma data:

# CT size
plot_data<-na.omit(data_gran.m[data_gran.m$variable=="CT_size",])

plot_animals<-levels(factor(plot_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(plot_data, aes(x=variable, y=as.numeric(value), color = as.factor(Animal_Code))) + 
  geom_quasirandom(width=0.1, size=2) +
  stat_summary(fun=median, geom="crossbar", width=0.3, color="black") +
  ylim(0, 8) +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(y="size (mm)", x="n=49 granulomas") +
  theme(panel.grid=element_blank(), axis.text.x=element_blank())

# FDG
plot_data<-na.omit(data_gran.m[data_gran.m$variable=="FDG_SUV",])

plot_animals<-levels(factor(plot_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(plot_data, aes(x=variable, y=as.numeric(value), color = as.factor(Animal_Code))) + 
  geom_quasirandom(width=0.1, size=2) +
  stat_summary(fun=median, geom="crossbar", width=0.3, color="black") +
  ylim(0, 29) +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(y="SUVR", x="n=41 granulomas") +
  theme(panel.grid=element_blank(), axis.text.x=element_blank())

# CFU/gran
plot_data<-na.omit(data_gran.m[data_gran.m$variable=="log_CFU",])

plot_animals<-levels(factor(plot_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(plot_data, aes(x=variable, y=as.numeric(value), color = as.factor(Animal_Code))) +
  geom_quasirandom(width=0.1, size=2) +
  stat_summary(fun=median, geom="crossbar", width=0.3, color="black") +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(y="CFU", x="n=52 granulomas") +
  theme(panel.grid=element_blank(), axis.text.x=element_blank())

# Granuloma Age
plot_data<-na.omit(data_gran.m[data_gran.m$variable=="gran_age",])

plot_animals<-levels(factor(plot_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(plot_data, aes(x=variable, y=as.numeric(value), color = as.factor(Animal_Code))) + 
  geom_quasirandom(width=0.15, size=2) +
  stat_summary(fun=median, geom="crossbar", width=0.3, color="black") +
  theme_bw() +
  ylim(0, 9) +
  scale_color_manual(values = color) +
  labs(y="Granuloma Age", x="n=51 granulomas") +
  theme(panel.grid=element_blank(), axis.text.x=element_blank())

##..Step 6: Evaluate at CFU correlations with features..##

# run individual correlations
corr.test(data_gran$gran_CFU, data_gran$gran_age, method = 'spearman')

# individual correlation plot
data_gran$log_CFU<-log10(as.numeric(data_gran$gran_CFU)+1)
data_gran <- data_gran %>% relocate(log_CFU, .before=CT_size)

ggscatterstats(data = data_gran, 
               x = log_CFU,
               y = CT_size,
               type = "spearman",
               xlab = "log(CFU+1)", 
               ylab = "Age (weeks)")

# plot all bi-axial relationships
corr_vars <- c("CT_size", "FDG_SUV", "gran_age")
corr_data <- melt(data_gran, id.vars = c('Animal_Code','sample','log_CFU'), measure.vars = corr_vars)

plot_animals<-levels(factor(corr_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(corr_data, aes(x=log_CFU, y=as.numeric(value)), color = as.factor(Animal_Code)) + 
  geom_point(aes(color = as.factor(Animal_Code))) + 
  geom_smooth(method='lm', formula= y~x) +
  scale_color_manual(values = color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = "None") +
  facet_wrap(.~variable, scale='free_y')

##..Step 7: Evaluate CFU differences based on histological status..##

plot_animals<-levels(factor(data_gran$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(data_gran, aes(x=as.factor(Necrosis_Score), y=CT_size)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(color = Animal_Code), size=1, alpha=1, width = 0.25)  +
  scale_color_manual(values = color) +
  stat_anova_test(aes(group = as.factor(Necrosis_Score)),  method = "one_way") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  theme(legend.position = "None") +
  labs(y="log_CFU", x='Necrosis Score') 

# assessing normality to decide if anova is appropriate for CFU analysis
ggqqplot(data_gran$log_CFU)
ggdensity(data_gran$log_CFU)
shapiro.test(data_gran$log_CFU)

anova <- aov(log_CFU ~ as.factor(Necrosis_Score), data = data_gran)
summary(anova)

##..Step 8: Plot granuloma data per animal..##

# append the animal code
data_animal <- left_join(data_animal, 
                         animal_color_key[,names(animal_color_key) %in% c('Animal','Animal_Code')],
                         by = c('Animal'))

# transform CFU to log10 CFU
data_animal$thoracic_CFU <- 1 + (log10(data_animal$thoracic_CFU))
data_animal$lung_CFU <- 1 + (log10(data_animal$lung_CFU))
data_animal$LN_CFU <- 1 + (log10(data_animal$LN_CFU))

# melt and subset
data_animal.m<-melt(data_animal, id.vars = c("Animal","Animal_Code"))
plot_data<-data_animal.m[data_animal.m$variable %in% c('Mtb_dose','Weeks_Necropsy','necropsy_score',
                                                       'thoracic_CFU','lung_CFU','LN_CFU'),]
plot_data$variable<-factor(plot_data$variable, levels = c('Mtb_dose','Weeks_Necropsy','necropsy_score',
                                                          'thoracic_CFU','lung_CFU','LN_CFU'))
plot_data<-plot_data[order(plot_data$variable),]

plot_animals<-levels(factor(plot_data$Animal_Code))
plot_animals<-mixedsort(plot_animals)
plot_data$Animal_Code<-factor(plot_data$Animal_Code, levels = plot_animals)
plot_data<-plot_data[order(plot_data$Animal_Code),]

plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal_Code<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(plot_data, aes(x=as.factor(Animal_Code), y=as.numeric(value)), color = as.factor(Animal_Code)) + 
  geom_col(aes(fill = as.factor(Animal_Code))) + 
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.ticks.x=element_blank()) +
  facet_wrap(.~variable, scale='free_y')

##..Step 9: Export updated metadata..##
write.csv(data_gran, 'study_cohort_metadata.csv', row.names = F) 
