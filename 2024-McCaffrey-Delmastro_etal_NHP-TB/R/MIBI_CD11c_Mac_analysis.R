# MIBI_CD11c_Mac_analysis.R
# Date created: 12/21/2023
# Author: Erin McCaffrey
#  
# Overview: Based on the correlation between CFU and the abudnace of CD11c+ Macrophages, 
# here we perform some additional analysis of the macrophage populations found
# in the TB granulomas

library(forcats)
library(viridis)
library(dplyr)
library(stringr)
library(ggpubr)
devtools::install_github("psyteachr/introdataviz")
library(introdataviz)

##..Step 1: Read in data ##

data<-read.csv('cell_stats_all_samples_meta_data.csv')
data<-droplevels(data[data$category == 'pheno_of_total',])
data<-tibble::rowid_to_column(data, "ID")
animal_color_key <- read.csv("./keys/animal_color_key.csv")

##..Step 2: Subset to only include macrophages/monocytes..##

macs <- c("CD11c+_Mac","CD14+CD11c+_Mac","CD14+_Mac_Mono","CD163+_Mac","CD206+_Mac","CD68+_Mac",
          "FN1+_Mac",'giant_cell')
data_mac <- droplevels(data[data$variable %in% macs,])

##..Step 3: Separate count and density data..##

count_data<-reshape2::dcast(data_mac, sample + log_CFU + burden + Animal_Code ~ variable, value.var = "n", fun.aggregate = sum)
density_data<-reshape2::dcast(data_mac, sample + log_CFU + burden + Animal_Code ~ variable, value.var = "cell_density", fun.aggregate = sum)

##..Step 4: Append total macrophages..##

count_data <- count_data %>%
  mutate(total_macs = rowSums(select(.,6:13)))

freq_data <- count_data %>% 
  mutate_at(vars(6:13),function(i)i/.$total_macs)

##..Step 5: Melt for easy plotting..##

freq_data.m <- reshape2::melt(freq_data, id.vars = c('sample', 'burden', 'log_CFU','Animal_Code'))
density_data.m <- reshape2::melt(density_data, id.vars = c('sample', 'burden', 'log_CFU','Animal_Code'))

##..Step 6: Plot..##

plot_data <- freq_data.m[!freq_data.m$variable %in% c('total_macs','gran_CFU'),]

ggplot(data = plot_data, aes(x = fct_reorder(as.factor(variable), value,.fun=median,.desc=TRUE), 
                             value, fill=burden)) +
  geom_boxplot(width = 0.75, alpha = .6, fatten = NULL, show.legend = TRUE, outliers = FALSE) +
  scale_color_manual(values = color) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.75)) +
  stat_compare_means(aes(group = burden), method= "wilcox.test", label = "p.format") +
  geom_jitter(aes(color = Animal_Code), width = 0.2, size = 2, show.legend = TRUE) +
  theme_bw() +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total Macrophages') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) 

## Plot a version faceted by subset with the animal data points shown ##
plot_animals<-levels(factor(plot_data$Animal_Code))
plot_colors<-droplevels(animal_color_key[animal_color_key$Animal_Code %in% plot_animals,])
plot_colors$Animal<-factor(plot_colors$Animal_Code, levels = plot_animals)
plot_colors<-plot_colors[order(plot_colors$Animal_Code),]
color<-as.vector(plot_colors$colour)

ggplot(data = plot_data, aes(x = burden, y = value)) +
  geom_boxplot(show.legend = TRUE, outliers = FALSE) +
  geom_jitter(aes(color = Animal_Code), width = 0.2, size = 2, show.legend = TRUE) +
  scale_color_manual(values = color) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Cell Type') + 
  labs(y = 'Frequency of Total Macrophages') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), 
        axis.text.x = element_text(angle=35,hjust=1)) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2)

##..Generate log2(FC) between the high versus low burden grans..##
summary_data <- freq_data[,!names(freq_data) %in% c('sample',
                                                    'log_CFU',
                                                    'total_macs')]
data_summary <- summary_data %>%
  group_by(burden) %>%
  summarize_all(list(mean=mean))

data_summary<-as.data.frame(t(data_summary[,-c(1)]))
names(data_summary)<-c('high','low')
data_summary$pheno_corrected <- rownames(data_summary)
data_summary$pheno_corrected <- stringr::str_remove(data_summary$pheno_corrected, "_mean")
data_summary$high<-as.numeric(as.character(data_summary$high))
data_summary$low<-as.numeric(as.character(data_summary$low))
data_summary$FC<-log2(data_summary$high/data_summary$low)

color_key <- read.csv("./keys/cell_color_key.csv")
plot_populations<-levels(factor(data_summary$pheno_corrected))
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

ggplot(data=data_summary[is.finite(data_summary$FC),], 
       aes(x=reorder(pheno_corrected, -FC), y=FC, fill=pheno_corrected)) +
  geom_bar(stat="Identity", width = 1) +
  scale_fill_manual(values = color) +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_minimal() + 
  theme(legend.position = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) + 
  labs(y="log2(mean freq high / mean freq low)") 


