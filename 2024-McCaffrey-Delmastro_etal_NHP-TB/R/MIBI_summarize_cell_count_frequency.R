# MIBI_summarize_cell_count_frequency.R
# Author: Erin McCaffrey
# 
# Overview: Script takes in the single cell data and returns a table summarizing
# the counts and frequencies of cell proportions per specimen. 
# 
# Packages needed: dplyr
# Input data: NHP_TB_MIBI-TOF_singlecelldata.csv

##..load necessary packages..##
library(dplyr)

##..read in data input..##
setwd("/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My Drive/Grad_School/AngeloLab/MIBIProjects/NHP_TB/Manuscript/public_data")
cell_table<-read.csv('NHP_TB_MIBI-TOF_singlecelldata.csv')


##..Summarize all cell phenotypes..##
cell_summary1 <- cell_table %>% 
  group_by(sample, .drop = F ) %>%
  count(pheno_corrected) %>%
  mutate(freq_of_total = prop.table(n)) %>%
  rename('variable' = 'pheno_corrected') %>%
  mutate(category = 'pheno_of_total') %>%
  bind_rows(.,
            cell_table %>%
              group_by(sample, sublineage) %>%
              tally() %>%
              mutate(freq_of_total = prop.table(n)) %>%
              rename('variable' = 'sublineage') %>%
              mutate(category = 'sublineage_of_total'),
            cell_table %>%
              group_by(sample, majorlineage) %>%
              tally() %>%
              mutate(freq_of_total = prop.table(n)) %>%
              rename('variable' = 'majorlineage') %>%
              mutate(category = 'majorlineage_of_total'))

##..Summarize immune cells only..##
cell_table_immune<-droplevels(cell_table[cell_table$majorlineage == 'immune',])

cell_summary2<-cell_table_immune %>% 
  group_by(sample) %>%
  count(pheno_corrected, .drop = F) %>%
  mutate(freq_of_total = prop.table(n)) %>%
  rename('variable' = 'pheno_corrected') %>%
  mutate(category = 'pheno_of_totalimmune') %>%
  bind_rows(.,
          cell_table_immune %>%
            group_by(sample, sublineage) %>%
            tally() %>%
            mutate(freq_of_total = prop.table(n)) %>%
            rename('variable' = 'sublineage') %>%
            mutate(category = 'sublineage_of_totalimmune'))

##..Combine summaries..##
cell_summary<-rbind(cell_summary1,cell_summary2)

##..Add in cellular area to compute density..##

area_data<-read.csv('./masks/all_samples_mask_data.csv')
cell_summary_area<-merge(cell_summary, area_data, by = c("sample"))
cell_summary_area$um2<-cell_summary_area$cellular_px*((400*400)/(1024*1024))
cell_summary_area$cell_density<-cell_summary_area$n/cell_summary_area$um2
cell_summary_area$freq_density<-cell_summary_area$freq_of_total/cell_summary_area$um2

##..Add necrosis:viable ratio..##
cell_summary_area$necrosis_viable_ratio<-cell_summary_area$necrosis_px/cell_summary_area$cellular_px

##..Add in metadata..##
meta_data<-read.csv('study_cohort_metadata.csv')
# drop<-c('sample1','sample2','sample19','sample30','sample31','sample36',
#         'sample39','sample45-2','sample46-2','sample46-3')
# cell_summary_area<-droplevels(cell_summary_area[!cell_summary_area$sample %in% drop,])
# cell_summary_area[cell_summary_area$sample == 'sample46-1',]$sample <- 'sample46'
cell_summary_meta_data<-merge(cell_summary_area, meta_data, by = "sample", all=TRUE)

##..Save..##
write.csv(cell_summary_meta_data, 'cell_stats_all_samples_meta_data.csv', row.names = F)



