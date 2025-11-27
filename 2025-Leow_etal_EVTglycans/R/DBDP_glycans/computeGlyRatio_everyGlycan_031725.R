# Get ratio of each glycan
# Author: Ke Leow
# Date: 03/17/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)

#--------------------------------
# Load MALDI data
#--------------------------------
#data before scaling/batch correction
maldi_data = read_csv("data/DBDP_transition/MALDI/library_matched/DBDP_transition_mean_gly_permask_wMeta.csv") 

#--------------------------------
# Calculate ratio of each glycan
#--------------------------------
data_ratio <- maldi_data %>% 
  gather(key = "glycans", value = "value", H2N2F1:H10N9F1) %>% 
  group_by(mask) %>% 
  mutate(total_gly_int = sum(value)) %>% 
  mutate(value = value/total_gly_int) %>% 
  mutate(glycans = str_c("R_", glycans)) %>% 
  select(mask, glycans,value) %>% 
  pivot_wider(names_from = glycans, values_from = value) %>% 
  inner_join(.,maldi_data, by = join_by(mask)) %>% 
  select(mask, Region, everything())


write.csv(data_ratio, 'data/DBDP_transition/MALDI/library_matched/DBDP_transition_glyRatio_everyGlycan_031725.csv', row.names = FALSE)
