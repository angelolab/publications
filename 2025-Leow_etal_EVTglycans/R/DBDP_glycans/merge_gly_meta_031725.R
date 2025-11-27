#--------------------------------
# merge avg dat with meta
# Author: Ke Leow
# Date: 03/17/25
#--------------------------------
library(tidyverse)

#metadata
meta <- read_csv("data/DBDP_transition/MALDI/library_matched/fiji/Overlay Elements of H5N4F1_annotated.csv")

#data: average intensities from extracted TIFFs
dat <- read_csv("data/DBDP_transition/MALDI/library_matched/DBDP_transition_mean_gly_permask.csv") 

dat_meta <- inner_join(dat[2:ncol(dat)], meta %>% select(Name, Region), by = c('mask' = 'Name')) 

write.csv(dat_meta, "data/DBDP_transition/MALDI/library_matched/DBDP_transition_mean_gly_permask_wMeta.csv", row.names = FALSE)

  