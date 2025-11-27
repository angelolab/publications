# Get ratio of each glycan and glycan types
# Author: Ke Leow
# Date: 03/18/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)

#--------------------------------
# Load MALDI data
#--------------------------------
#data before scaling/batch correction
maldi_data = read_csv("data/DBDP_transition/MALDI/library_matched/DBDP_transition_glyRatio_everyGlycan_031725.csv") 

glycans = maldi_data %>% 
  select(H2N2F1:H10N9F1) %>% colnames()
data_int <- maldi_data %>% select(R_H2N2F1:R_H10N9F1) %>% as.data.frame()
#--------------------------------
# Annotate glycans
#--------------------------------
my_gly_col <- tibble(composition = glycans) %>% 
  mutate(fucosylated = ifelse(grepl("F", composition), 1, 0)) %>% 
  mutate(sialylated = ifelse(grepl("S", composition), 1, 0)) %>% 
  mutate(fucosylated_sialylated = ifelse(grepl("S", composition) & grepl("F", composition), 1, 0)) %>% 
  mutate(highMannose = ifelse(grepl("H5N2\\b|H6N2\\b|H7N2\\b|H8N2\\b|H9N2\\b", composition), 1, 0)) %>% 
  mutate(hybrid = ifelse(grepl("H5N3|H6N3|H6N4", composition), 1, 0)) %>% 
  mutate(paucimannose = ifelse(grepl("H4N2|H3N2|H2N2", composition), 1, 0)) %>% 
  mutate(agalactosylated = ifelse(grepl("H3N3|H3N4|H3N5|H3N6", composition), 1, 0)) %>% 
  mutate(biantennary = ifelse(grepl("H5N4", composition), 1, 0)) %>% 
  mutate(triantennary = ifelse(grepl("H6N5", composition), 1, 0)) %>% 
  mutate(tetraantennary = ifelse(grepl("H7N6", composition), 1, 0)) %>% 
  mutate(highlybranched = ifelse(grepl("H7N6|H8N7|H9N8|H10N9", composition), 1, 0)) %>% 
  mutate(polylacnac = ifelse(grepl("H8N7|H9N8|H10N9", composition), 1, 0))%>% 
  mutate(bisecting = ifelse(grepl("H3N5|H4N5|H5N5|H3N6|H4N6|H5N6|H6N6|H3N7|H4N7|H5N7|H6N7|H7N7|H8N8", composition), 1, 0))%>%
  mutate(biant_partialgal = ifelse(grepl("H3N3|H3N4|H4N4", composition), 1, 0))%>% 
  mutate(triant_partialgal = ifelse(grepl("H3N5|H4N5|H5N5", composition), 1, 0))%>% 
  mutate(tetraant_partialgal = ifelse(grepl("H3N6|H4N6|H5N6|H6N6", composition), 1, 0))%>% 
  mutate(polylacnac_partialgal = ifelse(grepl("H4N7|H5N7|H6N7|H7N7|H8N8", composition), 1, 0))%>% 
  mutate(biantennary_fuc = ifelse(grepl("H5N4", composition) & grepl("F", composition), 1, 0)) %>% 
  mutate(triantennary_fuc = ifelse(grepl("H6N5", composition)& grepl("F", composition), 1, 0)) %>% 
  mutate(tetraantennary_fuc = ifelse(grepl("H7N6", composition)& grepl("F", composition), 1, 0)) %>% 
  mutate(highlybranched_fuc = ifelse(grepl("H7N6|H8N7|H9N8|H10N9", composition) & grepl("F", composition), 1, 0)) %>% 
  mutate(polylacnac_fuc = ifelse(grepl("H8N7|H9N8|H10N9", composition) & grepl("F", composition), 1, 0))%>% 
  mutate(bisecting_fuc = ifelse(grepl("H3N5|H4N5|H5N5|H3N6|H4N6|H5N6|H6N6|H3N7|H4N7|H5N7|H6N7|H7N7|H8N8", composition)& grepl("F", composition), 1, 0))%>%
  mutate(biant_partialgal_fuc = ifelse(grepl("H3N3|H3N4|H4N4", composition)& grepl("F", composition), 1, 0))%>% 
  mutate(triant_partialgal_fuc = ifelse(grepl("H3N5|H4N5|H5N5", composition)& grepl("F", composition), 1, 0))%>% 
  mutate(tetraant_partialgal_fuc = ifelse(grepl("H3N6|H4N6|H5N6|H6N6", composition)& grepl("F", composition), 1, 0))%>% 
  mutate(polylacnac_partialgal_fuc = ifelse(grepl("H4N7|H5N7|H6N7|H7N7|H8N8", composition)& grepl("F", composition), 1, 0))%>% 
  mutate(biantennary_sia = ifelse(grepl("H5N4", composition) & grepl("S", composition), 1, 0)) %>% 
  mutate(triantennary_sia = ifelse(grepl("H6N5", composition)& grepl("S", composition), 1, 0)) %>% 
  mutate(tetraantennary_sia = ifelse(grepl("H7N6", composition)& grepl("S", composition), 1, 0)) %>% 
  mutate(highlybranched_sia = ifelse(grepl("H7N6|H8N7|H9N8|H10N9", composition) & grepl("S", composition), 1, 0)) %>% 
  mutate(polylacnac_sia = ifelse(grepl("H8N7|H9N8|H10N9", composition) & grepl("S", composition), 1, 0))%>% 
  mutate(bisecting_sia = ifelse(grepl("H3N5|H4N5|H5N5|H3N6|H4N6|H5N6|H6N6|H3N7|H4N7|H5N7|H6N7|H7N7|H8N8", composition)& grepl("S", composition), 1, 0))%>%
  mutate(biant_partialgal_sia = ifelse(grepl("H3N3|H3N4|H4N4", composition)& grepl("S", composition), 1, 0))%>% 
  mutate(triant_partialgal_sia = ifelse(grepl("H3N5|H4N5|H5N5", composition)& grepl("S", composition), 1, 0))%>% 
  mutate(tetraant_partialgal_sia = ifelse(grepl("H3N6|H4N6|H5N6|H6N6", composition)& grepl("S", composition), 1, 0))%>% 
  mutate(polylacnac_partialgal_sia = ifelse(grepl("H4N7|H5N7|H6N7|H7N7|H8N8", composition)& grepl("S", composition), 1, 0))%>% 
  mutate(agalactosylated_fuc = ifelse(grepl("H3N3|H3N4|H3N5|H3N6", composition)& grepl("F", composition), 1, 0)) %>% 
  data.frame()
rownames(my_gly_col) = glycans

#--------------------------------
# Calculate ratio/total int of glycan types - take sum of glycans in each type
#--------------------------------
# Define glycan types (you can add more columns if needed)
glycan_types <- colnames(my_gly_col)[-1]

# Initialize a new dataframe to store summed intensities for each glycan type for each sample
summed_intensities <- maldi_data %>%
  select(1)  # Selecting sample IDs (assuming sample IDs are in the first column)

# Loop over each glycan type
for (glycan_type in glycan_types) {
  # Identify glycans belonging to the current glycan type (binary annotation 1)
  glycan_columns <- my_gly_col[glycan_type] == 1
  glycan_data_columns <- data_int[, glycan_columns]
  # Sum the intensities of these glycans for each sample
  summed_intensities[[paste0("R_", glycan_type)]] <- rowSums(glycan_data_columns)
}

maldi_data_ratio = summed_intensities %>% 
  inner_join(maldi_data) 

write.csv(maldi_data_ratio, 'data/DBDP_transition/MALDI/library_matched/DBDP_transition_glyRatio_glycanTypes_031825.csv', row.names = FALSE)
