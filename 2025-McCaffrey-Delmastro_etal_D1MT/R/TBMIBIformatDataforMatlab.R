# TBMIBIformatDataforMatlab.R
# Author: Erin McCaffrey
# Date created: 190731
# Overview: Script reads in the normalized expression matrix for all cohort samples (myco and sarcoid). 
# It then assigns numerical dummy variables for the FlowSOM clusters and the major lineage type
# so that data can be more readily indexed with Matlab's obnoxious inability to handle
# mixed data-type matrices that aren't a royal pain in the ass to subset with normal indexing 

##..Read in annotated data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/dataPerCell")
cell_data<-read.csv('NHP_cohort_data_norm_annotated.csv')

##..Dummy variables for FlowSOM clusters..##

cell_data <- cell_data %>% mutate(FlowSOM_ID=case_when(cell_data$cell_type == "B_cell" ~ 1,
                                                   cell_data$cell_type == "CD4_T" ~ 2,
                                                   cell_data$cell_type == "CD8_T" ~ 3, 
                                                   cell_data$cell_type == "endothelial" ~ 4,
                                                   cell_data$cell_type == "fibroblast" ~ 5,
                                                   cell_data$cell_type == "giant_cell" ~ 6,
                                                   cell_data$cell_type == "imm_other" ~ 7,
                                                   cell_data$cell_type == "Mac_Mono" ~ 8,
                                                   cell_data$cell_type == "mast" ~ 9,
                                                   cell_data$cell_type == "neutrophil" ~ 10,
                                                   cell_data$cell_type == "T_other" ~ 11,
                                                   cell_data$cell_type == "Treg" ~ 12,
                                                   cell_data$cell_type == "epithelial" ~ 13))

##..Dummy variables for major cell lineages..##

# These are already encoded with lintype_num as 1 (myeloid), 2 (lymph), 3 (granulocyte), 4 (other immune), 5 (non-immune)

##..Save as csv..##

write.csv(cell_data, "NHP_cohort_data_norm_annotated_matlab.csv", row.names = FALSE)



