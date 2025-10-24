# MIBI_spatial_niche_network.R
# Created by: Erin McCaffrey
# Date created: 03/20/24
#
# Overview: This script takes the niche output from Jolene and visualizes the data 
# as a network with: 1) one node per cell type 2) the size of the cell type
# indicating the number of niches it appears in 3) the edge weight between cell
# types based on the number of times those cells co-occur in a niche together

library(purrr)
library(ggplot2)
library(stringr)
library(ggnet)
library(igraph)
library(circlize)
library(dplyr)

##..Step 1: Read in the data..##

data<-read.csv('./spatial_analysis/QUICHE/v2/tb_annotated_table_tb_binary_updated_all-sig.csv')
sc_data <- read.csv('cell_cohort_data_metabolic_zones.csv')

##..Step 2: Append the single cell metabolic data..##

niche_data_sc <- merge(data, sc_data, by = c('tiled_label','sample'))

##..Optionally subset based on top enriched in low versus high..##
data_summary <- data %>%
  group_by(new_label) %>%
  summarize(median = median(logFC))

high_burden_niches <- data_summary[data_summary$median > 0, ]$new_label
low_burden_niches <- data_summary[data_summary$median < 0, ]$new_label

##..Step 3: Get a list of the unique niches and cell types..##

# unique_niches <- unique(data$new_label)
# unique_niches <- low_burden_niches
unique_niches <- high_burden_niches
unique_niches <- gsub("_", "", unique_niches)
unique_niches <- gsub("\\+","",unique_niches)

unique_cell_types <- unique(niche_data_sc$pheno_corrected)
unique_cell_types <- gsub("_", "", unique_cell_types)
unique_cell_types <- gsub("\\+","",unique_cell_types)
unique_cell_types <- unique_cell_types[!unique_cell_types == 'giantcell']

##..Step 4: Generate network data..##

cell_counts <- integer(length(unique_cell_types))
niche_network <- matrix(, nrow = length(unique_cell_types), ncol = length(unique_cell_types))

for(i in 1:length(unique_cell_types)) {
  
  cell_type <- unique_cell_types[i]
  
  # count instances of the cell type
  cell_type_count <- sum(str_count(unique_niches,coll(cell_type)))
  cell_counts[i] <- cell_type_count
  
  # subset to include niches with this cell type
  cell_type_niches <- unique_niches[grep(cell_type, unique_niches)]
  
  # count the number of occurrences with other cell types
  for(j in 1:length(unique_cell_types)) {
    
    cell_type_2 <- unique_cell_types[j]
    
    if(cell_type == cell_type_2) {
      niche_network[i, j] <- 0
    } else {
      cell_pair_niches <- cell_type_niches[grep(cell_type_2, cell_type_niches)]
      cell_pair_counts <- length(cell_pair_niches)
      niche_network[i, j] <- cell_pair_counts
    }
  }
}

cell_count_data <- as.data.frame(unique_cell_types)
cell_count_data$count <- cell_counts
colnames(cell_count_data) <- c('pheno_corrected', 'count')
cell_count_data$norm_count <- cell_count_data$count / max(cell_count_data$count)

rownames(niche_network) <- unique_cell_types
colnames(niche_network) <- unique_cell_types
niche_network

##..Step 5: Visualize network..##

g <- graph_from_adjacency_matrix(niche_network/max(niche_network),
                                 mode="undirected",
                                 weighted=TRUE,
                                 diag=FALSE
)

# set node size
node.size<-setNames(cell_count_data$norm_count*30, unique_cell_types)

# define layout
layout <- layout.fruchterman.reingold(g)

# set edge thickness
E(g)$width <- 5*E(g)$weight

# Remove vertices without edges
floating = which(degree(g)==0)
g2 = delete.vertices(g, floating)
l2 = layout[-floating,]
node.size2<-node.size[!(names(node.size) %in% names(floating))]

# get color key
color_key <- read.csv("./keys/cell_color_key_spatial.csv")
plot_populations<-names(node.size2)
plot_colors<-droplevels(color_key[color_key$Pheno %in% plot_populations,])
plot_colors$Pheno<-factor(plot_colors$Pheno, levels = plot_populations)
plot_colors<-plot_colors[order(plot_colors$Pheno),]
color<-as.vector(plot_colors$Hex)

##..Visualize network..##

plot(g2, 
     rescale = T,
     vertex.size=node.size2,
     layout =  l2,
     vertex.color = color,
     vertex.label.family="Helvetica",
     edge.color = "dimgrey",
     vertex.label.cex=0.5,
     vertex.label.color='black')