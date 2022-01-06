# TBMIBInetworkVis.R
# Author: Erin McCaffrey 
# Date created: 191030
# Overview: This script reads in a pairwise spatial enrichment matrix with a header of labels for the rows and columns. 
# First it filters out the negative interactions (negative z) and puts them in a separate matrix. It uses igraph to
# turn the matrices into graph objects and visualizes as a network.

library(igraph)

##..Import data..##
data<-read.csv("data/pooled-enrich_TB.csv")

##..Convert data to labeled matrix..##

markerTitles<-as.vector(colnames(data))
data_mat<-as.matrix(data)
rownames(data_mat)<-markerTitles
colnames(data_mat)<-markerTitles
data_mat[is.na(data_mat)]<-0
# drop NaKATPase and HLA 1
data_mat<-data_mat[-35:-36, -35:-36]

# Create a graph adjacency based on correlation distances between genes in  pairwise fashion.

g <- graph_from_adjacency_matrix(data_mat,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)


# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

# # Remove edges less than 0 (avoidances)
g <- delete_edges(g, E(g)[which(E(g)$weight<0)])

# Remove vertices without edges
g<-delete.vertices(simplify(g), degree(g)==0)

# plot
plot(g, rescale = TRUE)

# get communities with leading eigenvector
eg<-cluster_leading_eigen(g)
eg.clustering <- make_clusters(g, membership = eg$membership)

set.seed(4)
plot(eg.clustering,g, vertex.label.color= "white", vertex.label.cex=0.5, 
     vertex.label.family="Helvetica", edge.curved=0.25, edge.color = "dimgrey", vertex.size=18)


