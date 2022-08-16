# Test FlowSOM vs Leiden vs Phenograph time
# Author: Candace Liu
# Date: 8/16/22

library(data.table)
library(Seurat)
library(FlowSOM)
library(ggplot2)
library(Rphenograph)

set.seed(329)
all_num_cells = c(1000,10000,20000,50000,100000) #different number of observations to test
markers = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")

# Get all data
points = 1:12
sigma = 2
pixel_mat_dir = "pixel_mats" #where pixel matrices from imageToMatrix.py are located
dat_all_marks = rbindlist(lapply(points, function(x) {print(x)
                                                      one_tab = fread(file.path(pixel_mat_dir,paste0("Point", x, "_sigma", sigma,".csv")))
                                                      one_tab = one_tab[rowSums(one_tab[,..markers])>0, ] # remove empty pixels
                                                      return(one_tab)}))
dat = dat_all_marks[,..markers]

# Phenograph
phenograph_sample <- function(num_cells) {
  samp = dat[sample(.N, num_cells)]
  start_time = Sys.time()
  pheno_out = Rphenograph(samp)
  end_time = Sys.time()

  elapsed = as.numeric(difftime(end_time,start_time), units="secs")
  print(num_cells)
  print(elapsed)

  rm(samp)
  rm(pheno_out)

  return(elapsed[[1]])
}
test = lapply(100, phenograph_sample)
output = lapply(all_num_cells, phenograph_sample)
phenograph_dt = data.table(algo="phenograph",num_cells=all_num_cells, time=unlist(output))


# Leiden
leiden_sample <- function(num_cells) {
  samp = dat[sample(.N, num_cells)]
  start_time = Sys.time()
  # Get distance
  cell_dist = dist(samp)
  # Seurat's FindNeighbors function
  knn = FindNeighbors(cell_dist, k.param=10)
  # Seurat's FindClusters function (algorithm=4 is Leiden)
  leiden = FindClusters(knn$snn, algorithm=4, resolution=2)
  end_time = Sys.time()

  elapsed = as.numeric(difftime(end_time,start_time), units="secs")
  print(num_cells)
  print(elapsed)

  rm(samp)
  rm(cell_dist)
  rm(knn)
  rm(leiden)

  return(elapsed[[1]])
}
test = lapply(100, leiden_sample)
output = lapply(all_num_cells, leiden_sample)
leiden_dt = data.table(algo="leiden",num_cells=all_num_cells, time=unlist(output))


# FlowSOM
flowsom_sample <- function(num_cells) {
  samp = dat[sample(.N, num_cells)]
  start_time = Sys.time()
  fsom = SOM(as.matrix(samp))
  end_time = Sys.time()

  elapsed = as.numeric(difftime(end_time,start_time), units="secs")
  print(num_cells)
  print(elapsed)

  return(elapsed[[1]])
}
output = lapply(all_num_cells, flowsom_sample)
som_dt = data.table(algo="flowsom",num_cells=all_num_cells, time=unlist(output))


# Comparison
all_dt = rbindlist(list(phenograph_dt,leiden_dt,som_dt))
fwrite(all_dt,"som_leiden_phenograph_time_test.csv")
ggplot(all_dt, aes(x=num_cells, y=time, group=algo, color=algo)) +
  geom_point() +
  geom_line() +
  scale_fill_manual(values=c("#5E81AC","#A3BE8C","#BF616A")) +
  scale_color_manual(values=c("#5E81AC","#A3BE8C","#BF616A")) +
  theme_light() +
  xlab("Number of observations") +
  ylab("Time (seconds)")


