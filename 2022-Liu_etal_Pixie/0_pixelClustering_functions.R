# Functions for pixel clustering and cell clustering
# Author: Candace Lu
# Date: 8/15/22

library(FlowSOM)
library(data.table)
library(dplyr)
library(viridis)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(ConsensusClusterPlus)
library(arrow)

# Helper function for doing 99.9% normaliztion
quantile_norm <- function(x) {
  q = quantile(x[x!=0], 0.999)

  if (is.na(q)) { #NA if x is all 0s
    return(x)
  }

  if (q!=0) {
    return(x/q)
  } else {
    return(x/max(x))
  }
}


# Helper function for getting 99.9 percentiles
get_quantile_norm <- function(x) {
  q = quantile(x[x!=0], 0.999)

  if (is.na(q)) { #NA if x is all 0s
    return(0)
  }

  if (q!=0) {
    return(q)
  } else {
    return(max(x))
  }
}


# Get average percentile of all points
# Inputs:
#   points - list of all points
#   markers - list of markers to cluster on
#   pixel_mat_dir - directory where output of imageToMatrix.py is saved
#   sigma - sigma used for Gaussian blur in imageToMatrix.py
#   quantile - quantile for normalization
# Output:
#   data.table with all quantiles
get_point_avg_quantile <- function(points, markers, pixel_mat_dir, sigma, quantile = 0.999) {

  one_point <- function(point_num) {
    # Read data
    one_tab = fread(file.path(pixel_mat_dir,paste0("Point", point_num, "_sigma", sigma,".csv")))
    one_tab = one_tab[rowSums(one_tab[,..markers])>0, ] # remove empty pixels
  
    # Turn into frequency
    tab_markers = one_tab[,..markers]
    dat = tab_markers[,lapply(.SD, function(x) x/rowSums(tab_markers))]
  
    # Get 99.9 percentile
    all_quants = dat[,lapply(.SD, function(x) get_quantile_norm(x))]
  
    return(all_quants)
  }

  all_point_quant = rbindlist(lapply(points, function(x) one_point(x)))
  mean_quant = all_point_quant[,lapply(.SD,mean), .SDcols=markers]
  return(mean_quant)
}


# Use FlowSOM to cluster pixels, then consensus hierarchical cluster with and without z-score capping
# Inputs:
#   points - list of all points
#   markers - list of markers to cluster on
#   norm_vals - 99.9% normalization values
#   pixel_mat_dir - directory where output of imageToMatrix.py is saved
#   sigma - sigma used for Gaussian blur in imageToMatrix.py
#   k - number of pixel meta-clusters
#   cap - z-score cap to use when hierarchical clustering
#   xdim - width of the SOM grid
#   ydim - height of the SOM grid
#   npasses - number of passes through data for FlowSOM
#   alpha - learning rate for FlowSOM
#   mst - number of times to build mst
#   seed - for consensus clustering reproducibility
#   file_suffix - to append to file name
#   save_dat - whether to write data to file
#   save_som - whether to save SOM R object
# Output:
#   pixels x markers data.table, with columns for FlowSom cluster, hCluster, and hCluster_cap
pixelClustering_flowsom_consensus <- function(points, markers, norm_vals, pixel_mat_dir = "pixel_mats", sigma = 2, k = 20, cap = 3, xdim = 10, ydim = 10, npasses = 10, alpha = c(0.05,0.01), mst = 1, seed = 59, file_suffix = "", save_dat = FALSE, save_som = FALSE) {

  if (file_suffix != "") {
    file_suffix = paste0("_",file_suffix)
  }

  # Set order of normalization coefficients to match markers
  norm_vals = norm_vals[,..markers]

  # Read in pixel x marker data (data has been Guassian blurred in Python)
  tab = rbindlist(lapply(points, function(x) {print(x)
                                              one_tab = fread(file.path(pixel_mat_dir,paste0("Point", x, "_sigma", sigma,".csv")))
                                              one_tab = one_tab[rowSums(one_tab[,..markers])>0, ] # remove empty pixels
                                              return(select(one_tab, one_of(c('sample','x','y','label',markers))))}))
  # Pixel normalization
  tab_markers = tab[,..markers]
  dat = tab_markers[,lapply(.SD, function(x) x/rowSums(tab_markers))]

  # 99.9% marker normalization
  dat = dat[,Map(`/`,.SD,norm_vals)]

  # Cluster data using FlowSOM
  set.seed(seed)
  fSOM = SOM(data = as.matrix(dat), xdim = xdim, ydim = ydim, rlen = npasses, mst = mst, alpha = alpha)
  if (save_som == TRUE) {
    saveRDS(fSOM, file=paste0("pixelClustering_sigma",sigma,file_suffix,"_som.rds"))
  }
  dat$cluster = fSOM$mapping[,1]

  ## Get consensus hierarchical clusters with no z-score capping
  # Get mean of each cluster
  mean_dat = dat[, lapply(.SD, mean), by = cluster]
  mat_dat = data.frame(mean_dat[, ..markers])
  # Z-score the columns
  mat_dat = scale(mat_dat)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k, seed=seed)
  clust_to_hclust = consensus[[k]]$consensusClass
  names(clust_to_hclust) = mean_dat$cluster
  dat$hCluster = clust_to_hclust[as.character(dat$cluster)]

  ## Get hierarchical clusters with capping
  mat_dat = pmin(mat_dat, cap)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k, seed=seed)
  clust_to_hclust = consensus[[k]]$consensusClass
  names(clust_to_hclust) = mean_dat$cluster
  dat$hCluster_cap = clust_to_hclust[as.character(dat$cluster)]

  # Add back other columns
  to_add = colnames(tab)[colnames(tab) %in% c("sample","x","y","label")]
  out = cbind(tab[,..to_add], dat)
  out_coln = colnames(out)[!colnames(out) %in% markers]

  if (save_dat == TRUE) {
    fwrite(out, paste0("pixelClustering_sigma",sigma,file_suffix,".csv"))
  }
  fwrite(out[,..out_coln], paste0("pixelClustering_sigma",sigma,file_suffix,"_clusters.csv"))

  # Save table that maps cluster to hCluster
  clust_tab = data.table(cluster = mean_dat$cluster, hCluster_cap = clust_to_hclust[as.character(mean_dat$cluster)])
  fwrite(clust_tab, paste0("pixelClustering_sigma",sigma,file_suffix,"_clust_to_hclust.csv"))

  return(out)
}


# Makes plots for pixel clustering
# Input:
#   pixel_dat - pixel x marker data.table, with columns for cluster and hCluster
#   hclust_coln - column name in pixel_dat for hierarchal clusters (hCluster or hCluster_cap, probably want hCluster_cap)
#   colors - list of colors to use for annotations in heatmaps (need to be at least as long as number of hClusters)
#   markers - markers used for clustering
#   cap - z-score cap used when making hierarchical clusters
#   save_plots - whether to save heatmaps as pdf
#   file_suffix - file suffix for filenames
# Output:
#   list of 5 plots -
#     1. Histogram of # pixels per pixel cluster
#     2. Histogram of # pixels per pixel metacluster
#     3. Histogram of # clusters per pixel metacluster
#     4. Heatmap of pixel clusters x markers
#     5. Heatmap of pixel metaclusters x markers
pixelClustering_plots <- function(pixel_dat, hclust_coln, colors, markers, cap, save_plots=FALSE, file_suffix="") {
  # Get number of hierarchical clusters
  all_hclusts = sort(unique(pixel_dat[,get(hclust_coln)]))
  k_hclust = length(all_hclusts)
  # Annotation colors for each k_hclust to use in heatmaps
  mat_colors = colors[1:k_hclust]
  names(mat_colors) = paste0("pixel_h", all_hclusts)
  mat_colors = list(hclust = mat_colors)

  if (file_suffix != "") {
    file_suffix = paste0("_",file_suffix)
  }

  # 1. Histogram of counts in each pixel cluster
  counts = table(pixel_dat$cluster)
  df = data.frame(clust=as.numeric(names(counts)), count=as.vector(counts))
  p1 = ggplot(df, aes(x=clust, y=count)) +
         geom_bar(stat="identity") +
         ggtitle("Number of pixels in each cluster") +
         xlab("FlowSOM cluster") +
         ylab("Count") +
         theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())

  # 2. Histogram of counts in each pixel hierarchical cluster
  hcounts = table(pixel_dat[,get(hclust_coln)])
  df = data.frame(clust=paste0("pixel_h",as.numeric(names(hcounts))), count=as.vector(hcounts))
  df$clust = factor(df$clust, levels=df$clust)
  p2 = ggplot(df, aes(x=clust, y=count)) +
         geom_bar(stat="identity") +
         ggtitle("Number of pixels in each hierarchical cluster") +
         ylab("Count") +
         theme(axis.text.x = element_text(angle=90),
               axis.title.x = element_blank())

  # 3. Histogram of number of clusters in each hierarchical cluster
  clusts = unique(data.table(cluster = pixel_dat$cluster, hCluster = pixel_dat[,get(hclust_coln)]))
  count_clust = table(clusts$hCluster)
  df = data.frame(clust=paste0("pixel_h",names(count_clust)), count=as.vector(count_clust))
  df$clust = factor(df$clust, levels=df$clust)
  p3 = ggplot(df, aes(x=clust, y=count)) +
         geom_bar(stat="identity") +
         ggtitle("Number of clusters in each hierarchical cluster") +
         ylab("Count") +
         theme(axis.text.x = element_text(angle=90),
               axis.title.x = element_blank())

  # 4. Heatmap of pixel clusters x markers, average across pixel clusters
  mean_dat = pixel_dat[, lapply(.SD, mean), by = cluster]
  coln = c("cluster",markers,hclust_coln)
  fwrite(mean_dat[,..coln],paste0("pixelClustering_cluster_mean",file_suffix,".csv"))
  mean_dat = mean_dat[order(get(hclust_coln))]
  mat_dat = data.frame(mean_dat[, ..markers])
  rownames(mat_dat) = paste0("clust_", mean_dat$cluster)
  # Z-score and cap
  mat_dat = scale(mat_dat)
  mat_dat = pmin(mat_dat, cap)
  # Annotations
  mat_col = data.frame(hclust = paste0("pixel_h", mean_dat[,get(hclust_coln)]))
  rownames(mat_col) = paste0("clust_", mean_dat$cluster)
  # Determine breaks
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  # Make heatmap
  p4 = pheatmap(t(mat_dat),
                color = viridis(length(breaks)-1),
                breaks = breaks,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_colnames = FALSE,
                annotation_col = mat_col,
                annotation_colors = mat_colors,
                main = "Average across pixel clusters")

  if (save_plots==TRUE) {
    pheatmap(mat_dat,
             color = viridis(length(breaks)-1),
             breaks = breaks,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             annotation_row = mat_col,
             annotation_colors = mat_colors,
             main = "Average across pixel clusters",
             filename = paste0("pixelClustering_cluster",file_suffix,"_heatmap.pdf"),
             height = 14,
             width = 8)
  }

  # 5. Heatmap of pixel hierarchical cluster x markers, average across hierarchical clusters
  mean_dat = pixel_dat[, lapply(.SD, mean), by = eval(hclust_coln), .SDcols = markers]
  fwrite(mean_dat,paste0("pixelClustering_hCluster_mean",file_suffix,".csv"))
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat[,get(hclust_coln)])
  # Z-score the columns
  mat_dat = scale(mat_dat)
  # Make annotations
  mat_col = data.frame(hclust = paste0("pixel_h",mean_dat[,get(hclust_coln)]))
  rownames(mat_col) = paste0("clust_",mean_dat[,get(hclust_coln)])
  # Make labels
  hcounts_order = hcounts[as.character(mean_dat[,get(hclust_coln)])]
  hcounts_labs = paste0(paste0("pixel_h",names(hcounts_order),": "), hcounts_order)
  # Determine breaks
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  # Make heatmap
  p5 = pheatmap(t(mat_dat),
                color = viridis(length(breaks)-1),
                breaks = breaks,
                cluster_rows = FALSE,
                labels_col = hcounts_labs,
                annotation_col = mat_col,
                annotation_colors = mat_colors,
                main = "Average across pixel hierarchical clusters")

  if (save_plots==TRUE) {
    pheatmap(mat_dat,
             color = viridis(length(breaks)-1),
             breaks = breaks,
             cluster_cols = FALSE,
             labels_row = hcounts_labs,
             annotation_row = mat_col,
             annotation_colors = mat_colors,
             main = "Average across pixel hierarchical clusters",
             filename = paste0("pixelClustering_hCluster",file_suffix,"_heatmap.pdf"),
             height = 8,
             width = 10)
  }

  return(list(p1,p2,p3,p4[[4]],p5[[4]]))
}


# Count the number of pixel clusters per cell (need segmentation labels)
# Inputs:
#   pixel_dat - MUST contain "label" column (need segmentation labels in imageToMatrix.py)
#   coln - column name in pixel_dat to count ("cluster", "hCluster_cap")
#   cell_dat - data table with columns from segmentation data that you want to map to pixel clustering data, MUST contain "sample", "label", and "cell_size" columns
#   save_dat - whether to write counts to file
# Outputs:
#   cell x pixel cluster data.table with counts of each pixel cluster in each cell
pixelClusters_count <- function(pixel_dat, coln, cell_dat, file_suffix = "", save_dat = TRUE) {
  # Count number of each pixel cluster
  k_clust = length(unique(pixel_dat[,get(coln)]))
  count_tab = pixel_dat[, as.list(table(factor(get(coln), levels=1:k_clust))), by=.(sample,label)]
  # Change column names
  colnames(count_tab) = c("sample","label",paste0("pixel_", 1:k_clust))
  # Combine with cell clustering data
  clust_freq = count_tab[cell_dat, on=.(sample,label), nomatch=0]

  if (file_suffix != "") {
    file_suffix = paste0("_",file_suffix)
  }

  if (save_dat == TRUE) {
    fwrite(clust_freq, paste0(coln,"_freq",file_suffix,".csv"))
  }

  return(clust_freq)
}


# FlowSOM clustering of cells using pixel cluster frequency, consensus hierarchical clustering
# Input:
#   dat - table of count of each pixel cluster per cell (ouput of pixelClusters_count)
#   k_cell_hclust - number of cell metaclusters
#   cap - z-score cap to use for hierarchical clusters
#   quant_norm - quantile to normalize columns to (usually 99.9%, but including option here because 99.9% sometimes returns NA's)
#   xdim - width of the SOM grid
#   ydim - height of the SOM grid
#   npasses - number of passes through data for FlowSOM
#   alpha - learning rate for FlowSOM
#   mst - number of times to build mst
#   seed - seed for consensus clustering reproducibility
#   file_suffix - to append to file name
#   save_som - whether to save SOM R object
# Output:
#   point x pixel cluster data table with columns for cell cluster and cell metaclusters
cellClustering_flowsom_consensus <- function(dat, k_cell_hclust = 20, cap = 3, quant_norm = 0.999, xdim = 10, ydim = 10, npasses = 10, alpha = c(0.05,0.01), mst = 1, seed = 59, file_suffix = "", save_som = FALSE) {

  if (file_suffix != "") {
    file_suffix = paste0("_",file_suffix)
  }

  # Remove columns that are all 0
  dat = dat[,colSums(dat)!=0, with=FALSE]

  # Get column names of just pixel clusters
  coln = colnames(dat)[grepl("pixel_", colnames(dat), fixed=TRUE)]

  # 99.9% normalization
  dat_trans = dat[,..coln] / dat$cell_size
  dat_trans = dat_trans[, lapply(.SD, function(x) quantile_norm(x))]

  # Cluster data using FlowSOM
  set.seed(seed)
  fSOM = SOM(data = as.matrix(dat_trans), xdim = xdim, ydim = ydim, rlen = npasses, mst = mst, alpha = alpha)
  if (save_som == TRUE) {
    saveRDS(fSOM, file=paste0("cellClustering",file_suffix,"_som.rds"))
  }
  dat$pixelfreq_cluster = fSOM$mapping[,1]
  dat_trans$pixelfreq_cluster = fSOM$mapping[,1]

  ## Get hierarchical clusters
  # Average for each cluster
  mean_dat = dat_trans[, lapply(.SD, mean), by = pixelfreq_cluster, .SDcols=coln]
  # Convert expression values back to data frame
  mat_dat = data.frame(mean_dat[,..coln])
  # Z-score the columns
  mat_dat = scale(mat_dat)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k_cell_hclust, seed=seed)
  clust_to_hclust = consensus[[k_cell_hclust]]$consensusClass
  names(clust_to_hclust) = mean_dat$pixelfreq_cluster
  dat$pixelfreq_hclust = clust_to_hclust[as.character(dat$pixelfreq_cluster)]

  ## Get hierarchical clusters by capping z-score
  mat_dat = pmin(mat_dat, cap)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k_cell_hclust, seed=seed)
  clust_to_hclust = consensus[[k_cell_hclust]]$consensusClass
  names(clust_to_hclust) = mean_dat$pixelfreq_cluster
  dat$pixelfreq_hclust_cap = clust_to_hclust[as.character(dat$pixelfreq_cluster)]

  out_dat = dat[,c("sample","label","pixelfreq_cluster","pixelfreq_hclust","pixelfreq_hclust_cap")]
  fwrite(out_dat, paste0("cellClustering_kcell",k_cell_hclust,file_suffix,".csv"))

  # Save table that maps cluster to hCluster
  clust_tab = data.table(pixelfreq_cluster = mean_dat$pixelfreq_cluster, pixelfreq_hclust_cap = clust_to_hclust[as.character(mean_dat$pixelfreq_cluster)])
  fwrite(clust_tab, paste0("cellClustering_kcell",k_cell_hclust,file_suffix,"_clust_to_hclust.csv"))

  return(out_dat)
}


# Make heatmaps for cell clusters
# Inputs:
#   clust_freq -  table of count of each pixel cluster per cell (output of pixelClusters_count, MUST have columns "sample","label","cell_size") 
#   cell_clusts - cell cluster number for each cell (output of cellClustering_flowsom)
#   cell_seg_dat - cell data from segmentation, MUST have columns named "sample" for point/roi/fov # and "label" for label from segmentation
#   pixel_mean_dat - marker x pixel cluster table (needed for weighted averages)
#   pixel_dat_coln - MUST match column name you used to count frequencies in pixelClusters_count ("cluster","hCluster_cap"), needs to be in pixel_mean_dat
#   markers - markers used for clustering
#   cluster_rows - T/F, whether to cluster rows in the heatmap
#   show_rownames - T/F, whether to show rownames in the heatmap
#   save_plots - T/F, whether to save heatmaps (weighted heatmaps)
#   file_suffix - suffix for filenames
# Output:
#   list of heatmaps -
#     1. Pixel metacluster (frequency) x cell cluster
#     2. Pixel metacluster (frequency) x cell metacluster
#     3. Marker x cell cluster (integrated)
#     4. Marker x cell metacluster (integrated)
#     5. Marker x cell cluster (weighted)
#     6. Marker x cell metacluster (weighted)
cellClustering_plots <- function(clust_freq, cell_clusts, cell_seg_dat, pixel_mean_dat, pixel_dat_coln, markers, cap = 3, mat_colors, cluster_rows = TRUE, show_rownames = FALSE, save_plots = FALSE, file_suffix = "") {
  # Remove columns that are all 0
  clust_freq = clust_freq[,colSums(clust_freq)!=0, with=FALSE]

  # Get column names of just pixel clusters
  coln = colnames(clust_freq)[grepl("pixel_", colnames(clust_freq), fixed=TRUE)]
  # Make colors for heatmaps
  all_cell_colors = mat_colors[1:length(unique(cell_clusts$pixelfreq_hclust))]
  names(all_cell_colors) = paste0("cell_h", 1:length(unique(cell_clusts$pixelfreq_hclust)))

  # Combine cluster frequency and cluster labels, normalize data
  freq_dat = clust_freq[,..coln] / clust_freq$cell_size
  freq_dat = cbind(freq_dat, clust_freq[,c("sample","label")])
  cell_clust_dat = freq_dat[cell_clusts, on=.(sample,label)]

  if (file_suffix != "") {
    file_suffix = paste0("_",file_suffix)
  }

  # 1. Pixel metacluster x cell cluster heatmaps, averages across cell clusters
  mean_dat = cell_clust_dat[, lapply(.SD, mean), by = pixelfreq_cluster, .SDcols = c(coln,"pixelfreq_hclust","pixelfreq_hclust_cap")]
  mat_dat = data.frame(mean_dat[,..coln])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust), hclust_cap = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_cluster)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust)))], hclust_cap = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust_cap)))])
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p = pheatmap(t(mat_dat),
               color = viridis(length(breaks)-1),
               breaks = breaks,
               cluster_rows = cluster_rows,
               show_rownames = show_rownames,
               show_colnames = FALSE,
               annotation_col = mat_col,
               annotation_colors = mat_colors,
               annotation_legend = FALSE,
               main = "Average across cell clusters")
  # Cap
  mat_dat = pmin(mat_dat, cap)
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cap = pheatmap(t(mat_dat),
                   color = viridis(length(breaks)-1),
                   breaks = breaks,
                   cluster_rows = cluster_rows,
                   show_rownames = show_rownames,
                   show_colnames = FALSE,
                   annotation_col = mat_col,
                   annotation_colors = mat_colors,
                   annotation_legend = FALSE,
                   main = "Average across cell clusters, cap")
  plot1 = grid.arrange(p[[4]], p_cap[[4]], ncol=2)

  # 2. Pixel metacluster x cell metacluster heatmaps, average across cell metaclusters
  mean_dat = cell_clust_dat[, lapply(.SD, mean), by = pixelfreq_hclust, .SDcols = coln]
  mat_dat = data.frame(mean_dat[,..coln])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_hclust)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_hclust)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust)))])
  # Make labels
  hcounts = table(cell_clust_dat$pixelfreq_hclust)
  hcounts_order = hcounts[as.character(mean_dat$pixelfreq_hclust)]
  hcounts_labs = paste0(paste0("cell_h",names(hcounts_order),": "), hcounts_order)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p = pheatmap(t(mat_dat),
               color = viridis(length(breaks)-1),
               breaks = breaks,
               cluster_rows = cluster_rows,
               show_rownames = show_rownames,
               labels_col = hcounts_labs,
               annotation_col = mat_col,
               annotation_colors = mat_colors,
               main = "Average across cell hierarchical clusters")
  # Use averages across hierarchical clusters using cap
  mean_dat = cell_clust_dat[, lapply(.SD, mean), by = pixelfreq_hclust_cap, .SDcols = coln]
  mat_dat = data.frame(mean_dat[,..coln])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_hclust_cap)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_hclust_cap)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust_cap)))])
  # Make labels
  hcounts = table(cell_clust_dat$pixelfreq_hclust_cap)
  hcounts_order = hcounts[as.character(mean_dat$pixelfreq_hclust_cap)]
  hcounts_labs = paste0(paste0("cell_h",names(hcounts_order),": "), hcounts_order)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cap = pheatmap(t(mat_dat),
                   color = viridis(length(breaks)-1),
                   breaks = breaks,
                   cluster_rows = cluster_rows,
                   show_rownames = show_rownames,
                   labels_col = hcounts_labs,
                   annotation_col = mat_col,
                   annotation_colors = mat_colors,
                   main = "Average across cell hierarchical clusters, cap")
  plot2 = grid.arrange(p[[4]], p_cap[[4]], ncol=2)

  # 3. Marker x cell cluster, average across cell clusters (directly by integrating over cell)
  # Merge integrated segmentation data with cell clusters
  dat_clusts_only = cell_clust_dat[,c("pixelfreq_cluster", "pixelfreq_hclust", "pixelfreq_hclust_cap", "sample", "label")]
  dat_cell = dat_clusts_only[cell_seg_dat, on=c("sample","label"), nomatch=0]
  # Average for each cluster
  mean_coln = c(markers, "pixelfreq_hclust", "pixelfreq_hclust_cap")
  mean_dat = dat_cell[, lapply(.SD, mean), .SDcols = mean_coln, by = pixelfreq_cluster]
  mean_dat = mean_dat[order(pixelfreq_hclust)]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust), hclust_cap = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_cluster)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust)))], hclust_cap = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust_cap)))])
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell = pheatmap(t(mat_dat),
                    color = viridis(length(breaks)-1),
                    breaks = breaks,
                    cluster_rows = cluster_rows,
                    cluster_cols = FALSE,
                    show_colnames = FALSE,
                    annotation_col = mat_col,
                    annotation_colors = mat_colors,
                    annotation_legend = FALSE,
                    main = "Average across cell clusters (integrated)")
  # Make heatmap of marker x cell cluster, cap
  mean_dat = mean_dat[order(pixelfreq_hclust_cap)]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Z-score and cap
  mat_dat = scale(mat_dat)
  mat_dat = pmin(mat_dat, cap)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust), hclust_cap = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell_cap = pheatmap(t(mat_dat),
                        color = viridis(length(breaks)-1),
                        breaks = breaks,
                        cluster_rows = cluster_rows,
                        cluster_cols = FALSE,
                        show_colnames = FALSE,
                        annotation_col = mat_col,
                        annotation_colors = mat_colors,
                        annotation_legend = FALSE,
                        main = "Average across cell clusters (integrated), cap")
  plot3 = grid.arrange(p_cell[[4]], p_cell_cap[[4]], ncol=2)

  # 4. Marker x cell metacluster, average across cell metaclusters  (directly by integrating over cell)
  mean_dat = dat_cell[, lapply(.SD, mean), .SDcols = markers, by = pixelfreq_hclust]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_hclust)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_hclust)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust)))])
  # Make labels
  hcounts = table(dat_cell$pixelfreq_hclust)
  hcounts_order = hcounts[as.character(mean_dat$pixelfreq_hclust)]
  hcounts_labs = paste0(paste0("cell_h",names(hcounts_order),": "), hcounts_order)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell = pheatmap(t(mat_dat),
                    color = viridis(length(breaks)-1),
                    breaks = breaks,
                    cluster_rows = FALSE,
                    labels_col = hcounts_labs,
                    annotation_col = mat_col,
                    annotation_colors = mat_colors,
                    main = "Average across cell hierarchical clusters (integrated)")
  # Make heatmap of marker x cell hierarchical cluster, cap
  mean_dat = dat_cell[, lapply(.SD, mean), .SDcols = markers, by = pixelfreq_hclust_cap]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_hclust_cap)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_hclust_cap)
  # Make labels
  hcounts = table(dat_cell$pixelfreq_hclust_cap)
  hcounts_order = hcounts[as.character(mean_dat$pixelfreq_hclust_cap)]
  hcounts_labs = paste0(paste0("cell_h",names(hcounts_order),": "), hcounts_order)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell_cap = pheatmap(t(mat_dat),
                        color = viridis(length(breaks)-1),
                        breaks = breaks,
                        cluster_rows = FALSE,
                        labels_col = hcounts_labs,
                        annotation_col = mat_col,
                        annotation_colors = mat_colors,
                        main = "Average across cell hierarchical clusters (integrated), cap")
  plot4 = grid.arrange(p_cell[[4]], p_cell_cap[[4]], ncol=2)

  # 5. Marker x cell cluster, average across cell clusters (weighted, frequency of pixel cluster x average expression of pixel cluster)
  # Wrangle cell cluster info
  cell_info = clust_freq[,c("sample","label","cell_size")]
  # Wrangle pixel data
  mean_dat = pixel_mean_dat[order(get(pixel_dat_coln))]
  # Make sure pixel clusters are in same order in frequency table and in pixel data 
  coln_order = paste0("pixel_",mean_dat[,get(pixel_dat_coln)])
  coln_include = coln_order %in% coln
  coln_order = coln_order[coln_include]

  # Multiply # of each cluster by average expression of each pixel cluster
  weighted = data.table(as.matrix(clust_freq[,..coln_order]) %*% as.matrix(mean_dat[coln_include,..markers]))
  weighted_cell_info = cbind(weighted, cell_info)
  # Normalize by size
  dat_weighted = weighted_cell_info[,..markers] / weighted_cell_info$cell_size
  dat_weighted[,sample := weighted_cell_info$sample]
  dat_weighted[,label := weighted_cell_info$label]
  dat_weighted = dat_weighted[cell_clusts[,c("pixelfreq_cluster","pixelfreq_hclust","pixelfreq_hclust_cap","sample","label")], on=.(sample,label)]
  fwrite(dat_weighted,paste0("cellClustering_weightedMean",file_suffix,".csv"))

  # Average for each cluster
  mean_coln = c(markers, "pixelfreq_hclust", "pixelfreq_hclust_cap")
  mean_dat = dat_weighted[, lapply(.SD, mean), .SDcols = mean_coln, by = pixelfreq_cluster]
  fwrite(mean_dat,paste0("cellClustering_cluster_weightedMean",file_suffix,".csv"))
  mean_dat = mean_dat[order(pixelfreq_hclust)]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust), hclust_cap = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_cluster)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust)))], hclust_cap = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust_cap)))])
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell = pheatmap(t(mat_dat),
                    color = viridis(length(breaks)-1),
                    breaks = breaks,
                    cluster_rows = cluster_rows,
                    cluster_cols = FALSE,
                    show_colnames = FALSE,
                    annotation_col = mat_col,
                    annotation_colors = mat_colors,
                    annotation_legend = FALSE,
                    main = "Average across cell clusters")
  # Make heatmap of marker x cell cluster, cap
  mean_dat = mean_dat[order(pixelfreq_hclust_cap)]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Z-score and cap
  mat_dat = scale(mat_dat)
  mat_dat = pmin(mat_dat, cap)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust), hclust_cap = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_cluster)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell_cap = pheatmap(t(mat_dat),
                        color = viridis(length(breaks)-1),
                        breaks = breaks,
                        cluster_rows = cluster_rows,
                        cluster_cols = FALSE,
                        show_colnames = FALSE,
                        annotation_col = mat_col,
                        annotation_colors = mat_colors,
                        annotation_legend = FALSE,
                        main = "Average across cell clusters, cap")
  if (save_plots==TRUE) {
    # Annotations
    mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
    rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_cluster)
    mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust_cap)))])

    pheatmap(mat_dat,
             color = viridis(length(breaks)-1),
             breaks = breaks,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             annotation_row = mat_col,
             annotation_colors = mat_colors,
             main = "Average across cell clusters, cap",
             filename = paste0("cellClustering_cluster",file_suffix,"_heatmap.pdf"),
             height = 14,
             width = 8)
  }

  plot5 = grid.arrange(p_cell[[4]], p_cell_cap[[4]], ncol=2)

  # 6. Marker x cell metacluster, average across cell clusters (weighted, frequency of pixel cluster x average expression of pixel cluster)
  mean_dat = dat_weighted[, lapply(.SD, mean), .SDcols = markers, by = pixelfreq_hclust]
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_hclust)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_hclust)
  mat_colors = list(hclust = all_cell_colors[paste0("cell_h", sort(unique(mean_dat$pixelfreq_hclust)))])
  # Make labels
  hcounts = table(dat_weighted$pixelfreq_hclust)
  hcounts_order = hcounts[as.character(mean_dat$pixelfreq_hclust)]
  hcounts_labs = paste0(paste0("cell_h",names(hcounts_order),": "), hcounts_order)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell = pheatmap(t(mat_dat),
                    color = viridis(length(breaks)-1),
                    breaks = breaks,
                    cluster_rows = FALSE,
                    labels_col = hcounts_labs,
                    annotation_col = mat_col,
                    annotation_colors = mat_colors,
                    main = "Average across cell hierarchical clusters")
  # Make heatmap of marker x cell hierarchical cluster, cap
  mean_dat = dat_weighted[, lapply(.SD, mean), .SDcols = markers, by = pixelfreq_hclust_cap]
  fwrite(mean_dat,paste0("cellClustering_hCluster_weightedMean",file_suffix,".csv"))
  mat_dat = data.frame(mean_dat[,..markers])
  rownames(mat_dat) = paste0("clust_",mean_dat$pixelfreq_hclust_cap)
  # Z-score
  mat_dat = scale(mat_dat)
  # Annotations
  mat_col = data.frame(hclust = paste0("cell_h",mean_dat$pixelfreq_hclust_cap))
  rownames(mat_col) = paste0("clust_",mean_dat$pixelfreq_hclust_cap)
  # Make labels
  hcounts = table(dat_weighted$pixelfreq_hclust_cap)
  hcounts_order = hcounts[as.character(mean_dat$pixelfreq_hclust_cap)]
  hcounts_labs = paste0(paste0("cell_h",names(hcounts_order),": "), hcounts_order)
  # Make heatmap
  breaks = seq(min(mat_dat,na.rm=TRUE), max(mat_dat,na.rm=TRUE), length.out=100)
  p_cell_cap = pheatmap(t(mat_dat),
                        color = viridis(length(breaks)-1),
                        breaks = breaks,
                        cluster_rows = FALSE,
                        labels_col = hcounts_labs,
                        annotation_col = mat_col,
                        annotation_colors = mat_colors,
                        main = "Average across cell hierarchical clusters, cap")
  if (save_plots==TRUE) {
    pheatmap(mat_dat,
             color = viridis(length(breaks)-1),
             breaks = breaks,
             cluster_cols = FALSE,
             labels_row = hcounts_labs,
             annotation_row = mat_col,
             annotation_colors = mat_colors,
             main = "Average across cell hierarchical clusters, cap",
             filename = paste0("cellClustering_hCluster",file_suffix,"_heatmap.pdf"),
             height = 8,
             width = 10)
  }

  plot6 = grid.arrange(p_cell[[4]], p_cell_cap[[4]], ncol=2)

  return(list(plot1,plot2,plot3,plot4,plot5,plot6))
}


