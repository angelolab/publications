# Calculate cluster consistency score, pixel level
# 5 replicates must have been run using pixelClustering.R before running this script
# Author: Candace Liu
# Date: 8/15/22

library(data.table)
library(ggplot2)

name = "sigma2_passes10" #name of output specified in pixelClustering.R
reps = 1:5 #number of replicates
filenames = paste0("pixelClustering_",name,"_rep",reps,"_clusters.csv") #filename for all replicates
comparison_names = paste0("rep",reps)

clust_coln = "hCluster_cap" #name of column that contains cluster ids
threshold = 95 #threshold for calculating cluster consistency score
out_name = paste0(name,"_threshold",threshold)

compare <- function(name, filenames, comparison_names, clust_coln, percent) {

  for (i in 1:length(filenames)) {
    one_dat = fread(filenames[i])
    comp_name = comparison_names[i]
    setnames(one_dat,clust_coln,comp_name)
    cnames = c("sample","x","y",comp_name)
    if (i == 1) {
      all_clusts = one_dat[,..cnames]
    } else {
      all_clusts = all_clusts[one_dat[,..cnames], on=.(sample,x,y)]
    }
  }
  
  # Pairwise comparisons
  one_comparison <- function(coln1, coln2, group, ground_truth) {
    # Get total number of clusters
    total_num = length(unique(all_clusts[,get(coln1)]))

    all_totals = rep(0,total_num)
    all_counts = rep(0,total_num)

    for (i in 1:total_num) {
      dat = all_clusts[get(coln1) == i]
      # Get total number of pixels in cluster i
      all_totals[i] = nrow(dat)

      # Count number of pixels in other rep
      counts = dat[,.N,by=eval(coln2)]
      counts = counts[order(-N)]

      # Get cumulative sum
      counts[, cumsum := cumsum(N)]
      # Find number closest to some% of all pixels
      diffs = counts$cumsum - (round(percent*nrow(dat)))
      all_counts[i] = which(diffs == min(diffs[diffs>=0]))
   }

    dt = data.table(cluster=1:total_num, total=all_totals, count=all_counts, group=group, ground_truth=coln1)
    return(dt)
  }

  # Make list of all comparisons
  rep_tab = data.table(expand.grid(comparison_names, comparison_names, stringsAsFactors = FALSE))
  rep_tab[,group := paste0(Var1,"_",Var2)]
  rep_tab = rep_tab[Var1 != Var2]

  # Run all comparisons
  all_comparisons = rep_tab[,one_comparison(Var1,Var2,group), by=seq_len(nrow(rep_tab))]

  # Mean of all comparisons for one ground truth  
  gt_comparisons = all_comparisons[,lapply(.SD,mean), by=.(cluster,ground_truth), .SDcols="count"]

  # Get pixel-level purity scores
  for (comp in comparison_names) {
    total_num = length(unique(all_clusts[,get(comp)]))
    for (i in 1:total_num) {
      all_clusts[get(comp) == i, paste0(comp,"_purity") := gt_comparisons[ground_truth==comp & cluster==i, "count"]]
    }
  }

  fwrite(all_clusts, paste0("purity_",out_name,".csv"))
}

compare(name, filenames, comparison_names, clust_coln, threshold/100)

