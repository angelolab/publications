# Get correlation of pixel cluster frequency between replicate serial sections
# Author: Candace Liu
# Date: 8/15/22

library(data.table)
library(ggplot2)

# Point data (CIMAC data from Liu and Bosee et al., 2022), multiple serial sections per tissue core
point_dat = fread("cimac_points.csv")
num_points = length(unique(point_dat$point))

# Pixel clustering data (pixel clusters generated using pixelClustering.R)
pixel_dat = fread("pixelClustering_sigma2_clusters.csv") #path to output of pixelClustering.R
mapping = fread("manual_clusters.csv") #mapping of clusters to phenotype (manually done)
pixel_dat = pixel_dat[mapping, on=.(cluster)]
k_clust = length(unique(pixel_dat$phenotype))

# Count number of pixel clusters per sample
pixel_freq = pixel_dat[,as.list(table(factor(phenotype))), by=.(sample)]
pixel_freq = pixel_freq[order(sample)]
pixel_freq[,sample:=NULL]

# Correlation between serial sections of the same tissue core
corrs = cor(t(pixel_freq), method="spearman")
rownames(corrs) = paste0("point_",1:num_points)
colnames(corrs) = paste0("point_",1:num_points)

# Loop through all tissue cores
all_corrs = data.table(index=character(), corrs=numeric())
for (ind in unique(point_dat$index)) {
  samps = point_dat[index==ind]$point
  if (length(samps) > 1) {
    samps_corr = corrs[samps,samps]
    samps_dt = data.table(index=ind, corrs=samps_corr[upper.tri(samps_corr)])
    all_corrs = rbindlist(list(all_corrs,samps_dt))
  }
}

# Plot
means = all_corrs[,lapply(.SD,mean),by=index]
ggplot(means, aes(x=corrs)) +
  geom_histogram(binwidth=0.04, center=0.02) +
  scale_y_continuous(breaks=0:12) +
  coord_cartesian(xlim=c(0,1)) +
  xlab("Spearman correlation") +
  ylab("Number of cores") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

