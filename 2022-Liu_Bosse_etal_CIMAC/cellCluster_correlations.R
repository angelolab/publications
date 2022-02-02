# Liu and Bosse et al., Reproducible, high-dimensional imaging in archival human tissue by Multiplexed Ion Beam Imaging by Time-of-Flight (MIBI-TOF)
# Calculate Spearman correlation between serial sections of same core using cell phenotype frequencies, Figure 4B
# Also plot Spearman correlations as a function of slide distance, Supplementary Figure 8

library(data.table)
library(ggplot2)
library(stringr)

# Read in cell clustering data (generated in cellClustering.R)
cell_dat = fread("cell_table_size_normalized_clusters.csv")

# Replace names some mislabeled files
setDT(cell_dat)[,fov:=str_replace(fov,"201001_Slide23Stain2Run_shuffled_Point26_R1C3","201001_Slide23Stain2Run_shuffled_Point26_R0C0")]
setDT(cell_dat)[,fov:=str_replace(fov,"201001_Slide31Stain4Run_restart_Point15_R3C4","201001_Slide31Stain4Run_restart_Point15_R3C3")]
setDT(cell_dat)[,fov:=str_replace(fov,"201001_Slide31Stain4Run_restart_Point16_R3C5","201001_Slide31Stain4Run_restart_Point16_R3C4")]
setDT(cell_dat)[,fov:=str_replace(fov,"201007_Slide21Stain1Run_shuffled_4b_Point10_R1C4","201007_Slide21Stain1Run_shuffled_4b_Point10_R4C1")]

# Count number of cell clusters per sample
cell_freq = cell_dat[, .N, by=.(fov,phenotype)]
cell_freq = dcast(cell_freq, fov~phenotype, value.var="N")
cell_freq[is.na(cell_freq)] = 0

# Save order of samples and extract core index/slide number from sample name
cell_freq = cell_freq[order(fov)]
sample_info = data.table(fov_num=1:length(cell_freq$fov), fov=cell_freq$fov)
sample_info[,core := str_extract(fov,"(?!.*_)\\w+")]
sample_info[,slide := str_extract(fov,"(?<=Slide)\\d+")]
# Create dictionary mapping slide number to fov_num
point_to_slide = sample_info$slide
names(point_to_slide) = as.character(sample_info$fov_num)

# Correlation between each sample
cell_freq[,fov:=NULL]
corrs = cor(t(cell_freq),method="spearman")
num_points = nrow(corrs)
rownames(corrs) = 1:num_points
colnames(corrs) = 1:num_points

# Get correlation within each core and plot correlation as a function of slide distance
all_corrs = data.table(core=character(), corrs=numeric(), point1=integer(), point2=integer(), slide1=integer(), slide2=integer(), min_slide=integer(), max_slide=integer(), dist=integer())
for (ind in unique(sample_info$core)) {
  ind_dat = sample_info[core==ind]
  samps = ind_dat$fov_num
  if (length(samps) > 1) {
    samps_corr = corrs[samps,samps]
    samps_names = expand.grid(rownames(samps_corr),colnames(samps_corr),stringsAsFactors=FALSE)[as.vector(upper.tri(samps_corr)),]
    samps_dt = data.table(core=ind, corrs=samps_corr[upper.tri(samps_corr)], point1=as.integer(samps_names[,1]), point2=as.integer(samps_names[,2]))

    samps_dt[,slide1:=as.integer(point_to_slide[as.character(point1)])]
    samps_dt[,slide2:=as.integer(point_to_slide[as.character(point2)])]
    samps_dt[,min_slide:=min(slide1,slide2), by=1:nrow(samps_dt)]
    samps_dt[,max_slide:=max(slide1,slide2), by=1:nrow(samps_dt)]
    samps_dt[,dist:=max_slide-min_slide]

    # Only keep correlations that start at the first slide (which is slide number 21)
    keep_dat = samps_dt[slide1==21 | slide2==21]
    p = ggplot(keep_dat, aes(x=dist, y=corrs)) +
          geom_point() +
          geom_line() +
          ylim(c(0.8,1)) +
          scale_x_continuous(breaks=seq(2,10,2)) +
          xlab("Slide distance") +
          ylab("Spearman correlation") +
          ggtitle(ind) +
          theme_bw() +
          theme(panel.grid.minor.x = element_blank())
    print(p)

    # Save all correlations
    all_corrs = rbindlist(list(all_corrs,samps_dt))
  }
}

# Plot histogram of all means
means = all_corrs[,lapply(.SD,mean),by=core,.SDcols="corrs"]
ggplot(means, aes(x=corrs)) +
  geom_histogram(binwidth=0.05, center=0.025) +
  scale_y_continuous(breaks=0:15) +
  xlim(c(0,1)) +
  xlab("Spearman correlation") +
  ylab("Number of cores") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("cellClustering_mean_correlation.pdf", width=6, height=5)


