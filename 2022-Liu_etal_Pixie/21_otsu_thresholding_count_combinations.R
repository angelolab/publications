# Count combinations after Otsu thresholding
# Author: Candace Liu
# Date: 12/15/22

library(arrow)
library(data.table)

# Read in all points
points = 1:12
dat = rbindlist(lapply(points, function(x) {return(data.table(read_feather(file.path("pixel_mats_thresholded",paste0("Point",x,"_sigma2.feather")))))}))

# Get combinations
markers = c("CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31")
dat_markers = dat[,..markers]
dat_markers[,rownum:=1:.N,]
dat_markers$class = melt(dat_markers, "rownum")[,toString(variable[value==1]),rownum]$V1
dat_markers[,class:=gsub(", ","_",class)]
dat_markers[,rownum:=NULL]

# Get total number of positive markers per pixel
dat_markers$num_pos = rowSums(dat[,..markers])
dat = cbind(dat[,c("sample","x","y")], dat_markers[,c("class","num_pos")])

# Count combinations
counts_bysample = dat[,.N,by=.(sample,class,num_pos)]
counts_total = dat[,.N,by=.(class,num_pos)]
counts_numpos = dat[,.N,by=.(num_pos)]

fwrite(counts_bysample, "counts_bysample.csv")
fwrite(counts_total, "counts_total.csv")
fwrite(counts_numpos, "counts_numpos.csv")

