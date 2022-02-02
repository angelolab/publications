# Liu and Bosse et al., Reproducible, high-dimensional imaging in archival human tissue by Multiplexed Ion Beam Imaging by Time-of-Flight (MIBI-TOF)
# Create MPI/PPP scatterplots in Figure 3 and Supplementary Figure 4
# Plot CV within each core in Supplementary Figure 5

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(stringr)

# Read Ionpath data (normalized in mpi_ppp_ihc_regression.ipynb)
cimac_dat = fread("ionpath_norm_data.csv")

# Find mean of each ROI
fov_mean = cimac_dat[,lapply(.SD,mean),by=.(FOVName,Target),.SDcols=c("CalibratedPercentPixelPositive","CalibratedMeanFilteredIntensity")]
setnames(fov_mean,"CalibratedPercentPixelPositive","FOVMean_CalibratedPercentPixelPositive")
setnames(fov_mean,"CalibratedMeanFilteredIntensity","FOVMean_CalibratedMeanFilteredIntensity")
cimac_dat = fov_mean[cimac_dat,on=.(FOVName,Target)]
# Find sd of each FOV
fov_sd = cimac_dat[,lapply(.SD,sd),by=.(FOVName,Target),.SDcols=c("CalibratedPercentPixelPositive","CalibratedMeanFilteredIntensity")]
setnames(fov_sd,"CalibratedPercentPixelPositive","FOVSD_CalibratedPercentPixelPositive")
setnames(fov_sd,"CalibratedMeanFilteredIntensity","FOVSD_CalibratedMeanFilteredIntensity")
cimac_dat = fov_sd[cimac_dat,on=.(FOVName,Target)]
# Find cv of each FOV
cimac_dat[,FOVCV_CalibratedPercentPixelPositive := FOVSD_CalibratedPercentPixelPositive / FOVMean_CalibratedPercentPixelPositive]
cimac_dat[,FOVCV_CalibratedMeanFilteredIntensity := FOVSD_CalibratedMeanFilteredIntensity / FOVMean_CalibratedMeanFilteredIntensity]

# Remove dsDNA from plots
cimac_dat = cimac_dat[Target!="dsDNA"]

# Points to include in analysis
good_fovs = c('R1C2', 'R1C3', 'R1C5', 'R1C7', 'R1C9', 'R1C10',
             'R2C10', 'R2C11', 'R2C12', 'R3C2', 'R3C4', 'R6C7', 'R6C10',
             'R6C11', 'R7C6', 'R7C7', 'R7C10', 'R8C1', 'R8C10', 'R8C11')
cimac_dat = cimac_dat[FOVName %in% good_fovs]


## Create scatter plots for average MPI/PPI vs. MPI/PP
num_targets = length(unique(cimac_dat$Target))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# MPI
ggplot(cimac_dat, aes(x=FOVMean_CalibratedMeanFilteredIntensity, y=CalibratedMeanFilteredIntensity, color=Target)) +
  geom_abline(slope=1, intercept=0, color='gray') +
  geom_point() +
  scale_color_manual(values=getPalette(num_targets)) +
  scale_x_log10(limits=c(0.6,30)) +
  scale_y_log10(limits=c(0.6,30)) +
  annotation_logticks() +
  xlab("ROI Mean Positive Intensity, Average") +
  ylab("ROI Mean Positive Intensity") +
  theme_bw()
ggsave("mpi_plot.pdf",width=8,height=5)

# PPP
max_x = max(cimac_dat$FOVMean_CalibratedPercentPixelPositive)
max_y = max(cimac_dat$CalibratedPercentPixelPositive)
ggplot(cimac_dat, aes(x=FOVMean_CalibratedPercentPixelPositive/max_x*100, y=CalibratedPercentPixelPositive/max_y*100, color=Target)) +
  geom_abline(slope=1, intercept=0, color='gray') +
  geom_point() +
  scale_color_manual(values=getPalette(num_targets)) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  xlab("ROI Percent Positive Pixels, Average") +
  ylab("ROI Percent Positive Pixels") +
  theme_bw()
ggsave("ppp_plot.pdf",width=8,height=5)



## Boxplot of CVs
per_fov_dat = unique(cimac_dat[,c("FOVName","Target","FOVCV_CalibratedMeanFilteredIntensity")])
# Rename HLA1
setDT(per_fov_dat)[,Target:=str_replace(Target,"HLA class 1 A, B, and C, Na-K-ATPase alpha1","HLA1 + ATPase")]
ggplot(per_fov_dat, aes(x=Target, y=FOVCV_CalibratedMeanFilteredIntensity)) +
  geom_boxplot() +
  ggtitle("Mean Filtered Intensity") +
  ylab("ROI CV") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave("cv_boxplot.pdf", width=9, height=5)

