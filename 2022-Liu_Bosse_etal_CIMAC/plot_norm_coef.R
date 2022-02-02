# Liu and Bosse et al., Reproducible, high-dimensional imaging in archival human tissue by Multiplexed Ion Beam Imaging by Time-of-Flight (MIBI-TOF)
# Create normalization coefficient plots in Supplementary Figure 2

library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)

# Read in normalization coefficients (table created in mpi_ppp_ihc_regression.ipynb)
norm = fread("norm_coefs.csv")
norm_melt = melt(norm)
setnames(norm_melt,"V1","Run")
setnames(norm_melt,"variable","FOV")
setnames(norm_melt,"value","norm_coef")
norm_melt[,name:=paste0(Run,"_",FOV)]

# Points to include in analysis
good_fovs = c('R1C2', 'R1C3', 'R1C5', 'R1C7', 'R1C9', 'R1C10',
              'R2C10', 'R2C11', 'R2C12', 'R3C2', 'R3C4', 'R6C7', 'R6C10',
              'R6C11', 'R7C6', 'R7C7', 'R7C10', 'R8C1', 'R8C10', 'R8C11')
norm_melt = norm_melt[FOV %in% good_fovs]

## Make histogram of normalization coefficients
ggplot(norm_melt, aes(x=norm_coef)) +
  geom_histogram() +
  xlab("Normalization coefficient") +
  ylab("Count") +
  theme_bw()
ggsave("norm_coef_hist.pdf", height=4, width=6)


## Compare before and after normalization
# Read Ionpath data (normalized in mpi_ppp_ihc_regression.ipynb)
cimac_dat = fread("ionpath_norm_data.csv")
cimac_dat[,name:=paste0(Run,"_",FOVName)]
cimac_dat = cimac_dat[FOVName %in% good_fovs]
all_dat = norm_melt[,c("name","norm_coef")][cimac_dat, on=.(name)]

# Rename HLA1
setDT(all_dat)[,Target:=str_replace(Target,"HLA class 1 A, B, and C, Na-K-ATPase alpha1","HLA1 + ATPase")]

ggplot(all_dat, aes(x=norm_coef,y=PercentPixelPositive)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = ..rr.label..), geom = "label") +
  facet_wrap(~Target, nrow=3, scales='free') + 
  xlab("Normalization coefficient") +
  ylab("Percent Positive Pixel") +
  theme_bw()
ggsave("norm_ppp.pdf", height=6, width=12)
  
ggplot(all_dat, aes(x=norm_coef,y=CalibratedPercentPixelPositive)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = ..rr.label..), geom = "label") +
  facet_wrap(~Target, nrow=3, scales='free') +
  xlab("Normalization coefficient") +
  ylab("Calibrated Percent Positive Pixel") +
  theme_bw()
ggsave("norm_ppp_cal.pdf", height=6, width=12)

