# Compare cluster consistency scores across parameter choices, get p-values
# Author: Candace Liu
# Date: 8/15/22

library(data.table)
library(ggpubr)

names = c("passes10","nofreqnorm","no999") #name of files to compare
filenames = paste0("meanpurity_", c("sigma2_passes10","sigma2_passes10_nofreqnorm","sigma2_passes10_no999"), "_threshold95.csv")

# Get all data
one_file <- function(i) {
  name = names[i]
  filename = filenames[i]
  one_dat = fread(filename)
  one_dat[,group:=name]
  return(one_dat[,c("mean_purity","group")])
}
all_dat = rbindlist(lapply(1:length(names),one_file))

# Do comparison
compare_means(mean_purity ~ group, data=all_dat, method="wilcox.test", paired=FALSE)
all_dat[,group:=factor(group,levels=names)]

# Get stats
med = all_dat[,.(med = median(mean_purity)), by=group]
q25 = all_dat[,.(q25 = quantile(mean_purity,0.25)), by=group]
q75 = all_dat[,.(q75 = quantile(mean_purity,0.75)), by=group]
min = all_dat[,.(min = min(mean_purity)), by=group]
max = all_dat[,.(max = max(mean_purity)), by=group]

# Plot in box plot
plot_dat = med[q25,on=.(group)]
plot_dat = plot_dat[q75,on=.(group)]
plot_dat = plot_dat[min,on=.(group)]
plot_dat = plot_dat[max,on=.(group)]
plot_dat[,iqr:=q75-q25]

ggplot(plot_dat) +
  geom_boxplot(aes(x=reorder(group,med), ymin=q25, lower=q25, middle=med, upper=q75, ymax=q75), stat="identity",fill='lightgray') +
  ylim(c(1,5)) +
  theme_light() +
  theme(axis.title.x=element_blank())

