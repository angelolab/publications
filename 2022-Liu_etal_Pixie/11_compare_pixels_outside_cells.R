# Compare percent of pixels outside of cells in different tissues/datasets
# Author: Candace Liu
# Date: 8/15/22

library(data.table)
library(ggplot2)

dirs = c("ln","dcis","tnbc")
all_dat = rbindlist(lapply(dirs, function(x) {return(fread(file.path(x,"pixels_outside_cells.csv")))}))

# Get mean of each dataset
mean_dat = all_dat[,lapply(.SD,mean), by=group, .SDcols="perc"]

# Plot
ggplot(mean_dat, aes(x=reorder(group,perc), y=perc)) +
  geom_bar(stat="identity") +
  labs(x="", y="Average fraction of informative pixels outside of cells") +
  theme_light()

