# Volcano plot of DE glycans
# Adapted from: https://biostatsquid.com/volcano-plots-r-tutorial/
# Author: Ke Leow
# Date: 03/17/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggrepel)

#--------------------------------
# Load data 
#--------------------------------
df <- read_csv("data/DBDP_transition/MALDI/library_matched/DE_pval_glycans_DBDP.csv") %>% 
  #get -log10(pval)
  mutate(neglog10_pval = -log10(DB_DP_P.Value))

# set cutoff and groups to label as DE
logfc=0.25
adjpval=0.25
group1="DB"
group2="DP"
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NC"
#label up
df$diffexpressed[df$DB_DP_logFC > logfc & df$DB_DP_adj.P.Val < adjpval] <- group1
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$DB_DP_logFC < -logfc & df$DB_DP_adj.P.Val < adjpval] <- group2
#summarize number of DE features in each group
df %>% group_by(diffexpressed) %>% 
  summarise(de_features = n())

#add glycan annotations
glycan_data <- read_csv("data/glycan_peaklist_paperAnnotations_031825.csv") 
df<- inner_join(df,glycan_data %>% select(composition, Branches, type))
# df$type = factor(df$type, levels = c("Polylacnac","Tetraantennary","High Mannose","Other"))

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
topDE= df %>% filter(neglog10_pval > 1.3) %>% pull(Glycan)
  # pull(head(df[order(df$DB_DP_P.Value, decreasing = FALSE), "Glycan"], 30))
df$delabel <- ifelse(df$Glycan %in% topDE, df$composition, NA)
df$color_Branches <- ifelse(df$Glycan %in% topDE, df$Branches, NA) #color only labeled glycans
df$color_type <- ifelse(df$Glycan %in% topDE, df$type, NA) #color only labeled glycans

# Define custom colors
# custom_colors <- c("Partial.Gal"= "blue", "Polylacnac"= "lightblue","Tetraantennary"= "yellow","High Mannose"= "#2ca02c","Other"="grey")

#--------------------------------
# Plot Volcano
#--------------------------------
pdf("R_plots/DBDP_transition/volcano_DEglycan_colortypes_limma_alldonors.pdf", width = 8, height = 5)
ggplot(df, aes(x = DB_DP_logFC, y = neglog10_pval, col = color_type, label = delabel)) +
  geom_point(size = 3)+
  geom_text_repel(max.overlaps = Inf)+ # To show all labels 
  theme_classic(base_size = 12)+
  labs(x = "Log2(Fold Change",
       y = "-Log10(p-value)",
       legend = "Glycan type"
  )
  # xlab("log2fc") + 
  # ylab("-log10(pval)")
# theme(legend.position = "none")+
  # scale_color_manual(values = c("#5C0B18", "#2D4A72", "grey")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  # geom_vline(xintercept = c(-logfc, logfc), col = "gray", linetype = 'dashed') +
  # geom_hline(yintercept = 1.5, col = "gray", linetype = 'dashed') +
  # geom_text_repel(max.overlaps = Inf) # To show all labels 
dev.off()

