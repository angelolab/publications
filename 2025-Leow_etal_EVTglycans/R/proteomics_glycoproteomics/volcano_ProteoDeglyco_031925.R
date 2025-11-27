# Volcano plot of DeGlycoProteomics data 
# Author: Ke Leow
# Date: 03/02/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggrepel)

#--------------------------------
# Load data 
#--------------------------------
df <- read_csv('data/glycoproteomics/24-02-27_Angelo_Placenta_ProteomicsData.csv') 
deglyco <- read_csv('data/glycoproteomics/24-02-27_Angelo_Placenta_DeGlyco_Data_FirstTake.csv')%>% 
  mutate(glycopeptide=paste(Gene, `Peptide Sequence`, sep = "_")) 

#--------------------------------
# Pull comparison pair
#--------------------------------
input = "DB v.s. DP"

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NC"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$`Student's T-test Difference DP_DB` > 1 & df$`Student's T-test q-value DP_DB` < 0.05] <- "DP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$`Student's T-test Difference DP_DB` < -1 & df$`Student's T-test q-value DP_DB` < 0.05] <- "DB"
# flip position of DB and DP on x-axis (FC)
df$`Student's T-test Difference DB_DP` = -df$`Student's T-test Difference DP_DB`


# Filter to keep only the most significant entry per gene - and up in DB
# de_genes <- df %>%
#   group_by(Gene) %>%
#   slice_max(`-Log Student's T-test p-value DP_DB`, n = 1) %>%   # Keeps the entry with the lowest p-value for each gene
#   filter(diffexpressed == "DB") %>% 
#   ungroup() %>% pull(Genes)
# df$delabel <- ifelse(df$Genes %in% de_genes, df$Gene, NA)

# Find common genes (intersection)
common_genes <- intersect(deglyco$Gene, df$Genes)
df$deglyco <- ifelse(df$Genes %in% common_genes, "deglyco", NA)
df$deglyco <- factor(df$deglyco)

df$alpha_value <- ifelse(df$Genes %in% common_genes, 1, 0.5)

#volcano - without labels
pdf("R_plots/glycoproteomics/volcano_proteomics_DBDP.pdf", , width = 3, height = 3.5)
ggplot(df, aes(x = `Student's T-test Difference DB_DP`, y = `-Log Student's T-test p-value DP_DB`, col = diffexpressed)) +
  geom_point()+
  theme_classic(base_size = 13)+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(legend.position = "none")+
  xlab("Log2(Fold Change)") + 
  ylab("-Log10(p-value)")+
  xlim(-7.5, 7.5) +  # Set equal limits around 0
  scale_color_manual(values = c("#C41E3A", "#2D4A72", "grey")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = 1.3, col = "gray", linetype = 'dashed') + 
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank()        # Remove ticks
  )
dev.off()

#volcano - color if detected in deglyco
pdf("R_plots/glycoproteomics/volcano_proteomicsMappedGlycoproteins_DBDP.pdf", , width = 3, height = 3.5)
ggplot(df, aes(x = `Student's T-test Difference DB_DP`, y = `-Log Student's T-test p-value DP_DB`, col = deglyco, alpha = alpha_value)) +
  geom_point()+
  theme_classic(base_size = 13)+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(legend.position = "none")+
  xlab("Log2(Fold Change)") + 
  ylab("-Log10(p-value)")+
  xlim(-7.5, 7.5) +  # Set equal limits around 0
  scale_color_manual(values = c("#8B5B45FF", "grey")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = 1.3, col = "gray", linetype = 'dashed') + 
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank()        # Remove ticks
  )
dev.off()

#volcano with gene labels
# ggplot(df, aes(x = `Student's T-test Difference DB_DP`, y = `-Log Student's T-test p-value DP_DB`, col = diffexpressed)) +
#   geom_point()+
#   theme_classic(base_size = 15)+
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
#         legend.position = "none")+
#   # theme(legend.position = "none")+
#   xlab("Log2(Fold Change)") + 
#   ylab("-Log10(Pvalue)")+
#   xlim(-7.5, 7.5) +  # Set equal limits around 0
#   scale_color_manual(values = c("#5C0B18", "#2D4A72", "grey")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
#   geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = 1.3, col = "gray", linetype = 'dashed') +
#   geom_label_repel(aes(label = delabel), na.rm = TRUE,
#                    box.padding = 0.4, point.padding = 0.3,
#                    segment.color = 'grey50', max.overlaps = 20,
#                    fill = "white", color = "black", size = 3)

#--------------------------------
# Plot number of differential DB/DP in proteo v.s. deglyco
#--------------------------------
proteo_count <- df %>% 
  select(diffexpressed) %>% 
  mutate(dataset = "Proteins")
df %>% group_by(diffexpressed) %>% summarise(n =n())

deglyco_count <- df %>% 
  filter(deglyco == "deglyco") %>%
  select(diffexpressed) %>% 
  mutate(dataset = "Glycoproteins")
deglyco_count %>% group_by(diffexpressed) %>% summarise(n =n())

# Combine the datasets
combined_data <- bind_rows(proteo_count, deglyco_count)

# Count occurrences of DB, DP, NC in each dataset
count_data <- combined_data %>%
  group_by(dataset, diffexpressed) %>%
  summarise(count = n(), .groups = "drop") %>% 
  filter(diffexpressed != "NC")

count_data$dataset <- factor(count_data$dataset, levels = c("Proteins", "Glycoproteins"))
count_data$diffexpressed <- factor(count_data$diffexpressed, levels = c("DP", "DB"))

pdf("R_plots/glycoproteomics/barplot_counts_proteomicsMappedGlycoproteins_DBDP.pdf", , width = 4, height = 4)
ggplot(count_data, aes(x = dataset, y = count, fill = diffexpressed)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Dataset", y = "Count", fill = "Differential") +
  scale_fill_manual(values = c("#2D4A72","#C41E3A"))+
  theme_classic(base_size = 15)+
  theme(legend.position="top")
dev.off()

# #--------------------------------
# # Plot ratio of DB/DP in proteo v.s. deglyco
# #--------------------------------
# # Count occurrences of DB and DP in each dataset
# count_data <- combined_data %>%
#   filter(diffexpressed %in% c("DB", "DP")) %>%
#   group_by(dataset, diffexpressed) %>%
#   summarise(count = n(), .groups = "drop") %>%
#   pivot_wider(names_from = diffexpressed, values_from = count, values_fill = 0) %>%
#   mutate(ratio = DB / DP)  # Calculate ratio
# 
# count_data$dataset <- factor(count_data$dataset, levels = c("Proteins", "Glycoproteins"))
# 
# # Plot the ratio
# ggplot(count_data, aes(x = dataset, y = ratio)) +
#   geom_bar(stat = "identity", width = 0.6) +
#   labs(x = "Dataset", y = "Ratio of DB/DP ") +
#   theme_classic(base_size = 15) +
#   theme(legend.position = "none")
# 
# #--------------------------------
# # Combine with deglyco-correlation results
# #--------------------------------
# cor_df <- read_csv('data/glycoproteomics/corr_ProteoDeglyco_byPeptides.csv')
# 
# # Get median correlation coeff per gene and plot
# cor_median <- cor_df %>% 
#   group_by(Gene) %>% 
#   summarise(median_coefficient = median(Correlation))
# median(cor_median$median_coefficient)
# 
# df_map_cor <- full_join(df, cor_median, by = join_by(Genes == Gene))
# 
# df_map_cor %>% filter(deglyco == "deglyco" & median_coefficient <0.5) %>% 
#   group_by(diffexpressed) %>% summarise(n =n())
# 
# 
# 
# write.csv(df_map_cor, "data/glycoproteomics/proteomics_DE_mappedDeglyco_020425.csv", row.names = F)
