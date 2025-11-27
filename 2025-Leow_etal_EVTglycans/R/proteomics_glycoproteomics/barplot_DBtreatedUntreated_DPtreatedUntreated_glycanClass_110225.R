# Plot percentage of glycan classes for comparisons between treated/untreated DB and treated/untreated DP
# Author: Ke Leow
# Date: 02/10/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(RColorBrewer)
#--------------------------------
# Load data
#--------------------------------
data <- read_csv('data/glycoproteomics/24-12-06_E00006_Angelo_LFQ_DB_DP_PL_Updated_PartialGalNew_021025.csv') %>%
  mutate(glycan.class = case_when(
    Partially.Galactosylated.new == TRUE ~ "Partial Gal",
    Polylacnac == TRUE ~ "PolyLacNAc",
    Tetraantennary == TRUE ~ "Tetraantennary",
    # Triantennary == TRUE ~ "Triantennary",
    # Biantennary == TRUE ~ "Biantennary",
    `Glycan Classification` == "High Mannose" ~ "High Mannose",
    TRUE ~ "Other"
  )) 
data$glycan.class = factor(data$glycan.class, levels = c("Partial Gal", "PolyLacNAc","Tetraantennary","High Mannose","Other"))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
data$diffexpressed_treated <- "NC"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
data$diffexpressed_treated[data$treated_DB_vs_DP_limma_logFC > 1 & data$`treated_DB_vs_DP_limma_adj.P.Val` < 0.05] <- "DB"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data$diffexpressed_treated[data$treated_DB_vs_DP_limma_logFC < -1 & data$`treated_DB_vs_DP_limma_adj.P.Val` < 0.05] <- "DP"
#summarize number of DE features in each group
data %>% group_by(diffexpressed_treated) %>% 
  summarise(DP_DB = n())

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
data$diffexpressed_untreated <- "NC"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
data$diffexpressed_untreated[data$untreated_DB_vs_DP_limma_logFC > 1 & data$`untreated_DB_vs_DP_limma_adj.P.Val` < 0.05] <- "DB"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data$diffexpressed_untreated[data$untreated_DB_vs_DP_limma_logFC < -1 & data$`untreated_DB_vs_DP_limma_adj.P.Val` < 0.05] <- "DP"
#summarize number of DE features in each group
data %>% group_by(diffexpressed_untreated) %>% 
  summarise(DP_DB = n())

#--------------------------------
# Plot ratio of DE glycopeptides by glycan classes - FN1
#--------------------------------
# data_plot <- data %>% 
#   filter(Gene == "FN1") %>% 
  



# Define function to classify regulation
classify_regulation <- function(logFC, pval) {
  if (!is.na(logFC) & !is.na(pval)) {
    if (logFC > 1 & pval < 0.05) return("treated")
    if (logFC < -1 & pval < 0.05) return("untreated")
  }
  return(NA)
}

# Apply classification for treated condition
data <- data %>%
  mutate(
    DB = mapply(classify_regulation, DBtreatedDBuntreated_limma_logFC, DBtreatedDBuntreated_limma_adj.P.Val),
    DP = mapply(classify_regulation, DPtreatedDPuntreated_limma_logFC, DPtreatedDPuntreated_limma_adj.P.Val)
  )

# Count up/downregulated glycopeptides per glycan class and tissue type
counts <- data %>%
  pivot_longer(cols = c(DB, DP), 
               names_to = "Condition", values_to = "Regulation") %>%
  filter(!is.na(Regulation)) %>%
  group_by(Condition, glycan.class, Regulation) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Condition, Regulation) %>%
  mutate(Proportion = Count / sum(Count)) %>% 
  mutate(Condition = ifelse(Condition == "DB", "decB", "decP"))

counts$Regulation <- factor(counts$Regulation, levels = c("untreated", "treated"))
counts$Condition <- factor(counts$Condition, levels = c("decP", "decB"))

# Define custom colors
# custom_colors <- c("Partial Gal"= "#DD4124FF", "PolyLacNAc"= "#0F85A0FF","Tetraantennary"= "#EDD746FF","High Mannose"= "#2ca02c","Other"="grey")
custom_colors <- c("Partial Gal"= "#FF7518", "PolyLacNAc"= brewer.pal(5, "GnBu")[5],"Tetraantennary"= brewer.pal(5, "GnBu")[4],"High Mannose"= "#2ca02c","Other"="grey")

# Plot
pdf("R_plots/glycoproteomics/stackedbar_treated_untreated_glycans_110225.pdf", , width = 4, height = 4)
ggplot(counts, aes(x = Regulation, y = Proportion, fill = glycan.class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~  Condition) +
  scale_fill_manual(values = custom_colors) +
  labs(x = NULL,
       y = "Proportion",
       fill = "Glycan Class") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank()        # Remove ticks
  )
dev.off()
