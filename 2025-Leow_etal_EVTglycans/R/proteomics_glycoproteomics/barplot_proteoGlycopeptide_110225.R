# Barplot to compare global proteo and glycopeptide DE
# Author: Ke Leow
# Date: 02/07/25

#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)

#--------------------------------
# Load data 
#--------------------------------
proteomics_data <- read_csv('data/glycoproteomics/24-02-27_Angelo_Placenta_ProteomicsData.csv') %>% 
   mutate(`Student's T-test Difference DB_DP` = -`Student's T-test Difference DP_DB`)
glycopeptide_data <- read_csv('data/glycoproteomics/24-12-06_E00006_Angelo_LFQ_DB_DP_PL_Updated_PartialGalNew_021025.csv') %>% 
  mutate(Glycopeptide = paste(Modified.Sequence, composition, sep = "_"))

plot_gene_comparison <- function(gene, proteomics_data, glycopeptide_data) {
  # Extract protein-level data for the gene
  prot_data <- proteomics_data %>% filter(Genes == gene)
  
  # Extract glycopeptide-level data for the gene
  glyco_data <- glycopeptide_data %>% filter(Gene == gene)
  
  # Add a type column
  prot_data <- prot_data %>% mutate(Type = "Protein", Glycopeptide = Genes)
  glyco_data <- glyco_data %>% mutate(Type = "Glycopeptide")
  
  # Combine data
  plot_data <- 
    bind_rows(prot_data %>%
                select(Glycopeptide,
                       `Student's T-test Difference DB_DP`,
                       `Student's T-test q-value DP_DB`,
                       Type) %>%
                rename("logFC" = "Student's T-test Difference DB_DP",
                       "adj.P.Val" = "Student's T-test q-value DP_DB" ), 
              glyco_data %>% select(Glycopeptide, 
                                    treated_DB_vs_DP_limma_logFC, 
                                    treated_DB_vs_DP_limma_adj.P.Val,
                                    Partially.Galactosylated.new) %>% 
                rename("logFC" = "treated_DB_vs_DP_limma_logFC",
                       "adj.P.Val" = "treated_DB_vs_DP_limma_adj.P.Val" )
    )
  
  
  # Define colors based on significance
  plot_data$Significance <- ifelse(plot_data$adj.P.Val < 0.05, "p < 0.05", "p ≥ 0.05")
  
  # Ensure "Protein" appears first in the x-axis order
  plot_data$Glycopeptide <- factor(plot_data$Glycopeptide,
                                   levels = c(unique(glyco_data$Glycopeptide), gene))
  return(plot_data)
}
#--------------------------------
# Plot
#--------------------------------
gene = "LGALS3BP"
pdf("R_plots/glycoproteomics/barplotLGALS3BP_proteoGlycopeptide_logfc.pdf", , width = 6, height = 6)
plot_data <- plot_gene_comparison(gene, proteomics_data, glycopeptide_data)
# Plot
ggplot(plot_data, aes(y = Glycopeptide, x = logFC, fill = Significance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red","grey")) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "black") +  # Dotted lines at x = 1 and -1
  theme_classic() +
  labs(title = gene,
       x = NULL,
       x = "Log2 Fold Change",
       fill = "Partial.Gal") +  # Legend label
  theme(axis.text.y = element_text(hjust = 1))
dev.off()

# gene = "FN1"
# plot_data <- plot_gene_comparison(gene, proteomics_data, glycopeptide_data)
# # Plot
# ggplot(plot_data, aes(y = Glycopeptide, x = logFC, fill = Partially.Galactosylated.new)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("grey", "red")) +
#   geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "black") +  # Dotted lines at x = 1 and -1
#   theme_classic() +
#   labs(title = gene,
#        x = NULL,
#        x = "Log2 Fold Change",
#        fill = "Partial.Gal") +  # Legend label
#   theme(axis.text.y = element_text(hjust = 1))
# 
# gene = "ITGA5"
# plot_data <- plot_gene_comparison(gene, proteomics_data, glycopeptide_data)
# # Plot
# ggplot(plot_data, aes(y = Glycopeptide, x = logFC, fill = Partially.Galactosylated.new)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("grey", "red")) +
#   geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "black") +  # Dotted lines at x = 1 and -1
#   theme_classic() +
#   labs(title = gene,
#        x = NULL,
#        x = "Log2 Fold Change",
#        fill = "Partial.Gal") +  # Legend label
#   theme(axis.text.y = element_text(hjust = 1))
# 
# gene = "LAMP2"
# plot_data <- plot_gene_comparison(gene, proteomics_data, glycopeptide_data)
# # Plot
# ggplot(plot_data, aes(y = Glycopeptide, x = logFC, fill = Partially.Galactosylated.new)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("grey", "red")) +
#   geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "black") +  # Dotted lines at x = 1 and -1
#   theme_classic() +
#   labs(title = gene,
#        x = NULL,
#        x = "Log2 Fold Change",
#        fill = "Partial.Gal") +  # Legend label
#   theme(axis.text.y = element_text(hjust = 1))
# 
