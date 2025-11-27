# Get list of DB and DP glycoproteins from treated Lacnacase data
# Author: Ke Leow
# Date: 05/27/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(VennDiagram)

#--------------------------------
# Load data
#--------------------------------
data <- read_csv('data/glycoproteomics/24-12-06_E00006_Angelo_LFQ_DB_DP_PL_Updated_PartialGalNew_021025.csv') %>%
  mutate(glycan.class = case_when(
    Partially.Galactosylated.new == TRUE ~ "Partially Galactosylated",
    Polylacnac == TRUE ~ "Polylacnac",
    Tetraantennary == TRUE ~ "Tetraantennary",
    Triantennary == TRUE ~ "Triantennary",
    Biantennary == TRUE ~ "Biantennary",
    `Glycan Classification` == "High Mannose" ~ "High Mannose",
    TRUE ~ "Other Complex/Hybrid"
  ))

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
# Get list of lacnac glycopeptides in treated
#--------------------------------
# Filter for DB
# Filter for partial gal, tetra or polylacnac
treated_DB_HB <- data %>% filter(diffexpressed_treated == "DB") %>% 
  filter(glycan.class %in% c("Partially Galactosylated", "Polylacnac", "Tetraantennary")) %>% 
  select(1:8, composition, glycan.class, treated_DB_vs_DP_limma_logFC, treated_DB_vs_DP_limma_P.Value,treated_DB_vs_DP_limma_adj.P.Val)

treated_DB_HB_genes <- treated_DB_HB %>% pull(Gene)%>% unique()

# write.csv(treated_DB_HB, "data/glycoproteomics//DBDP_highlyBranchedGlycopeptides_treated_052725.csv", row.names = F)

#--------------------------------
# Get list of lacnac glycopeptides in untreated
#--------------------------------
# Filter for DB
# Filter for tetra or polylacnac
untreated_DB_HB <- data %>% filter(diffexpressed_untreated == "DB") %>% 
  filter(glycan.class %in% c("Polylacnac", "Tetraantennary")) %>% 
  select(1:8, composition, glycan.class, treated_DB_vs_DP_limma_logFC, treated_DB_vs_DP_limma_P.Value,treated_DB_vs_DP_limma_adj.P.Val)

untreated_DB_HB_genes <- untreated_DB_HB %>% pull(Gene)%>% unique()

# write.csv(untreated_DB_HB, "data/glycoproteomics//DBDP_highlyBranchedGlycopeptides_untreated_052725.csv", row.names = F)

#--------------------------------
# Summarize treated and untreated lists
#--------------------------------
# Union of all genes
all_genes <- union(treated_DB_HB_genes, untreated_DB_HB_genes)

# Determine the group
group <- sapply(all_genes, function(gene) {
  in_treated <- gene %in% treated_DB_HB_genes
  in_untreated <- gene %in% untreated_DB_HB_genes
  if (in_treated & in_untreated) {
    "both"
  } else if (in_treated) {
    "treated"
  } else {
    "untreated"
  }
})

# Create data frame
gene_df <- data.frame(Gene = all_genes, DB_DP_comparison = group) %>% 
  arrange(DB_DP_comparison)


# write.csv(gene_df, "data/glycoproteomics//DBDP_highlyBranchedGenes_treated_untreated_summary_052725.csv", row.names = F)

#--------------------------------
# Plot venn
#--------------------------------
#suppress log if needed
flog.threshold(ERROR, name = "VennDiagramLogger")

#compare with Glyco list
venn.plot <- venn.diagram(
  x = list(treated= treated_DB_HB_genes, untreated = untreated_DB_HB_genes),
  # category.names = c("Lacnacase", 
  #                    "Deglyco"),
  filename = NULL,
  output = TRUE,
  # Customizing fonts to use Arial
  cat.fontfamily = "sans",
  fontfamily = "sans",
  # Increase font size for numbers inside the circles
  cex = 2,  # Adjust as needed (larger number = larger font)
  # Increase font size for category labels
  cat.cex = 2,  # Adjust as needed
  # Adjust label positions to avoid overlap with circles
  cat.pos = c(0, 0),  # Adjust position angles for labels
  # cat.dist = c(0.15, 0.15),  # Increase distance of labels from circles
  scaled = FALSE  # <- disables area-proportional scaling
)

pdf("R_plots/glycoproteomics/venn_dbHB_treated_untreated.pdf", width = 6, height = 6)
# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)
dev.off()
