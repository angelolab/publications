# Map DE glycoproteins to Andi's transcriptomics data
# Author: Ke Leow
# Date: 02/18/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(VennDiagram)

#--------------------------------
# Load data
#--------------------------------
#get EVT glycoproteins
data <- read_csv('data/Andi_EVT_RoserTranscriptomics/Ke_Roser_unfiltered_deg.csv')

evt_filter <- data %>%
  dplyr::filter(pvals_adj < 0.05, logfoldchanges > 1, cell_type %in% c("pEVT", "eEVT", "iEVT"))  

evt_genes <-evt_filter %>%
  # filter(pvals_adj < 0.1) %>% 
  group_by(gene) %>%
  summarise(
    count = n(),
    cell_types = paste(unique(cell_type), collapse = ", ")
  )

#get EVT lacnacs
lacnac <- read_csv("data/glycoproteomics/DBDP_highlyBranchedGenes_treated_untreated_summary_052725.csv")

evt_lacnac <- lacnac %>% 
  inner_join(evt_genes,., by = join_by(gene == Gene))

# write.csv(evt_lacnac, "data/glycoproteomics/evt_HB_25_052725.csv", row.names = F)

evt_lacnac_genes <- evt_lacnac %>% pull(gene)

#--------------------------------
# Plot venn
#--------------------------------
#suppress log if needed
flog.threshold(ERROR, name = "VennDiagramLogger")

#compare with Glyco list
venn.plot <- venn.diagram(
  x = list(`DB-HB`= lacnac$Gene, EVT = evt_genes$gene),
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

pdf("R_plots/glycoproteomics/venn_EVT_dbHB.pdf", width = 6, height = 6)
# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)
dev.off()