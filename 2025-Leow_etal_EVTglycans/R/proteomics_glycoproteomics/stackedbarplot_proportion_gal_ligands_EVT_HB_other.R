# map TCGA degs downloaded from GEPIA2
# Author: Ke Leow
# Date: 04/03/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
## 1. Helper: get galectin ligands from one STRING interaction table ----
get_galectin_ligands <- function(interactions_file, galectin_genes) {
  string <- readr::read_tsv(interactions_file, show_col_types = FALSE)
  
  string_long <- string %>%
    dplyr::rename(node1 = `#node1`) %>% 
    dplyr::mutate(
      galectin = ifelse(node1 %in% galectin_genes, node1,
                        ifelse(node2 %in% galectin_genes, node2, NA)),
      gene = ifelse(node1 %in% galectin_genes, node2,
                    ifelse(node2 %in% galectin_genes, node1, NA))
    ) %>%
    dplyr::filter(!is.na(galectin), !is.na(gene)) %>%
    dplyr::filter(!gene %in% galectin_genes)   # drop LGALS themselves
  
  # unique ligand genes
  unique(string_long$gene)
}

## 2. Helper: get denominator (non-LGALS genes) from string_node_degrees.tsv ----
get_total_genes <- function(folder, galectin_genes) {
  node_file <- file.path(folder, "string_node_degrees.tsv")
  nodes <- readr::read_tsv(node_file, show_col_types = FALSE)
  
  # Adjust this if the column name differs (e.g. "preferredName", "#node")
  # For many STRING exports this is "node" or "preferredName"
  if ("node" %in% names(nodes)) {
    gene_col <- "node"
  } else if ("preferredName" %in% names(nodes)) {
    gene_col <- "preferredName"
  } else if ("#node" %in% names(nodes)) {
    gene_col <- "#node"
  } else {
    stop("Could not find a node column in string_node_degrees.tsv")
  }
  
  all_genes <- nodes[[gene_col]]
  # Exclude LGALS genes from denominator
  non_galectin_genes <- setdiff(unique(all_genes), galectin_genes)
  non_galectin_genes
}


#--------------------------------
# Get counts of galectin ligands
#--------------------------------
# Define galectin gene symbols
galectin_genes <- c("LGALS1", "LGALS2", "LGALS3", "LGALS4", 
                    "LGALS7", "LGALS7B", "LGALS8", "LGALS9", "LGALS12", "LGALS13", "LGALS14", "LGALS16")

## 3. Define the four lists and paths ----
interaction_files <- c(
  "EVT-HB"    = "STRING_interactions/EVT_HB_LGALS/string_interactions_short (2).tsv",
  "other" = "STRING_interactions/19genes_other_HB_LGALS/string_interactions_short (2).tsv"
)

## 4. Loop over lists to get counts ----
ratio_df <- purrr::imap_dfr(interaction_files, function(f, list_name) {
  folder <- dirname(f)
  
  # denominator: all non-LGALS genes from string_node_degrees.tsv
  non_galectin_genes <- get_total_genes(folder, galectin_genes)
  n_total <- length(non_galectin_genes)
  
  # numerator: galectin ligands (restricted to non-LGALS & present in this list)
  ligands <- get_galectin_ligands(f, galectin_genes)
  ligands <- intersect(ligands, non_galectin_genes)
  n_ligand <- length(ligands)
  
  tibble(
    list = list_name,
    n_total = n_total,
    n_ligand = n_ligand,
    n_non_ligand = n_total - n_ligand,
    prop_ligand = n_ligand / n_total
  )
})

#--------------------------------
# Plot
#--------------------------------
plot_df <- ratio_df %>%
  tidyr::pivot_longer(
    cols = c(n_ligand, n_non_ligand),
    names_to = "type",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    type = dplyr::recode(
      type,
      n_ligand = "Galectin ligand",
      n_non_ligand = "Non-ligand"
    ),
    list = factor(list, levels = c("other-HB", "EVT-HB"))
  )

pdf("R_plots/glycoproteomics/stackedbar_gal_ligand_EVT_HB_other.pdf", , width = 2, height = 4)
ggplot(plot_df, aes(x = list, y = count, fill = type)) +
  geom_col(position = "fill") +  # stacked to 1 → proportions
  scale_fill_manual(
    values = c(
      "Non-ligand" = "grey70",
      "Galectin ligand" = "#54278F"  # purple
    )
  ) +
  # scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL,
    # title = "Proportion of galectin ligands in each glycoprotein list"
  ) +
  theme_classic(base_size = 12)+
  theme(legend.position = "none", 
        # axis.text.x = element_text(angle = 30, hjust = 1)
        )
dev.off()