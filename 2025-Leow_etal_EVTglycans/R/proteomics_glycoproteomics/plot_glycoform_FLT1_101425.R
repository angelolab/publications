# Plot glycoforms
# Author: Ke Leow
# Date: 10/13/25
#--------------------------------
# Load packages/functions
#--------------------------------
library(ggplot2)
library(tidyverse)
library(ggforce)
library(grid)   # for unit()
#--------------------------------
# Get protein sequence info
#--------------------------------
gene = "FLT1"
seq_length <- 1338

# Domain regions
domains <- tribble(
  ~domain, ~start, ~end,
  "",    26,     747, #ECD
  "",    759,    780, #TM
  "Kinase domain",   827,    1158 #KD
)

# Glycosites
glycosites <- tribble(
  ~pos,  ~label,
  164,  "",
  196,  "", 
  323,  ""
)

# Highlighted subdomains with required angle column
highlight <- tribble(
  ~name,     ~start, ~end,  ~angle,
  "", 151,     214,   0, #D2
  "",     240,    340,   0 #D3 - move slightly right
)

# --- Plot ---
pdf(paste0("R_plots/glycoproteomics/",gene,"_sequence_alldomains.pdf"), width = 6, height = 3)
ggplot() +
  # Protein backbone
  geom_segment(aes(x = 1, xend = seq_length, y = 0, yend = 0),
               color = "grey40", linewidth = 2) +
  
  # Grey domain rectangles
  geom_rect(
    data = domains,
    aes(xmin = start, xmax = end, ymin = -0.07, ymax = 0.07),
    fill = "grey80", color = NA
  ) +
  
  # Domain labels
  geom_text(
    data = domains,
    aes(x = (start + end)/2, y = 0, label = domain),
    color = "black", size = 3
  ) +
  
  # Highlight FNIII-5 and EDB with ovals
  # geom_ellipse(
  #   data = highlight,
  #   aes(x0 = (start + end)/2, y0 = 0,
  #       a = (end - start)/2, b = 0.12,
  #       angle = angle),
  #   fill = "#2c7fb8", color = "white", fill = NA
  # ) +
  geom_rect(
    data = highlight,
    aes(xmin = start, xmax = end, ymin = -0.07, ymax = 0.07),
    fill = "grey40", color = NA, alpha = 0.6
  )+

  # Glycosite lines
  geom_segment(data = glycosites,
               aes(x = pos, xend = pos, y = -0.1, yend = 0.1),
               linewidth = 1, color = "black") +
  
  # Glycosite labels
  geom_text(data = glycosites,
            aes(x = pos, y = 0.13, label = label),
            size = 3.8, vjust = 0) +
  
  # Sequence endpoints
  annotate("text", x = 1, y = -0.1, label = "1", hjust = 0, size = 3.3) +
  annotate("text", x = seq_length, y = -0.1, label = seq_length, hjust = 1, size = 3.3) +
  
  # Style
  # scale_linetype_manual(values = c("FNIII 5" = "solid", "EDB" = "dashed")) +
  coord_cartesian(xlim = c(0, seq_length + 60), ylim = c(-0.55, 0.55), expand = FALSE) +
  labs(title = NULL, fill = NULL) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )
dev.off()

# # Plot FNIII repeat region
# # Parameters
# seq_length <- 2477
# glycosites <- data.frame(
#   pos   = c(1007, 1291),
#   label = c("N1007", "N1291")
# )
# 
# # 18 FNIII repeats, each 90 aa, starting at aa 610
# domains <- tibble::tibble(i = 1:18) |>
#   mutate(start = 610 + (i - 1) * 90,
#          end   = start + 90 - 1)
# 
# # Plot
# pdf("R_plots/glycoproteomics/FN1_sequence.pdf", , width = 6, height = 3)
# ggplot() +
#   # Protein backbone
#   geom_segment(aes(x = 1, xend = seq_length, y = 0, yend = 0),
#                linewidth = 2, color = "grey40", lineend = "round") +
#   # FNIII domain bars
#   geom_rect(data = domains,
#             aes(xmin = start, xmax = end, ymin = -0.05, ymax = 0.05),
#             fill = "grey70",color="white", color = NA) +
# 
#   # Glycosite lines
#   geom_segment(data = glycosites,
#                aes(x = pos, xend = pos, y = -0.1, yend = 0.1),
#                linewidth = 1, color = "black") +
#   # Glycosite labels
#   geom_text(data = glycosites,
#             aes(x = pos, y = 0.13, label = label),
#             size = 3.3, vjust = 0) +
#   # Axis labels for start/end
#   annotate("text", x = 1, y = -0.05, label = "1", hjust = 0, size = 3.3) +
#   annotate("text", x = seq_length, y = -0.05, label = as.character(seq_length),
#            hjust = 1, size = 3.3) +
#   coord_cartesian(xlim = c(0, seq_length + 40), ylim = c(-0.45, 0.55), expand = FALSE, clip = "off") +
#   labs(x = "Amino acid position", y = NULL,
#        # title = "Fibronectin (FN1) sequence with FNIII domains and key N-glycosites"
#        ) +
#   theme_void() +
#   theme(
#     plot.margin = margin(10, 30, 10, 30),
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5)
#   )
# dev.off()

#--------------------------------
# Plot glycan dots
#--------------------------------

glycan_table <- read_csv(paste0("data/glycoproteomics//DB_glycopeptides_",gene,"_pos.csv"))|>
  mutate(glycan_class_relabel = case_when(
    glycan_class %in% c("Partial Gal", "PolyLacNAc", "Tetraantennary") ~ "highly branched",
    glycan_class == "High Mannose" ~ "High Mannose",
    TRUE ~ "Other"
  ))

# === Count per site and ensure all classes present ===
glycan_summary <- glycan_table |>
  count(position, glycan_class_relabel, name = "count") |>
  group_by(position) |>
  mutate(prop = count / sum(count))|>
  mutate(glycan_class_relabel = factor(glycan_class_relabel, levels = c("highly branched", "High Mannose", "Other")))

# === Pie chart per glycosite ===
pdf(paste0("R_plots/glycoproteomics/",gene,"_glycantypes.pdf"), width = 6, height = 3)
ggplot(glycan_summary, aes(x = "", y = prop, fill = glycan_class_relabel)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            size = 8) +
  scale_fill_manual(values = c("High Mannose" = "#2ca02c",
                               "highly branched" = "#FF7518",
                               "Other" = "#0868AC")) +
  facet_wrap(~position, nrow = 1) +
  theme_void(base_size = 14) +
  theme(strip.text.x = element_blank(),   # removes glycosite labels
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 13, face = "bold")
  ) 
dev.off()


# # === Summarize and ensure all glycosites have all 3 glycan classes ===
# glycan_summary <- glycan_table |>
#   count(position, glycan_class_relabel, name = "count") |>
#   complete(position, glycan_class_relabel = c("High Mannose", "highly branched", "Other"),
#            fill = list(count = 0)) |>
#   arrange(position, glycan_class_relabel)
# 
# 
# # === Position glycan dots below each glycosite ===
# glycan_summary <- glycan_summary |>
#   mutate(glycan_class_relabel = factor(glycan_class_relabel,
#                                        levels = c("highly branched","High Mannose",  "Other")),
#          y = -0.6,
#          x_offset = as.numeric(glycan_class_relabel) - mean(as.numeric(glycan_class_relabel)),
#          x_plot = position + x_offset * 40)
# 
# # === Plot ===
# ggplot() +
#   # FNIII domain bars
#   geom_rect(data = domains,
#             aes(xmin = start, xmax = end, ymin = -0.15, ymax = 0.15),
#             fill = "grey70", alpha = 0.9) +
#   # Protein backbone
#   geom_segment(aes(x = 1, xend = seq_length, y = 0, yend = 0),
#                linewidth = 2, color = "grey40", lineend = "round") +
#   # Glycosite lines
#   geom_segment(data = glycosites,
#                aes(x = pos, xend = pos, y = -0.35, yend = 0.35),
#                linewidth = 1, color = "#1f77b4") +
#   # Glycosite labels
#   geom_text(data = glycosites,
#             aes(x = pos, y = 0.4, label = label),
#             size = 3.8, vjust = 0) +
#   # Glycan dots (always three per site)
#   geom_point(data = glycan_summary,
#              aes(x = x_plot, y = y, fill = glycan_class_relabel),
#              size = 5, shape = 21, color = "black") +
#   geom_text(data = glycan_summary,
#             aes(x = x_plot, y = y, label = count),
#             size = 3, color = "black") +
#   scale_fill_manual(values = c("High Mannose" = "#e78ac3",
#                                "highly branched" = "#8da0cb",
#                                "Other" = "#fc8d62")) +
#   # Start and end labels
#   annotate("text", x = 1, y = -0.3, label = "1", hjust = 0, size = 3.3) +
#   annotate("text", x = seq_length, y = -0.3, label = as.character(seq_length),
#            hjust = 1, size = 3.3) +
#   coord_cartesian(xlim = c(0, seq_length + 40), ylim = c(-1.0, 0.55), expand = FALSE) +
#   labs(x = "Amino acid position", y = NULL,
#        title = "Fibronectin (FN1) with FNIII domains, N-glycosites, and glycan composition summary") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title.y  = element_blank(),
#     axis.text.y   = element_blank(),
#     axis.ticks.y  = element_blank(),
#     panel.grid    = element_blank(),
#     legend.title  = element_blank(),
#     legend.position = "bottom"
#   )
