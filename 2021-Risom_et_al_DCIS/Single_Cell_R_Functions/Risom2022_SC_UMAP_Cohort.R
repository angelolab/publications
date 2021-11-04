# UMAP DCIS
# 200403 for DCIS cohort

## Generate umap embedding for single-cell data

#UMAP - generated using the transformed, unscaled parameters. After UMAP generation, All parameters will be scaled 99.9th percentile for visualizaiton.
#(Scaling can be performed before running UMAP - up to you on which approach)

#```{r create umap for each run}

source("~/Risom2022_R_scripts/Risom2022_SC_FUNCTIONS")
library(umap) # umap library
setwd("~/RDATA/DCIS/200813_CLEANTABLE/")

# ALLdata<-read.csv("200909_CellTable_Fullinfo_Status.csv")

#drop other
Alldataplot <- read.csv("201016_CellTable.csv")
dcisconcurrent <- c("DCIS","concurrent")
ALLdata_nonorm <- droplevels(Alldataplot[Alldataplot$Tissue_Type %in% dcisconcurrent, ])

ALLdata_noother <- droplevels(ALLdata_nonorm[!ALLdata_nonorm$celllineage %in% 'other', ])

# i removed hh3 since this should provide much separation
#INCLUDE CD3, CD4, CD8, CD11c, CD14, CD20, CD31, CD36, CD44, CD45, CD68, CK5, CK7, ECAD, FAP, FOXP3, HER2, HH3, HLADRDPDQ, MPO, PanKRT, SMA, VIM, Tryptase, HIF1a, AR, GZMB, P63, CD56, ER, P, Nuc, HH3, cellSize
#DONTINCLUDE COX2, GLUT1, IDO1, Ki67, MMP9, PD1, PDL1, pS6, COLI
cluster.markers <- c(
  
  # Immune
  "CD3",
  "CD4",
  "CD8",
  "CD11c",
  "CD14",
  "CD20",
  # "CD44",
  "CD45",
  "CD68",
  "Tryptase",
  "MPO",
  # "GZMB",
  "VIM",
  "HLADRDPDQ",
  "Nuc",
  # "cell_size",
  
  
  
  # Cancer
  "CK5norm", 
  "CK7",
  "ECAD",
  "PanKRT",
  #"CD44",
  # "P63",
  "ER",
  "AR",
  "HER2",
  # 
  
  # Stroma
  "FAP",
  "SMA",
  "CD36",
  "CD31"
  #"cellSize"
  
)

non.cluster <- c(
  "COX2",
  "GZMB",
  "IDO1",
  "GLUT1",
  
  "MMP9",
  "P63",
  "PD1",
  "PDL1",
  "pS6"
)

surface.markers <- c(
  "CD3",
  "CD4",
  "CD8",
  "CD11c",
  "CD14",
  "CD20",
  "CD31",
  #"CD36",
  "CD44",
  "CD45",
  "CD68"
)

tumor.markers <- c(
  "CK5", 
  "CK7",
  "ECAD",
  "HER2",  
  "PanKRT",
  "SMA",
  "VIM",
  "CD44"
)

non.exp.param <- c(
  "cellLabelInImage",
  "cellSize",
  "X140empty",
  "TiusseNum",
  "Tissue",
  "Background",
  "C",
  "Na",
  "P",
  "Ca40",
  "TotalTissueEv",
  "umap1",
  "umap2",
  "somCluster",
  "metaCluster",
  "metaClusterAnnotation"
  
)


# UMAP is created with unscaled (but transformed) parameters

# sc.data <- create.uMAP.csv(sc.data.scale.transform, cluster.markers) # takes ~1 min 20 seconds
umapHR4 <- uwot::umap(ALLdata_noother[,cluster.markers], 
                   n_neighbors = 30, 
                   min_dist = 0.2, spread = 3, 
                   verbose = TRUE
)

colnames(umapHR4) <- c("umap1", "umap2")
umapDF <- as.data.frame(umapHR4)
plot(umapDF)
# umap
# dev.off()
ALLdata_noother.umap <- cbind(ALLdata_noother, umapDF)
# sc.data.scale <- copy(sc.data)
# sc.data.scale <- scaleData(sc.data.scale, non.exp.param)
ALLdata_noother.umap
params.for.final.fig <- c(
  "PanKRT",
  "CK7",
  "CK5norm",
  "ER",
  "HER2",
  "AR",
  "Ki67",
  "IDO1",
  "PDL1",
  "FAP",
  "CD36",
  "GLUT1"
)


# ALLdata_noother.umap<- read.csv("201017_Cohort_allcellsNoother_UMAP.csv")
#define the tissues you want to umap with
# keep.tissues <- c("normal","case","ctrl","concurrent","ipscis","ipsinv","continv","contcis")
keep.tissues <- c("normal","DCIS","concurrent","IBC","lobular")

ALLdata_noother.umap.goodtissue <- ALLdata_noother.umap[ALLdata_noother.umap$Tissue_Type %in% keep.tissues,]


# write.csv(ALLdata_noother.umap.goodtissue, file="20107_Cohort_allcellsNoother_UMAP.csv",row.names = FALSE)

ALLdata_noother.umap.goodtissue.dt<-as.data.table(ALLdata_noother.umap.goodtissue)

# for(i in 1:length(params.for.final.fig)) {
 # ALLdata_noother.umap.goodtissue[ALLdata_noother.umap.goodtissue$i >1] <- 1
# }


# ALLdata_noother.umap.goodtissue < read.csv("")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# wes<-gg_color_hue(16)


### Print UMAP figures
# {r printing umap figs, echo=F}

# Umap overlaid with tissue identity - ps
print_umap_by_tumor(ALLdata_noother.umap.goodtissue)
ggsave("201017_umap_DCISonly_Cohort_tissue_overlay.ps", height = 6, width = 8, units = "in")

# Umap overlaid with compartment identity - ps
print_umap_by_compartment(ALLdata_noother.umap.goodtissue)
ggsave("201017_umap_Cohort_COMPARTMENT_noLegend_overlay.png", height = 6, width = 8, units = "in")


# Umap overlaid with phenotype identity - png
print_umap_by_sublineage(ALLdata_noother.umap)
ggsave("201017_umap_DCISolny_sublineage_overlay3.png", height = 6, width = 8, units = "in")

# Umap overlaid with phenotype identity - png
print_umap_by_phenotype(ALLdata_noother.umap.goodtissue)
ggsave("umap_phenotype_overlay.png", height = 6, width = 8, units = "in")

# Umap overlaid with lineage identity - png
print_umap_by_lineage(ALLdata_noother.umap.goodtissueHR)
ggsave("umap_lineage_overlay.png", height = 6, width = 8, units = "in")



# Umap overlaid with tissue identity - png
print_umap_by_status(ALLdata_noother.umap.goodtissue)
ggsave("umap_status_overlay.png", height = 6, width = 8, units = "in")



# facet of measured parameters for final fig - png
printUmapFacet(ALLdata_noother.umap.goodtissue.dt, params.for.final.fig, dotSize = 0.01)
ggsave("201017_umap_DCISonly_Cohort_facet_final_FIG1F.png", height = 8, width = 6.5, units = "in")






# facet of measured parameters for final fig - ps
printUmapFacet(sc.data.scale, params.for.final.fig, dotSize = 0.2)
ggsave("umap_facet_final.ps", height = 8, width = 6, units = "in")

# facet of measured parameters for final fig with no legend- ps
printUmapFacet(sc.data.scale, params.for.final.fig,legendKey = F, dotSize = 0.2)
ggsave("umap_facet_final_no_legend.ps", height = 4, width = 6, units = "in")

# facet of measured parameters for final fig with no legend- png
printUmapFacet(sc.data.scale, params.for.final.fig,legendKey = F, dotSize = 0.2)
ggsave("umap_facet_final_no_legend.png", height = 4, width = 5, units = "in")

# printUmapFacet(sc.data.scale, surface.markers)
# ggsave("umap_surface_overlay_scale_facet.ps", height = 5, width = 5, units = "in")
# 
# printUmapFacet(sc.data.scale, non.cluster)
# ggsave("umap_non_cluster_overlay_scale_facet.ps", height = 5, width = 5, units = "in")
# 
# printUmapFacet(sc.data.scale, cluster.markers)
# ggsave("umap_cluster_markers_overlay_scale_facet.ps", height = 5, width = 5, units = "in")
# 
# printUmapFacet(sc.data.scale, tumor.markers)
# ggsave("umap_tumor_markers_overlay_scale_facet.ps", height = 5, width = 5, units = "in")

# printUmapFacet(sc.data, non.cluster)
# printUmapFacet(sc.data.scale, surface.markers)

# print_biaxial_facet(sc.data, "SMA", "CD45", "HH3")

## print graphs

#```{r print annotated cluster graphs, echo=F}

# prints total number of cells for each annotation group - ps
print_cell_event_number_annotation(sc.data)
ggsave("annotated_cluster_cell_events.ps", height = 5, width = 6, units = "in")

# prints total number of cells for each annotation group - png
print_cell_event_number_annotation(sc.data)
ggsave("annotated_cluster_cell_events.png", height = 5, width = 6, units = "in")

# prints umap with color overlay highlighting annotation groups - png
print_umap_by_annotation(sc.data)
ggsave("umap_metacluster_overlay_with_annotation.png", height = 5, width = 7, units = "in")

# prints umap with color overlay highlighting annotation groups - ps
print_umap_by_annotation(sc.data)
ggsave("umap_metacluster_overlay_with_annotation.ps", height = 5, width = 7, units = "in")

# prints facet umap with color overlay highlighting annotation groups - png
print_umap_by_annotation_facet(sc.data)
ggsave("umap_metacluster_overlay_with_annotation_facet.png", height = 5, width = 7, units = "in")

# prints facet umap with color overlay highlighting annotation groups - ps
print_umap_by_annotation_facet(sc.data)
ggsave("umap_metacluster_overlay_with_annotation_facet.ps", height = 5, width = 7, units = "in")


umapDF <- copy(sc.data)
umapDF <- umapDF[metaClusterAnnotation!="none",]
umapDF <- as.data.frame(umapDF)

umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]


ggplot(data = umapDF) +
  aes(x = umap1, y = umap2, color = metaClusterAnnotation) + facet_wrap(~Tissue) +
  geom_point(size = 0.5) +
  theme_minimal() + #geom_jitter() +
  scale_color_discrete() +
  #scale_colour_manual(values =wes) +
  guides(color=guide_legend(nrow = 12, ncol = 1, byrow = T)) + theme(legend.position = "none") + 
  labs(title = "UMAP overlaid with cluster annotation", color = "Annotation\nGroups")
ggsave(filename = "UMAP_cluster_overlay_tissue_facet.ps")


ggplot(sc.data) +
  aes(x = metaClusterAnnotation) +
  geom_bar(fill = "#0c4c8a") +
  theme_minimal() +
  facet_wrap(vars(Tissue)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                                   axis.title.x = element_blank(), legend.key = element_blank()
  )
ggsave(filename = "cluster_event_per_tissue_facet.png")


## print graphs

#```{r print annotated cluster graphs, echo=F}

# prints total number of cells for each annotation group - ps
print_cell_event_number_annotation(sc.data)
ggsave("annotated_cluster_cell_events.ps", height = 5, width = 6, units = "in")

# prints total number of cells for each annotation group - png
print_cell_event_number_annotation(sc.data)
ggsave("annotated_cluster_cell_events.png", height = 5, width = 6, units = "in")

# prints umap with color overlay highlighting annotation groups - png
print_umap_by_annotation(sc.data)
ggsave("umap_metacluster_overlay_with_annotation.png", height = 5, width = 7, units = "in")

# prints umap with color overlay highlighting annotation groups - ps
print_umap_by_annotation(sc.data)
ggsave("umap_metacluster_overlay_with_annotation.ps", height = 5, width = 7, units = "in")

# prints facet umap with color overlay highlighting annotation groups - png
print_umap_by_annotation_facet(sc.data)
ggsave("umap_metacluster_overlay_with_annotation_facet.png", height = 5, width = 7, units = "in")

# prints facet umap with color overlay highlighting annotation groups - ps
print_umap_by_annotation_facet(sc.data)
ggsave("umap_metacluster_overlay_with_annotation_facet.ps", height = 5, width = 7, units = "in")


umapDF <- copy(sc.data)
umapDF <- umapDF[metaClusterAnnotation!="none",]
umapDF <- as.data.frame(umapDF)

umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]


ggplot(data = umapDF) +
  aes(x = umap1, y = umap2, color = metaClusterAnnotation) + facet_wrap(~Tissue) +
  geom_point(size = 0.5) +
  theme_minimal() + #geom_jitter() +
  scale_color_discrete() +
  #scale_colour_manual(values =wes) +
  guides(color=guide_legend(nrow = 12, ncol = 1, byrow = T)) + theme(legend.position = "none") + 
  labs(title = "UMAP overlaid with cluster annotation", color = "Annotation\nGroups")
ggsave(filename = "UMAP_cluster_overlay_tissue_facet.ps")


ggplot(sc.data) +
  aes(x = metaClusterAnnotation) +
  geom_bar(fill = "#0c4c8a") +
  theme_minimal() +
  facet_wrap(vars(Tissue)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                                   axis.title.x = element_blank(), legend.key = element_blank()
  )
ggsave(filename = "cluster_event_per_tissue_facet.png")

```



