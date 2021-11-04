
###########################################
##..Install packages and open libraries..##
###########################################

source("~/Risom2022_R_scripts/Risom2022_SC_FUNCTIONS")

##..Open necessary libraries..##

library(flowCore)           
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(FlowSOM)
library(plyr)
library(dplyr)
library(viridis)
library(data.table) # to create data table objects
library(wesanderson) # color pallate
library(scales) # for non-scientific notation
library(umap) # umap library
library(forcats) # reverse ordering of factors
library(gdata)
library(ggpubr)
################################# Sam NORMALIZE #######################################################

#setwd("/Volumes/MIBI Files/remote/dataPerCell")
setwd("~/RDATA/DCIS/200428_FRANKENSTEIN5/")

plot.data<-read.csv("200514_CellTable_Fullinfo.csv")


# TUMOR HORMONE CLUST #####################################################################################################
#########################################################################################################################################################
########################################################################################################################################################
#################################..FlowSOM Cluster round 1: Tumor Cells..##########################################################################################################
########################################################################################################################################################
##################################################################################################################################################################################################
##########################################################################################################################################################################################################

tumor<-c("tumor") #can append if wanting to include non-immune clusters that look contaminated
datatumor<-plot.data[plot.data$celllineage %in% tumor,]
# data_stroma<-sc.data.trimbad[sc.data.trimbad$celllineage %in% fibroblast,]
# immunectrls<-c("tonsil","lymphnode")
# data_epi_trim<-droplevels(data_epi[!data_epi$Tissue %in% immunectrls,])

ff_tumor <- flowFrame(exprs = data.matrix(datatumor[datatumor$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"), -c(63:75)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_4=c("ER","AR","HER2")

# sc.data.scale.transform.norm[,4:55] <- data.frame(t(t(sc.data[,4:55]) / as.numeric(percentile.vector)))
# #sc.data.scale.transform.norm<-sc.data.scale.transform.norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)
# 
# #TrimBadPoints
# badpoints<-c(2101,2301,6101,5311,5501)
# sc.data.scale.transform.norm.trimbad<-droplevels(sc.data.scale.transform.norm[!sc.data.scale.transform.norm$Point_Num %in% badpoints,]
##..Run FlowSOM random seed for reproducibility..##

set.seed(121)
out_fSOM_tumor <- FlowSOM::ReadInput(ff_tumor, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_tumor <- FlowSOM::BuildSOM(out_fSOM_tumor, colsToUse = clusterChannels_4, xdim=3, ydim=3)
out_fSOM_tumor <- FlowSOM::BuildMST(out_fSOM_tumor)
labels_tumor <- out_fSOM_tumor$map$mapping[,1]

out_fSOM_tumor <- UpdateNodeSize(out_fSOM_tumor, reset = TRUE)
#tiff("190203_FS-round2-nodes.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_tumor, view = "grid", markers = clusterChannels_4)
#dev.off()
out_fSOM_tumor <- UpdateNodeSize(out_fSOM_tumor)
#tiff("190203_FS-round2-MST.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_tumor, view = "MST", markers = clusterChannels_4)
FlowSOM::PlotStars(out_fSOM_tumor, view="tSNE",markers = clusterChannels_4)
dev.off()

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_tumor<-out_fSOM_tumor[["map"]][["mapping"]]
fs_clusters_tumor<-fs_clusters_tumor[,1]
data_fs_clusters_tumor <- cbind(datatumor, fs_clusters_tumor)

# go through all clusters and calculate mean for every channel
hm_allclusters_tumor <- matrix(, nrow = length(unique(fs_clusters_tumor)), ncol = length(clusterChannels_4))
for(i in 1:length(unique(fs_clusters_tumor))) {
  temp_mat <- data_fs_clusters_tumor[data_fs_clusters_tumor[,"fs_clusters_tumor"] == i, clusterChannels_4]
  hm_allclusters_tumor[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_tumor) <- paste("cluster", 1:length(unique(fs_clusters_tumor)), sep = "")
colnames(hm_allclusters_tumor) <- clusterChannels_4  
hm_allclusters_tumor

# plot heatmap of all clusters
#tiff("plots/draft_figs/190606_FS-myeloid-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/190606_FS-myeloid-hmap.eps",width=15, height=11)
pdf(file = "200514_TumorHormone4", height = 3, width = 4)  
heatmap.2(hm_allclusters_tumor, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
          col = inferno(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))#, margins = c(8,14))
dev.off()
cluster=1
table(data_fs_clusters_tumor$fs_clusters_tumor)
table(data_fs_clusters_tumor[data_fs_clusters_tumor$fs_clusters_tumor==cluster,]$Point_Num)
table(data_fs_clusters_tumor[data_fs_clusters_tumor$fs_clusters_tumor==cluster,]$Tissue)
table(sc.data.trimbad$Tissue)
table(data_fs_clusters_tumor$fs_clusters_tumor)
write.csv(data_fs_clusters_tumor, file="200514_DCIScohort_FLOWSOMtumorHORMONE.csv",row.names = FALSE)
# 
# #### NAME EACH CLUSTER ################################
# ##..Add cell type to each event..##
# TUMOR_clust<-c(1,2,3,4,5,6)
# # MYOEP2_clust<-c(22,18)
# # TUMOR_clust<-c(10,5,1,2,6,9,16,15,11,14,13)
# # TUMOR2_clust<-c(5,10)
# epi_pheno<-fs_clusters_tumor
# epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR_clust,"TUMOR")
# # epi_pheno<-replace(epi_pheno,epi_pheno %in% MYOEP2_clust,"MYOEP2")
# # epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR_clust,"TUMOR")
# # epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR2_clust,"TUMOR2")
# data_fs_clusters_tumor$sublineage<-epi_pheno
# # # data_immune[mast_rows,]$stroma_pheno<-"mast"
# # # data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
# # write.csv(data_fs_clusters_epi, file="200303_DCIScohort_FLOWSOM_tumor16_Annotated.csv",row.names = FALSE)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
ER_clust<-c(1,4)
ER_HER2_clust<-c(2)
ER_AR_clust<-c(3)
AR_HER2_clust<-c(6)
AR_clust<-c(5)
HER2_clust<-c(9,8)
NEG_clust<-c(7)

# MYOEP2_clust<-c(22,18)
# TUMOR_clust<-c(10,5,1,2,6,9,16,15,11,14,13)
# TUMOR2_clust<-c(5,10)
epi_pheno<-fs_clusters_tumor
epi_pheno<-replace(epi_pheno,epi_pheno %in% ER_clust,"ER")
epi_pheno<-replace(epi_pheno,epi_pheno %in% ER_HER2_clust,"ER_HER2")
epi_pheno<-replace(epi_pheno,epi_pheno %in% ER_AR_clust,"ER_AR")
epi_pheno<-replace(epi_pheno,epi_pheno %in% AR_HER2_clust,"AR_HER2")
epi_pheno<-replace(epi_pheno,epi_pheno %in% AR_clust,"AR")
epi_pheno<-replace(epi_pheno,epi_pheno %in% HER2_clust,"HER2")
epi_pheno<-replace(epi_pheno,epi_pheno %in% NEG_clust,"NEG")

data_fs_clusters_tumor$HRpheno<-epi_pheno
# # data_immune[mast_rows,]$stroma_pheno<-"mast"
# # data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_tumor, file="200513_DCIScohort_FLOWSOM_HORMONE_phenotype.csv",row.names = FALSE)

#############################################  ASSESS EPI CLUSTER FREQ BY CONDITIONS  #############################################  

# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

##..Import data..##

datatumor<-read.csv("200513_DCIScohort_FLOWSOM_HORMONE_phenotype.csv")

##..Subset the cell data and sample identifiers..##

cell_data_tumor<-datatumor %>% select(Point_Num, HRpheno)

##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_tumor <- as.data.frame(table(cell_data_tumor$Point_Num, cell_data_tumor$HRpheno))
names(Freqs_tumor) <- c("SampleID","HRpheno","count")

##..Get overall totals on a per sample basis..##

cell_totals_tumor<-aggregate(Freqs_tumor$count, by=list(Category=Freqs_tumor$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_tumor$SampleID)) {
  frequencies <- Freqs_tumor[Freqs_tumor$SampleID==i,"count"] / cell_totals_tumor[cell_totals_tumor$Category==i,2]
  Freqs_tumor[Freqs_tumor$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_tumor<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_tumor$Point_Num)
# filter annotation data by PointNum
annotation_data_tumor_cohort <- droplevels(annotation_data_tumor[annotation_data_tumor$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_tumor_cohort$PointNumber <- factor(annotation_data_tumor_cohort$PointNumber, levels=pointNums)
annotation_data_cohort_tumor<-annotation_data_tumor_cohort[order(annotation_data_tumor_cohort$PointNumber),]
# cast the annotations to the frequency data
anno_data_tumor<-rep(annotation_data_tumor_cohort,1)
annotated_Freqs_tumor<-cbind(Freqs_tumor, anno_data_tumor)



# cell totals per sample
write.csv(cell_totals_tumor, file="200513_HRpheno_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_tumor, file="200513_HRpheno_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_tumor, file="200513_HRpheno_freqs_per_sample_annotated.csv",row.names = FALSE)
# write out single cell data epi

# write out single cell data stroma
# write.csv(data_fs_clusters_stroma, file="SingleCellData_T_annotated.csv",row.names = FALSE)

################################ PLOT Epi Phenos IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)
library(ggpubr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(7)


data_tumor<-read.csv("200513_HRpheno_freqs_per_sample_annotated.csv")
plot_data<-data_tumor
#plot_data<-droplevels(data_epi[data_epi$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_tumor[data_tumor$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
#plot_data<-droplevels(data_epi[data_epi$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$HRpheno),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$HRpheno <- factor(plot_data$HRpheno, levels=cluster_order)
plot_data<-plot_data[order(plot_data$HRpheno),]

##..Plot all clusters across all samples..##
pdf(file = "200513_HRphenoFreqsAllbreast", height =4, width = 7) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(HRpheno),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(HRpheno))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

##..Plot frequency broken down by Point Number..##

pdf(file = "200513_HRphenoAcrossallbreast", height =12, width = 40) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(HRpheno))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 40, hjust=1)) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
cell_bar
dev.off()

##..Plot each cluster frequency as box across conditions..##

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.x = element_text(angle = 20))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = HRpheno)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~HRpheno, scales = "free")
# session

# status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
#   geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(HRpheno))) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~HRpheno, scales = "free")
# status

pdf(file = "200513_HRphenoAcrossCaseCtrl", height =6, width = 6) 
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon','concurrent','contcis','ipscis','continv','ipsinv','normal'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(HRpheno))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~HRpheno, scales = "free")
status
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])

pdf(file = "200513_HRphenoAcrossCaseCtrl_Invasiveonly", height =6, width = 6) 
status<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('tonsil','lymphnode','placenta','colon','concurrent','contcis','ipscis','continv','ipsinv','normal'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(HRpheno))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~HRpheno, scales = "free")
status
dev.off()


# p <- plot_data_invasive[plot_data_invasive$Status %in% c('case','ctrl'),]
# pplot <- ggboxplot(p, x = "Status" , y= "frequency")
# pplot
# pplot + stat_compare_means( method = "t.test")

pdf(file = "200513_HRphenoAcrossOutcome", height =6, width = 8) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(HRpheno))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~HRpheno, scales = "free")
Tissue
dev.off()


#### STATS ##

# my_comparisons <- list( c('normal', 'DCIS'), c('normal', 'NewCIS'), c('normal','NewInv'), c('normal','concurrent'), c('DCIS','concurrent'), c('DCIS','NewInv'))
# 
# TissueBox<-  ggboxplot(plot_data[plot_data$HRpheno %in% c('HRpheno_CK7'),], x = "Tissue", y = "frequency" )  +
#   theme +
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
# TissueBox

##### SAM VIOLIN PLOT ~~~~
# g <- ggplot(h9.trilin.dt.sub) +
#   aes(x = sampleID, y = Puromycin, fill = lineage) + 
#   blank_theme +
#   geom_violin(scale = "width") + 
#   stat_summary( # plots a dot for the median and line for +/-SD
#     fun.data=data_summary, 
#     geom="pointrange", 
#     size = 0.5) +
#   # labs(title = Title) +
#   # ylab("Normalized\nexpression") +
#   geom_hline(yintercept = median(h9.trilin.dt.sub[ # add horizontal line at the median
#     sampleID %in% c("hPSC_Day0_H9P36"),
#     ]$Puromycin),
#     linetype =2) +
#   xlab("") + 
#   stat_compare_means( # perform one way anova
#     method = "anova", 
#     label.y = 1.2, 
#     label.x = 1.2, 
#     size = 2
#   ) + 
#   stat_compare_means(
#     label = "p.signif",
#     hide.ns = F, # perform wilcox test
#     method = "wilcox.test",
#     ref.group = "hPSC_Day0_H9P36",
#     label.y = 1.07,
#     size = 2
#   ) +
#   scale_x_discrete(limits = lineage.order) #+ coord_flip()
# print(g)


