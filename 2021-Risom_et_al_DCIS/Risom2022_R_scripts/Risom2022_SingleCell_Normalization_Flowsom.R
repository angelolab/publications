# Date created: 200226
# Overview: Script imports concatenated data set of asinh transformed sc data. Runs FlowSOM on the transformed and
# normalized data.


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
################################# Sam NORMALIZE #######################################################
                  
#setwd("/Volumes/MIBI Files/remote/dataPerCell")
setwd("~/Desktop/DCIS/Segmentation/200303_CLEANED_RESEGMENT/")

# project.dir = "~/Desktop/DCIS/Segmentation/"
csv.dir = "~/Single_CSVs/"
info.dir = "~/Info/"
info.csv <- read.csv(file = paste0(info.dir, "METADATA.csv"))


Status <- c(
  "tonsil",
  "lymphnode",
  "normal",
  "ctrl",
  "case",
  "concurrent",
  "ipsinv",
  "continv",
  "ipscis",
  "contcis",
  "placenta",
  "colon"
)

Tissue_Type <- c(
  "Tonsil",
  "Lymph Node",
  "Normal Breast",
  "DCIS",
  "IBC",
  "ILC",
  "Colon",
  "Placenta"
  
)

Tissue <- c(
  "tonsil",
  "lymphnode",
  "normal",
  "DCIS",
  "NewCIS",
  "placenta",
  "colon",
  "concurrent",
  "NewInv"
  
)


table(sc.data$Point_Num)
hist(sc.data$Point_Num)

# add meta info
sc.data <- as.data.table(sc.data)
info.csv <- as.data.table(info.csv)
keep <- c("PointNumber","CohortNumber", "Tissue_Type", "Status", "Tissue", "TMAD_ID","TMAD_Patient","Years_Since_First_CIS","Days_To_First_Event","Case_Identifier","RAHBT_Case_Group") #what columns to keep and merge to the sc.data

sc.data <- merge(sc.data, info.csv[, keep, with=F], by.x="Point_Num", by.y="PointNumber", all.x=T) #merges info.csv and sc.data only where 'x' in sc.data (e.g. Point_Num) matches y in other csv


table(info.csv$Tissue)
table(info.csv$Status)
table(info.csv$Case_Identifier)
hist(sc.data$CohortNumber)

                
##..Save concatenated dataframes..##
write.csv(sc.data, file="SingleCellData.csv",row.names = FALSE)

non.params <- c(
  "event_num",
  "Point_Num",
  "Status",
  "Tissue",
  "Tissue_Type",
  "major_axis_length",
  "minor_axis_length",
  "cell_size",
  "area",
  "label",
  "perimeter",
  "eccentricity",
  "TMAD_ID",
  "TMAD_Patient",
  "Years_Since_First_CIS",
  "Days_To_First_Event",
  "Case_Identifier",
  "RAHBT_Case_Group")

# measure.params <- colnames(sc.data)[!colnames(sc.data) %in% non.params]

########## single cell Linscale and transform // Leave out with combined csv #######

# #take all the column names that are opposite of (!) non.params.
# #all_sizenorm_scale<-all_sizenorm #scaled by 1000 
# #all_sizenorm_scale[,4:50]<-all_sizenorm_scale[,4:50]*1000
# sc.data.scale <- copy(sc.data)
# sc.data.scale[,3:55]<-sc.data.scale[,3:55]*100
# write.csv(sc.data.scale, file="20120_SingleCellData5pxall_SizenormScaled.csv",row.names = FALSE)
# ##### must go to functions and run the asinTranform function written by David before calling it here#####
# #LOAD THE FULL 'IONPATH WP FUNCT TYLER SCS' FUNCTION SET'
# # ^ SEE THE ' source("~/Desktop/Rtools/functions/IONpath_wp_funct_Tyler_SCs_csv.R") ' at the top
# 
# sc.data.scale.transform <- copy(sc.data.scale)
# sc.data.scale.transform <- asinTransform(sc.data.scale.transform, fa=non.params)
# write.csv(sc.data.scale.transform, file="20120_SingleCellData5pxall_SizenormScaledTrans.csv",row.names = FALSE)


###..Percent normalization and scale data 0-99th percentile..##

v <- 1:1000
v
quantile_value <- 0.999
quantile(v, quantile_value)

# calculating percentile for 1 vector
percentile.vector <- apply(sc.data[,4:55], 2, function(x) quantile(x, quantile_value, names = F)) #columns 4-56 have the expression data, leave out area etc
percentile.vector
sc.data.scale.transform.norm<-sc.data
sc.data.scale.transform.norm[,4:55] <- data.frame(t(t(sc.data[,4:55]) / as.numeric(percentile.vector)))
#sc.data.scale.transform.norm<-sc.data.scale.transform.norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)

#TrimBadPoints
badpoints<-c(2101,2301,6101,5311,5501,2329,2205,4411) #6203 gradient,bg
sc.data.scale.transform.norm.trimbad<-droplevels(sc.data.scale.transform.norm[!sc.data.scale.transform.norm$Point_Num %in% badpoints,])
table(sc.data.scale.transform.norm.trimbad$Point_Num)
hist(sc.data.scale.transform.norm.trimbad$CohortNumber)
table(sc.data.scale.transform.norm.trimbad$Tissue)

#write csv
write.csv(sc.data.scale.transform.norm.trimbad, file="200323_Frankenstein4_SizenormLinscaledTransNorm999trimbad.csv",row.names = FALSE) #save the annotated asinh-99th percentile scaled data

####TestPlots###
# g <- ggplot(sc.data.scale.transform, aes(x = CD45, y = PanKRT)) + geom_point(size = 1, alpha = 1/50) + labs(title = "scaled_tranformed")
# print(g)
g <- ggplot(sc.data.scale.transform.norm.trimbad, aes(x = CD45, y = PanKRT)) + geom_point(size = 1, alpha = 1/100) + labs(title = "scaled_transform_norm_trimbad")
print(g)
# g <- ggplot(sc.data.scale, aes(x = CD45, y = PanKRT)) + geom_point(size = 1, alpha = 1/50) + labs(title = "scaled")
# print(g)
g <- ggplot(sc.data, aes(x = CD45, y = PanKRT)) + geom_point(size = 1, alpha = 1/30) + labs(title = "raw")
print(g)


#Look at Nuc distributions
h <- ggplot(sc.data.scale.transform.norm, aes(x =CD4, y = CD8)) + geom_point(size = 1, alpha = 1/20) + labs(title = "NucHH3 vs Area")
print(h)
#Look at HH3 vs P
s <- ggplot(sc.data.scale.transform.norm, aes(x = HH3, y = P)) + geom_point(size = 1, alpha = 1/100) + labs(title = "Area by Session")
print(s)
#Look at HH3 vs Time
s <- ggplot(sc.data.scale.transform.norm, aes(x = CohortNumber, y = Nuc)) + geom_point(size = 1, alpha = 1/100) + labs(title = "NucHH3 by Session")
print(s)
#Look at Carbon vs Session
s <- ggplot(sc.data.scale.transform.norm, aes(x = CohortNumber, y = ECADKRTCD45GLUT1CD44)) + geom_point(size = 1, alpha = 1/100) + labs(title = "Membrane by Session")
print(s)
#Add heat into these scales?? In sc.data you should see the heat is all on axis, avgs are very low lik below 1

################################################### ERIN FLOWSOM ##########################################################
                #######################################################################
                ###..FlowSOM Cluster round 1: Immune, Endo, Epithelial, Fibroblast..###
                #######################################################################

sc.data.scale.transform.norm.trimbad<-read.csv("200323_Frankenstein4_SizenormLinscaledTransNorm999trimbad.csv")


# sc.data.scale.transform.norm.trimbad<-read.csv("200303_Frankenstein2_SizenormLinscaledTransNorm999trimbad.csv")
ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewCIS','NewInv','normal'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column

clusterChannels_1<-c('CD45','SMA','CK7','CK5','VIM','CD31','PanKRT','ECAD','Tryptase','MPO','CD20','CD3','CD8','CD4','CD14','CD68','FAP','CD36','CD11c','HLADRDPDQ','P63','CD44') #add ER, AR, HER2?

##..Run FlowSOM random seed for reproducibility..##

set.seed(781)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusterChannels_1, xdim=10, ydim=10)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
fs_clusters <- out_fSOM$map$mapping[,1]

out_fSOM <- UpdateNodeSize(out_fSOM, reset = TRUE)
FlowSOM::PlotStars(out_fSOM, view = "grid", markers = clusterChannels_1)
devout_fSOM <- UpdateNodeSize(out_fSOM)
FlowSOM::PlotStars(out_fSOM, view = "MST", markers = clusterChannels_1)
FlowSOM::PlotStars(out_fSOM, view="tSNE",markers = clusterChannels_1)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
data_fs_clusters <- as.data.frame(cbind(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewCIS',"NewInv",'normal'),], fs_clusters))

# Print number of cells per cluster..##
table(data_fs_clusters$fs_clusters)
table(data_fs_clusters$fs_clusters, data_fs_clusters$Tissue)
##..Look at tissue comp of cluster..##
cluster <- 100
table(data_fs_clusters[data_fs_clusters$fs_clusters==cluster,]$Point_Num)


hist(data_fs_clusters[data_fs_clusters$fs_clusters==cluster,]$SMA)

# go through all clusters and calculate mean for every channel
hm_allclusters <- matrix( nrow = length(unique(fs_clusters)), ncol = length(clusterChannels_1))
for(i in 1:length(unique(fs_clusters))) {
  temp_mat <- data_fs_clusters[data_fs_clusters[,"fs_clusters"] == i, clusterChannels_1]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- paste("cluster", 1:length(unique(fs_clusters)), sep = "")
colnames(hm_allclusters) <- clusterChannels_1  
hm_allclusters

# plot heatmap of all clusters
#tiff("plots/draft_figs/190605_FlowSOM-lineage-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/190605_FlowSOM-lineage-hmap.eps",width=15, height=11)
pdf(file = "200424_CK7_CK5_FLOWSOMheatmap100clust_Frankenstein4_RESEG_noimmuneaddnewevent_tall", height = 18, width = 12) 
heatmap.2(hm_allclusters, 
           scale = "none",
           Colv = T, Rowv = T,
           hclustfun = hclust,
           dendrogram = c("both","row","column","none"),
           trace = "none",
           #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
           #col = colorRampPalette(brewer.pal(9,"Blues"))(100),
           col = viridis(256),
           density.info = 'none',
           key.title = '',
           lhei = c(1,7),
           cexRow = 0.6, cexCol = 0.9, margins = c(8,14),
          breaks=seq(0, 1, length.out=257))
dev.off() 

# # Pull out sparse pops just in case they are absorbed in clustering
# 
# neut_row_idx<-which(data_fs_clusters$fs_clusters == 105)
# neut_rows<-rownames(data_fs_clusters[neut_row_idx,])
# 
# mast_row_idx<-which(data_fs_clusters$fs_clusters == 22)
# mast_rows<-rownames(data_fs_clusters[mast_row_idx,])
# 

##..Meta-cluster..##

# # try the suggested automatic metaclustering method for a hint for k
# auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 30)
# max(auto_meta)
# 
# # do a manual metaclustering
# chosen_k=20
# set.seed(711)
# out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = chosen_k)
# meta_results <- out_meta[fs_clusters]
# 
# ##..Visualize FlowSOM metaclusters output on heatmap..##
# 
# # make combined expression matrix
# data_meta_clusters30 <- as.data.frame(cbind(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewCIS','NewInv','normal','tonsil','lymphnode'),], meta_results))
# 
# 
# # go through all clusters and calculate mean for every channel
# hm_metaclusters30 <- matrix(, nrow = chosen_k, ncol = length(clusterChannels_1))
# for(i in 1:chosen_k) {
#    temp_mat <- data_meta_clusters30[data_meta_clusters30[,"meta_results"] == i, clusterChannels_1]
#    hm_metaclusters30[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
#  }
# 
# # rename
# rownames(hm_metaclusters30) <- paste("cluster", 1:chosen_k, sep = "")
# colnames(hm_metaclusters30) <- clusterChannels_1
# hm_metaclusters30

# make a metacluster heatmap

#tiff("plots/draft_figs/190531_lineage_metaclust-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
# pdf(file = "200303_Metaclusterheatmap_clustFrankenstein4_20_RESEG", height = 12, width = 16) 
# heatmap.2(hm_metaclusters30,
#            scale = "none",
#            Colv = T, Rowv = T,
#            hclustfun = hclust,
#            dendrogram = c("both","row","column","none"),
#            trace = "none",
#            col = colorRampPalette(brewer.pal(9,"Blues"))(100),
#            #col = magma(256),
#            density.info = 'none',
#            key.title = '',
#            lhei = c(1,7),
#            cexRow = 0.6, cexCol = 0.9, margins = c(8,14),
#            breaks=seq(0, 1, length.out=101))
# dev.off()

# check the amount of cells in each cluster
table(data_fs_clusters$fs_clusters)
table(data_fs_clusters$fs_clusters, data_fs_clusters$Tissue)

pdf(file = "200323_ClustPerTissue", height = 20, width = 20) 
clust_count_pertissue<-as.data.frame(table(data_fs_clusters$fs_clusters, data_fs_clusters$Tissue))
ggplot(data=clust_count_pertissue, aes(x=Var2, y=Freq, fill=Var2)) +
  geom_bar(stat='identity') +
  facet_wrap(~Var1, scale='free_y')
dev.off()
##..Look at tissue comp of cluster..##
cluster <- 13
table(data_fs_clusters[data_fs_clusters$fs_clusters==cluster,]$Tissue)
table(data_fs_clusters[data_fs_clusters$fs_clusters==cluster,]$Point_Num)
hist(data_fs_clusters[data_fs_clusters$fs_clusters==cluster,]$MPO)

#export clusters
write.csv(data_fs_clusters, file="200323_DCIScohortFrankenstein4_ScaleTransNormTrimbad_FLOWSOM100.csv",row.names = FALSE) #save the annotated asinh-99th percentile scaled data
##..Add cell type to each event..##
stroma_clust<-c(90,89,70,59,79,58,47,57,77,48)
endo_clust<-c(67,78,68)
immune_clust<-c(40,24,14,50,69,37,66,28,38,16,6,27,39,29,20,9,10,19,30,17,18,8,26,25,7,35,23,80,4,49,5)
tumor_clust<-c(31,71,41,61,81,72,32,42,83,21,82,91,92,94,33,11,87,12,3,15,1,2,34,88,98,99,95,97,86,75,51,93,22,13,96)
myoep_clust<-c(53,63,52,62,44,54,43,55,56,46,73,84,85,65,64,74,36,45,76)
unknown_clust<-c(100,60)
#stroma_clust<-c(11,15,7,22,19)
cell_type<-fs_clusters

cell_type<-replace(cell_type,cell_type %in% stroma_clust,"stroma")
cell_type<-replace(cell_type,cell_type %in% endo_clust,"endo")
cell_type<-replace(cell_type,cell_type %in% tumor_clust,"tumor")
cell_type<-replace(cell_type,cell_type %in% immune_clust,"immune")
cell_type<-replace(cell_type,cell_type %in% myoep_clust,"myoep")
cell_type<-replace(cell_type,cell_type %in% unknown_clust,"unknown")
#cell_type<-replace(cell_type,cell_type %in% endo_clust,"endothelial")
#cell_type<-replace(cell_type,cell_type %in% fibro_clust,"fibroblast")
#cell_type<-replace(cell_type,cell_type %in% fibro_clust,"myoep")

##..Add cell type to each event..##
data_fs_clusters$lineage<-cell_type
#data_gran_norm$lineage<-data_gran$lineage
write.csv(data_fs_clusters, file="200323_DCIScohortFrankenstein4_ScaleTransNormtrimbad_FLOWSOM100_ANNOTATED.csv",row.names = FALSE) #save the annotated asinh-99th percentile scaled data

 #############################################  ASSESS CLUSTER FREQ BY CONDITIONS  #############################################  
#############################################   #############################################   #############################################  
# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

library(dplyr)
library(viridis)
library(forcats)
library(reshape2)
library(tidyr)
library(mefa)

##..Import data..##

# setwd("~/Desktop/DCIS/Segmentation/200226_CLEANED_Frankenstein/")
data<-read.csv("200323_DCIScohortFrankenstein4_ScaleTransNormtrimbad_FLOWSOM100_ANNOTATED.csv")

##..Subset the cell data and sample identifiers..##

cell_data<-data %>% select(Point_Num, lineage)

##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs <- as.data.frame(table(cell_data$Point_Num, cell_data$lineage))
names(Freqs) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

cell_totals<-aggregate(Freqs$count, by=list(Category=Freqs$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs$SampleID)) {
  frequencies <- Freqs[Freqs$SampleID==i,"count"] / cell_totals[cell_totals$Category==i,2]
  Freqs[Freqs$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data<-read.csv("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/Metadata/200120_COHORT_METADATA.csv") 
annotation_data<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data$Point_Num)
# filter annotation data by PointNum
annotation_data_cohort <- droplevels(annotation_data[annotation_data$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_cohort$PointNumber <- factor(annotation_data_cohort$PointNumber, levels=pointNums)
annotation_data_cohort<-annotation_data_cohort[order(annotation_data_cohort$PointNumber),]
# cast the annotations to the frequency data
anno_data<-rep(annotation_data_cohort,1)
annotated_Freqs<-cbind(Freqs, anno_data)

##..Save as a csv later use and plotting..##

# cell totals per sample
write.csv(cell_totals, file="200323_Frankenstein4_FLOWSOMlineage_cell_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs, file="200323_Frankenstein2_FLOWSOMlineage_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs, file="200323_Frankenstein2_FLOWSOMlineage_freqs_per_sample_annotated.csv",row.names = FALSE)

################################ PLOT METACLUSTERS IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(6)

##..Import data..##
# setwd("~/Desktop/DCIS/Segmentation/Frankenstein5px_single_cell_output/")
# data_meta_clusters30<-read.csv("200303_DCIScohortFrankenstein4_ScaleTransNormtrimbad_FLOWSOMMetacluster20_ANNOTATED.csv")
data<-read.csv("200323_Frankenstein2_FLOWSOMlineage_freqs_per_sample_annotated.csv")
plot_data<-data
plot_data<-droplevels(data[data$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewCIS','NewInv','concurrent'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$cell_type),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$cell_type <- factor(plot_data$cell_type, levels=cluster_order)
plot_data<-plot_data[order(plot_data$cell_type),]

##..Plot all clusters across all samples..##
pdf(file = "200323_Flowsom100_Cellbox", height = 12, width = 16) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(cell_type),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(cell_type))) + 
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
pdf(file = "200323_FLOWSOM100_Cellbar", height = 12, width = 40)
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(cell_type))) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
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
               axis.text.x = element_text(angle = 20, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = cell_type)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~cell_type, scales = "free")
# session
pdf(file = "200323_FLOWSOM100_CaseCtrl", height = 15, width = 15)
status<-ggplot(plot_data[!plot_data$Status %in% c('continv','ipsinv','contcis','ipscis','tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status
dev.off()

pdf(file = "200323_FLowwsom100_FreqNEWEVENT", height = 15, width = 15)
reccur<-ggplot(plot_data[!plot_data$Status %in% c('ctrl', 'case', 'tonsil','lymphnode','placenta','colon','concurrent','normal'),]) +
  geom_boxplot(aes(x=Status, y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
reccur
dev.off()

pdf(file = "200323_FLowwsom100_FreqTissue", height = 15, width = 15)
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS', 'NewCIS', 'concurrent','NewInv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
Tissue
dev.off()
  
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status


#FIBROCLUSTERING ################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
               ################ ################################################################################
                               ##..FlowSOM Cluster round 2:  Stroma cells..##
               ##################### #############################################################################
################ ################################################################################################ ################################################################################################ ####################################################################
################ ################################################################################################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
data_fs_clusters<-read.csv("200323_DCIScohortFrankenstein4_ScaleTransNormtrimbad_FLOWSOM100_ANNOTATED.csv")
library(dplyr)
stroma<-c("stroma") #can append if wanting to include non-immune clusters that look contaminated
data_stroma<-data_fs_clusters[data_fs_clusters$lineage %in% stroma,]
ff_stroma <- flowFrame(exprs = data.matrix(data_stroma[data_stroma$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"),
                                       -c(63:73)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_2=c("VIM","FAP","SMA","CD36")
# clusterChannels_2=c("Tryptase","VIM","CD44","CD45","HLADRDPDQ","CD3","CD20","CD4","CD8", "FAP", "CD11c", "CD36", "MPO", "CD68", "CD31","CD14", "SMA", "GLUT1")
##..Run FlowSOM random seed for reproducibility..##

set.seed(102)
out_fSOM_stroma <- FlowSOM::ReadInput(ff_stroma, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_stroma <- FlowSOM::BuildSOM(out_fSOM_stroma, colsToUse = clusterChannels_2, xdim=3, ydim=2)
out_fSOM_stroma <- FlowSOM::BuildMST(out_fSOM_stroma)
labels_stroma <- out_fSOM_stroma$map$mapping[,1]

out_fSOM_stroma <- UpdateNodeSize(out_fSOM_stroma, reset = TRUE)
#tiff("plots/190203_FS-round2-nodes.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_stroma, view = "grid", markers = clusterChannels_2)
#dev.off()
out_fSOM_stroma <- UpdateNodeSize(out_fSOM_stroma)
#tiff("plots/190203_FS-round2-MST.tiff", units="in", width=15, height=11, res=300)
# pdf(file = "200303_Stroma49_mst", height = 12, width = 12)
FlowSOM::PlotStars(out_fSOM_stroma, view = "MST", markers = clusterChannels_2)
# dev.off()

# pdf(file = "200303_Stroma49_tsne", height = 12, width = 12)
FlowSOM::PlotStars(out_fSOM_stroma, view="tSNE",markers = clusterChannels_2)
 # dev.off()

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_stroma<-out_fSOM_stroma[["map"]][["mapping"]]
fs_clusters_stroma<-fs_clusters_stroma[,1]
data_fs_clusters_stroma <- cbind(data_stroma, fs_clusters_stroma)

# go through all clusters and calculate mean for every channel
hm_allclusters_stroma <- matrix(, nrow = length(unique(fs_clusters_stroma)), ncol = length(clusterChannels_2))
for(i in 1:length(unique(fs_clusters_stroma))) {
  temp_mat <- data_fs_clusters_stroma[data_fs_clusters_stroma[,"fs_clusters_stroma"] == i, clusterChannels_2]
  hm_allclusters_stroma[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_stroma) <- paste("cluster", 1:length(unique(fs_clusters_stroma)), sep = "")
colnames(hm_allclusters_stroma) <- clusterChannels_2  
hm_allclusters_stroma

# plot heatmap of all clusters
#tiff("plots/draft_figs/190603_FS-imm-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/90606_FS-imm-hmap.eps",width=15, height=11)

pdf(file = "200401_Frankenstein4_StromaFLOWSOM6", height =3, width = 4) 
heatmap.2(hm_allclusters_stroma, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
          #col = colorRampPalette(brewer.pal(9,"RdPu"))(100),
          col = inferno(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))#, margins = c(8,14))
dev.off()

table(data_fs_clusters_stroma$fs_clusters_stroma)
table(data_fs_clusters_stroma$fs_clusters_stroma, data_fs_clusters_stroma$Point_Num)

##..Look at tissue comp of cluster..##
cluster <- 3
table(data_fs_clusters_stroma[data_fs_clusters_stroma$fs_clusters_stroma==cluster,]$Point_Num)
hist(data_fs_clusters_stroma[data_fs_clusters_stroma$fs_clusters_stroma==cluster,]$VIM)

write.csv(data_fs_clusters_stroma, file="200323_Fibro6.csv",row.names = FALSE)

# #### Granular Stroma Phenotyupe Clusters ################################
# ##..Add cell type to each event..##

FIBROBLAST_clust<-c(1,2,3,4,5,6)
stroma_sublin<-fs_clusters_stroma
stroma_sublin<-replace(stroma_sublin,stroma_sublin %in% FIBROBLAST_clust,"FIBROBLAST")
# stroma_sublin<-replace(stroma_sublin,stroma_sublin %in% ENDOTHELIUM_clust,"ENDOTHELIUM")
# stroma_sublin<-replace(stroma_sublin,stroma_sublin %in% STROMAOTHER_clust,"STROMAOTHER")
data_fs_clusters_stroma$sublineage<-stroma_sublin

write.csv(data_fs_clusters_stroma, file="2003401_DCIScohort_FLOWSOMstroma5_AnnotatedSublineage.csv",row.names = FALSE)
#
#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##

FIBRO1_clust<-c(5,6)
CAF1_clust<-c(3)
CAF2_clust<-c(4)
CAF3_clust<-c(2)
NORMFIBRO_clust<-c(1)

# ENDO4_clust<-c(89)
stroma_pheno<-fs_clusters_stroma
stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% FIBRO1_clust,"FIBRO1")
stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% CAF1_clust,"CAF1")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% STROMAOTHER2_clust,"STROMAOTHER2")
stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% CAF2_clust,"CAF2")
stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% CAF3_clust,"CAF3")
stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% NORMFIBRO_clust,"NORMFIBRO")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% FIBRO5_clust,"FIBRO5")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% FIBRO6_clust,"FIBRO6")
# # stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% FIBRO7_clust,"FIBRO7")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% ENDO1_clust,"ENDO1")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% ENDO2_clust,"ENDO2")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% ENDO3_clust,"ENDO3")
# stroma_pheno<-replace(stroma_pheno,stroma_pheno %in% ENDO4_clust,"ENDO4")
data_fs_clusters_stroma$phenotype<-stroma_pheno
# data_ALL[mast_rows,]$stroma_pheno<-"mast"
# data_ALL[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_stroma, file="200401_DCIScohort_FLOWSOM_Fibro6_phenotype.csv",row.names = FALSE)

table(data_fs_clusters_stroma$Point_Num)
# ################### DROP IMMUNE CTRLS ############################
# immunectrls<-c("tonsil","lymphnode")
# data_stromapheno_trim<-droplevels(data_fs_clusters_stroma[!data_fs_clusters_stroma$Tissue %in% immunectrls,])
# write.csv(data_stromapheno_trim, file="200304_DCIScohort_FLOWSOMstroma100_AnnotatedPhenotype_Dropimmunectrl.csv",row.names = FALSE)
# #############################################  ASSESS IMMUNE CLUSTER FREQ BY CONDITIONS  #############################################  

# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

##..Import data..##

data<-read.csv("200401_DCIScohort_FLOWSOM_Fibro6_phenotype.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_stroma<-data %>% select(Point_Num, fs_clusters_stroma)
cell_data_stroma<-data %>% select(Point_Num, phenotype)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_stroma <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$phenotype))
# Freqs_stroma <- as.data.frame(table(cell_data_stroma$Point_Num, cell_data_stroma$fs_clusters_stroma))
names(Freqs_stroma) <- c("SampleID","phenotype","count")

##..Get overall totals on a per sample basis..##

cell_totals_stroma<-aggregate(Freqs_stroma$count, by=list(Category=Freqs_stroma$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_stroma$SampleID)) {
  frequencies <- Freqs_stroma[Freqs_stroma$SampleID==i,"count"] / cell_totals_stroma[cell_totals_stroma$Category==i,2]
  Freqs_stroma[Freqs_stroma$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_stroma<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_stroma<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_stroma$Point_Num)
# filter annotation data by PointNum
annotation_data_stroma_cohort <- droplevels(annotation_data_stroma[annotation_data_stroma$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_stroma_cohort$PointNumber <- factor(annotation_data_stroma_cohort$PointNumber, levels=pointNums)
annotation_data_cohort<-annotation_data_stroma_cohort[order(annotation_data_stroma_cohort$PointNumber),]
# cast the annotations to the frequency data
# anno_data_stroma<-rep(annotation_data_stroma_cohort,length(unique(cell_data_stroma$fs_clusters_stroma))
anno_data_stroma<-rep(annotation_data_stroma_cohort,1)  #the rep of #phenotypes was not good!
annotated_Freqs_stroma<-cbind(Freqs_stroma, anno_data_stroma)



# cell totals per sample
write.csv(cell_totals_stroma, file="Stroma_pheno_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_stroma, file="Stroma_pheno_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_stroma, file="Stroma_pheno_freqs_per_sample_annotated.csv",row.names = FALSE)

read.csv("Stroma_pheno_totals.csv")
data <- read.csv("Stroma_pheno_freqs_per_sample.csv")
table(data$cell_type)
table



################################ PLOT Stroma Phenos IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(6)


data_stroma<-read.csv("Stroma_pheno_freqs_per_sample_annotated.csv")
plot_data<-data_stroma
#plot_data<-droplevels(data_stroma[data_stroma$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_stroma[data_stroma$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
# plot_data<-droplevels(data_stroma[data_stroma$Status %in% c("case", "ctrl", "concurrent", "ipsinv", "continv", "contcis", "ipscis", "normal"),])
#plot_data<-droplevels(data_stroma[data_stroma$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$phenotype),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$cell_type <- factor(plot_data$phenotype, levels=cluster_order)
plot_data<-plot_data[order(plot_data$phenotype),]

##..Plot all clusters across all samples..##
pdf(file = "200401_StromaFreqsAllBreast", height =5, width =10) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(cell_type),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(cell_type))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 40, hjust=1)) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

##..Plot frequency broken down by Point Number..##

pdf(file = "200401_StromaClustAcrossPatients", height =12, width = 35) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(cell_type))) + 
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
               axis.text.x = element_text(angle = 40, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = cell_type)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~cell_type, scales = "free")
# session

pdf(file = "200401_StromaTypeAcrossDCISoutcome", height =6, width = 7) 
status<-ggplot(plot_data[!plot_data$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(cell_type))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])



pdf(file = "200401_StromaTypeAcrossDCISoutcome_InvasiveOnly", height =6, width = 7) 
statusinvasive<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(cell_type))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
statusinvasive
dev.off()

library(ggsignif)

pdf(file = "200401_StromaTypeAcrossOutcome", height =6, width = 7) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
Tissue
dev.off()

# ENDO CLUSTERING################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################
##..FlowSOM Cluster round 2:  Endo cells..##
##################### #############################################################################
################ ################################################################################################ ################################################################################################ ####################################################################
################ ################################################################################################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
################ ################################################################################################ ################################################################################
data_fs_clusters<-read.csv("200323_DCIScohortFrankenstein4_ScaleTransNormtrimbad_FLOWSOM100_ANNOTATED.csv")
library(dplyr)
endo<-c("endo") #can append if wanting to include non-immune clusters that look contaminated
data_endo<-data_fs_clusters[data_fs_clusters$lineage %in% endo,]
ff_endo <- flowFrame(exprs = data.matrix(data_endo[data_endo$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"),
                                                       -c(63:73)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_2=c("CD31","FAP","SMA","CD36")
# clusterChannels_2=c("Tryptase","VIM","CD44","CD45","HLADRDPDQ","CD3","CD20","CD4","CD8", "FAP", "CD11c", "CD36", "MPO", "CD68", "CD31","CD14", "SMA", "GLUT1")
##..Run FlowSOM random seed for reproducibility..##

set.seed(106)
out_fSOM_endo <- FlowSOM::ReadInput(ff_endo, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_endo <- FlowSOM::BuildSOM(out_fSOM_endo, colsToUse = clusterChannels_2, xdim=2, ydim=2)
out_fSOM_endo <- FlowSOM::BuildMST(out_fSOM_endo)
labels_endo <- out_fSOM_endo$map$mapping[,1]

out_fSOM_endo <- UpdateNodeSize(out_fSOM_endo, reset = TRUE)
#tiff("plots/190203_FS-round2-nodes.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_endo, view = "grid", markers = clusterChannels_2)
#dev.off()
out_fSOM_endo <- UpdateNodeSize(out_fSOM_endo)
#tiff("plots/190203_FS-round2-MST.tiff", units="in", width=15, height=11, res=300)
# pdf(file = "200303_endo49_mst", height = 12, width = 12)
FlowSOM::PlotStars(out_fSOM_endo, view = "MST", markers = clusterChannels_2)
# dev.off()

# pdf(file = "200303_endo49_tsne", height = 12, width = 12)
FlowSOM::PlotStars(out_fSOM_endo, view="tSNE",markers = clusterChannels_2)
# dev.off()

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_endo<-out_fSOM_endo[["map"]][["mapping"]]
fs_clusters_endo<-fs_clusters_endo[,1]
data_fs_clusters_endo <- cbind(data_endo, fs_clusters_endo)

# go through all clusters and calculate mean for every channel
hm_allclusters_endo <- matrix(, nrow = length(unique(fs_clusters_endo)), ncol = length(clusterChannels_2))
for(i in 1:length(unique(fs_clusters_endo))) {
  temp_mat <- data_fs_clusters_endo[data_fs_clusters_endo[,"fs_clusters_endo"] == i, clusterChannels_2]
  hm_allclusters_endo[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_endo) <- paste("cluster", 1:length(unique(fs_clusters_endo)), sep = "")
colnames(hm_allclusters_endo) <- clusterChannels_2  
hm_allclusters_endo

# plot heatmap of all clusters
#tiff("plots/draft_figs/190603_FS-imm-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/90606_FS-imm-hmap.eps",width=15, height=11)

pdf(file = "200401_Frankenstein4_endoFLOWSOM4", height =3, width = 4) 
heatmap.2(hm_allclusters_endo, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
          #col = colorRampPalette(brewer.pal(9,"RdPu"))(100),
          col = inferno(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))#, margins = c(8,14))
dev.off()
table(data_fs_clusters_endo$fs_clusters_endo)
table(data_fs_clusters_endo$fs_clusters_endo, data_fs_clusters_endo$Point_Num)

##..Look at tissue comp of cluster..##
cluster <- 3
table(data_fs_clusters_endo[data_fs_clusters_endo$fs_clusters_endo==cluster,]$Point_Num)
hist(data_fs_clusters_endo[data_fs_clusters_endo$fs_clusters_endo==cluster,]$VIM)

write.csv(data_fs_clusters_endo, file="200401_Endo4.csv",row.names = FALSE)

# #### Granular endo Phenotyupe Clusters ################################
# ##..Add cell type to each event..##

ENDO_clust<-c(1,2,3,4)
endo_sublin<-fs_clusters_endo
endo_sublin<-replace(endo_sublin,endo_sublin %in% ENDO_clust,"ENDO")
# endo_sublin<-replace(endo_sublin,endo_sublin %in% ENDOTHELIUM_clust,"ENDOTHELIUM")
# endo_sublin<-replace(endo_sublin,endo_sublin %in% endoOTHER_clust,"endoOTHER")
data_fs_clusters_endo$sublineage<-endo_sublin

write.csv(data_fs_clusters_endo, file="200401_DCIScohort_FLOWSOM_Endo4_AnnotatedSublineage.csv",row.names = FALSE)
#
#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##

ENDO1_clust<-c(1)
ENDO2_clust<-c(3)
ENDO3_clust<-c(4)
ENDO4_clust<-c(2)


# ENDO4_clust<-c(89)
endo_pheno<-fs_clusters_endo
endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO1_clust,"ENDO1")
endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO2_clust,"ENDO2")
endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO3_clust,"ENDO3")
endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO4_clust,"ENDO4")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% CAF1_clust,"CAF1")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% endoOTHER2_clust,"endoOTHER2")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% CAF2_clust,"CAF2")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% CAF3_clust,"CAF3")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% NORMFIBRO_clust,"NORMFIBRO")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% FIBRO5_clust,"FIBRO5")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% FIBRO6_clust,"FIBRO6")
# # endo_pheno<-replace(endo_pheno,endo_pheno %in% FIBRO7_clust,"FIBRO7")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO1_clust,"ENDO1")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO2_clust,"ENDO2")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO3_clust,"ENDO3")
# endo_pheno<-replace(endo_pheno,endo_pheno %in% ENDO4_clust,"ENDO4")
data_fs_clusters_endo$phenotype<-endo_pheno
# data_ALL[mast_rows,]$endo_pheno<-"mast"
# data_ALL[neut_rows,]$endo_pheno<-"neutrophil"
write.csv(data_fs_clusters_endo, file="200401_DCIScohort_FLOWSOM_Endo4_phenotype.csv",row.names = FALSE)


# ################### DROP IMMUNE CTRLS ############################
# immunectrls<-c("tonsil","lymphnode")
# data_endopheno_trim<-droplevels(data_fs_clusters_endo[!data_fs_clusters_endo$Tissue %in% immunectrls,])
# write.csv(data_endopheno_trim, file="200304_DCIScohort_FLOWSOMendo100_AnnotatedPhenotype_Dropimmunectrl.csv",row.names = FALSE)
# #############################################  ASSESS IMMUNE CLUSTER FREQ BY CONDITIONS  #############################################  

# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

##..Import data..##

data<-read.csv("200401_DCIScohort_FLOWSOM_Endo4_phenotype.csv")

##..Subset the cell data and sample identifiers..##

# cell_data_endo<-data %>% select(Point_Num, fs_clusters_endo)
cell_data_endo<-data %>% select(Point_Num, phenotype)
##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_endo <- as.data.frame(table(cell_data_endo$Point_Num, cell_data_endo$phenotype))
# Freqs_endo <- as.data.frame(table(cell_data_endo$Point_Num, cell_data_endo$fs_clusters_endo))
names(Freqs_endo) <- c("SampleID","phenotype","count")

##..Get overall totals on a per sample basis..##

cell_totals_endo<-aggregate(Freqs_endo$count, by=list(Category=Freqs_endo$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_endo$SampleID)) {
  frequencies <- Freqs_endo[Freqs_endo$SampleID==i,"count"] / cell_totals_endo[cell_totals_endo$Category==i,2]
  Freqs_endo[Freqs_endo$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_endo<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_endo<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_endo$Point_Num)
# filter annotation data by PointNum
annotation_data_endo_cohort <- droplevels(annotation_data_endo[annotation_data_endo$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_endo_cohort$PointNumber <- factor(annotation_data_endo_cohort$PointNumber, levels=pointNums)
annotation_data_cohort<-annotation_data_endo_cohort[order(annotation_data_endo_cohort$PointNumber),]
# cast the annotations to the frequency data
# anno_data_endo<-rep(annotation_data_endo_cohort,length(unique(cell_data_endo$fs_clusters_endo))
anno_data_endo<-rep(annotation_data_endo_cohort,1)
annotated_Freqs_endo<-cbind(Freqs_endo, anno_data_endo)



# cell totals per sample
write.csv(cell_totals_endo, file="endo_pheno_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_endo, file="endo_pheno_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_endo, file="endo_pheno_freqs_per_sample_annotated.csv",row.names = FALSE)



################################ PLOT ENDO Phenos IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(6)


data_endo<-read.csv("endo_pheno_freqs_per_sample_annotated.csv")
plot_data<-data_endo
#plot_data<-droplevels(data_endo[data_endo$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_endo[data_endo$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
# plot_data<-droplevels(data_endo[data_endo$Status %in% c("case", "ctrl", "concurrent", "ipsinv", "continv", "contcis", "ipscis", "normal"),])
#plot_data<-droplevels(data_endo[data_endo$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$phenotype),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$cell_type <- factor(plot_data$phenotype, levels=cluster_order)
plot_data<-plot_data[order(plot_data$phenotype),]

##..Plot all clusters across all samples..##
pdf(file = "200401_endoFreqsAllBreast", height =5, width =10) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(phenotype),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(phenotype))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 40, hjust=1)) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

##..Plot frequency broken down by Point Number..##

pdf(file = "200401_endoClustAcrossPatients", height =12, width = 35) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(phenotype))) + 
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
               axis.text.x = element_text(angle = 40, hjust=1))

# change x to be whatever you want to compare across

# session<-ggplot(plot_data) +
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = cell_type)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~cell_type, scales = "free")
# session

pdf(file = "200401_endoTypeAcrossDCISoutcome", height =6, width = 7) 
status<-ggplot(plot_data[!plot_data$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])

pdf(file = "200401_endoTypeAcrossDCISoutcome_Invasive", height =6, width = 7) 
status<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('normal','lymphnode','placenta','colon','tonsil','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  # geom_signif(comparisons = list(c('ctrl','case')), 
  # map_signif_level=TRUE) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

library(ggsignif)

pdf(file = "200401_endoTypeAcrossOutcome", height =6, width = 7) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
Tissue
dev.off()



# TUMOR CLUST #####################################################################################################
#########################################################################################################################################################
########################################################################################################################################################
#################################..FlowSOM Cluster round 3: Tumor Cells..##########################################################################################################
########################################################################################################################################################
##################################################################################################################################################################################################
##########################################################################################################################################################################################################

tumor<-c("tumor") #can append if wanting to include non-immune clusters that look contaminated
datatumor<-data_fs_clusters[data_fs_clusters$lineage %in% tumor,]
# immunectrls<-c("tonsil","lymphnode")
# data_epi_trim<-droplevels(data_epi[!data_epi$Tissue %in% immunectrls,])

ff_tumor <- flowFrame(exprs = data.matrix(datatumor[datatumor$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"), -c(63:75)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_4=c("CK7","PanKRT", "CK5","ECAD","VIM","CD44")

# sc.data.scale.transform.norm[,4:55] <- data.frame(t(t(sc.data[,4:55]) / as.numeric(percentile.vector)))
# #sc.data.scale.transform.norm<-sc.data.scale.transform.norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)
# 
# #TrimBadPoints
# badpoints<-c(2101,2301,6101,5311,5501)
# sc.data.scale.transform.norm.trimbad<-droplevels(sc.data.scale.transform.norm[!sc.data.scale.transform.norm$Point_Num %in% badpoints,]
##..Run FlowSOM random seed for reproducibility..##

set.seed(669)
out_fSOM_tumor <- FlowSOM::ReadInput(ff_tumor, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_tumor <- FlowSOM::BuildSOM(out_fSOM_tumor, colsToUse = clusterChannels_4, xdim=3, ydim=2)
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
#dev.off()

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
pdf(file = "200401_tumorFLOWSOM6", height = 3, width = 4)  
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
table(data_fs_clusters_tumor$fs_clusters_tumor)
write.csv(data_fs_clusters_tumor, file="200313_DCIScohort_FLOWSOMtumor.csv",row.names = FALSE)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
TUMOR_clust<-c(1,2,3,4,5,6)
# MYOEP2_clust<-c(22,18)
# TUMOR_clust<-c(10,5,1,2,6,9,16,15,11,14,13)
# TUMOR2_clust<-c(5,10)
epi_pheno<-fs_clusters_tumor
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR_clust,"TUMOR")
# epi_pheno<-replace(epi_pheno,epi_pheno %in% MYOEP2_clust,"MYOEP2")
# epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR_clust,"TUMOR")
# epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR2_clust,"TUMOR2")
data_fs_clusters_tumor$sublineage<-epi_pheno
# # data_immune[mast_rows,]$stroma_pheno<-"mast"
# # data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
# write.csv(data_fs_clusters_epi, file="200303_DCIScohort_FLOWSOM_tumor16_Annotated.csv",row.names = FALSE)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
TUMOR1_clust<-c(4)
TUMOR2_clust<-c(5)
TUMOR3_clust<-c(2)
TUMOR4_clust<-c(1)
TUMOR5_clust<-c(6)
TUMOR6_clust<-c(3)
# MYOEP2_clust<-c(22,18)
# TUMOR_clust<-c(10,5,1,2,6,9,16,15,11,14,13)
# TUMOR2_clust<-c(5,10)
epi_pheno<-fs_clusters_tumor
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR1_clust,"TUMOR1")
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR2_clust,"TUMOR2")
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR3_clust,"TUMOR3")
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR4_clust,"TUMOR4")
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR5_clust,"TUMOR5")
epi_pheno<-replace(epi_pheno,epi_pheno %in% TUMOR6_clust,"TUMOR6")

data_fs_clusters_tumor$phenotype<-epi_pheno
# # data_immune[mast_rows,]$stroma_pheno<-"mast"
# # data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_tumor, file="200401_DCIScohort_FLOWSOM_tumor6_phenotype.csv",row.names = FALSE)

#############################################  ASSESS EPI CLUSTER FREQ BY CONDITIONS  #############################################  

# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

##..Import data..##

datatumor<-read.csv("200401_DCIScohort_FLOWSOM_tumor6_phenotype.csv")

##..Subset the cell data and sample identifiers..##

cell_data_tumor<-datatumor %>% select(Point_Num, phenotype)

##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_tumor <- as.data.frame(table(cell_data_tumor$Point_Num, cell_data_tumor$phenotype))
names(Freqs_tumor) <- c("SampleID","phenotype","count")

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
write.csv(cell_totals_tumor, file="Tumor_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_tumor, file="Tumor_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_tumor, file="Tumor_freqs_per_sample_annotated.csv",row.names = FALSE)
# write out single cell data epi
write.csv(data_fs_clusters_tumor, file="SingleCellData_Tumor_annotated.csv",row.names = FALSE)
# write out single cell data stroma
# write.csv(data_fs_clusters_stroma, file="SingleCellData_T_annotated.csv",row.names = FALSE)

################################ PLOT Epi Phenos IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(6)


data_tumor<-read.csv("Tumor_freqs_per_sample_annotated.csv")
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
cluster_order<-levels(fct_reorder(as.factor(plot_data$phenotype),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$phenotype <- factor(plot_data$phenotype, levels=cluster_order)
plot_data<-plot_data[order(plot_data$phenotype),]

##..Plot all clusters across all samples..##
pdf(file = "200401_TumorFreqsAllbreast", height =8, width = 10) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(phenotype),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(phenotype))) + 
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

pdf(file = "200401_TumorAcrossallbreast", height =12, width = 40) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(phenotype))) + 
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
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = phenotype)) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~phenotype, scales = "free")
# session

# status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
#   geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(phenotype))) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~phenotype, scales = "free")
# status

pdf(file = "200401_TumorAcrossCaseCtrl", height =6, width = 7) 
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon','concurrent','contcis','ipscis','continv','ipsinv','normal'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])

pdf(file = "200401_TumorAcrossCaseCtrl_Invasiveonly", height =6, width = 7) 
status<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('tonsil','lymphnode','placenta','colon','concurrent','contcis','ipscis','continv','ipsinv','normal'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

pdf(file = "200401_TumorAcrossOutcome", height =6, width = 7) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
Tissue
dev.off()

#MYOEP CLUST ####################################################################################################################################################################################################################################
##################################################################################################################
########################################################################################################################################################
#############################################################################################################################################################################################
##################################################################################################################
#####################################################################################################
#########################################################################################################################################################
########################################################################################################################################################
#################################..FlowSOM Cluster round 4: MYOEP Cells..##########################################################################################################
########################################################################################################################################################
##################################################################################################################################################################################################
##########################################################################################################################################################################################################

myoep<-c("myoep") #can append if wanting to include non-immune clusters that look contaminated
data_myoep<-data_fs_clusters[data_fs_clusters$lineage %in% myoep,]
# immunectrls<-c("tonsil","lymphnode")
# data_epi_trim<-droplevels(data_epi[!data_epi$Tissue %in% immunectrls,])

ff_myoep <- flowFrame(exprs = data.matrix(data_myoep[data_myoep$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"), -c(63:75)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_5=c("CK7", "CK5","VIM","CD44","P63","SMA")

# sc.data.scale.transform.norm[,4:55] <- data.frame(t(t(sc.data[,4:55]) / as.numeric(percentile.vector)))
# #sc.data.scale.transform.norm<-sc.data.scale.transform.norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)
# 
# #TrimBadPoints
# badpoints<-c(2101,2301,6101,5311,5501)
# sc.data.scale.transform.norm.trimbad<-droplevels(sc.data.scale.transform.norm[!sc.data.scale.transform.norm$Point_Num %in% badpoints,]
##..Run FlowSOM random seed for reproducibility..##

set.seed(677)
out_fSOM_myoep <- FlowSOM::ReadInput(ff_myoep, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_myoep <- FlowSOM::BuildSOM(out_fSOM_myoep, colsToUse = clusterChannels_5, xdim=3, ydim=2)
out_fSOM_myoep <- FlowSOM::BuildMST(out_fSOM_myoep)
labels_myoep <- out_fSOM_myoep$map$mapping[,1]

out_fSOM_myoep <- UpdateNodeSize(out_fSOM_myoep, reset = TRUE)
#tiff("190203_FS-round2-nodes.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_myoep, view = "grid", markers = clusterChannels_5)
#dev.off()
out_fSOM_myoep <- UpdateNodeSize(out_fSOM_myoep)
#tiff("190203_FS-round2-MST.tiff", units="in", width=15, height=11, res=300)
# FlowSOM::PlotStars(out_fSOM_myoep, view = "MST", markers = clusterChannels_5)
# FlowSOM::PlotStars(out_fSOM_myoep, view="tSNE",markers = clusterChannels_5)
#dev.off()

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_myoep<-out_fSOM_myoep[["map"]][["mapping"]]
fs_clusters_myoep<-fs_clusters_myoep[,1]
data_fs_clusters_myoep <- cbind(data_myoep, fs_clusters_myoep)

# go through all clusters and calculate mean for every channel
hm_allclusters_myoep <- matrix(, nrow = length(unique(fs_clusters_myoep)), ncol = length(clusterChannels_5))
for(i in 1:length(unique(fs_clusters_myoep))) {
  temp_mat <- data_fs_clusters_myoep[data_fs_clusters_myoep[,"fs_clusters_myoep"] == i, clusterChannels_5]
  hm_allclusters_myoep[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_myoep) <- paste("cluster", 1:length(unique(fs_clusters_myoep)), sep = "")
colnames(hm_allclusters_myoep) <- clusterChannels_5  
hm_allclusters_myoep

# plot heatmap of all clusters
#tiff("plots/draft_figs/190606_FS-myeloid-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/190606_FS-myeloid-hmap.eps",width=15, height=11)
pdf(file = "200401_myoepFLOWSOM6", height = 3, width = 4)  
heatmap.2(hm_allclusters_myoep, 
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
table(data_fs_clusters_myoep$fs_clusters_myoep)
write.csv(data_fs_clusters_myoep, file="200401_DCIScohort_FLOWSOMmyoep.csv",row.names = FALSE)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
MYOEP_clust<-c(1,2,3,4,5,6)
# MYOEP2_clust<-c(22,18)
myoepi_pheno<-fs_clusters_myoep
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP_clust,"MYOEP")
# myoepi_pheno<-replace(epi_pheno,epi_pheno %in% MYOEP2_clust,"MYOEP2")
data_fs_clusters_myoep$sublineage<-myoepi_pheno
# # data_immune[mast_rows,]$stroma_pheno<-"mast"
# # data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
# write.csv(data_fs_clusters_epi, file="200303_DCIScohort_FLOWSOM_tumor16_Annotated.csv",row.names = FALSE)

#### NAME EACH Phenotype ################################
##..Add cell type to each event..##
MYOEP1_clust<-c(3)
MYOEP2_clust<-c(6)
MYOEP3_clust<-c(1)
MYOEP4_clust<-c(4)
MYOEP5_clust<-c(5)
MYOEP6_clust<-c(2)

cluster=4
table(data_fs_clusters_myoep$fs_clusters_myoep)
table(data_fs_clusters_myoep[data_fs_clusters_myoep$fs_clusters_myoep==cluster,]$Point_Num)
# MYOEP2_clust<-c(22,18)
myoepi_pheno<-fs_clusters_myoep
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP1_clust,"MYOEP1")
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP2_clust,"MYOEP2")
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP3_clust,"MYOEP3")
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP4_clust,"MYOEP4")
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP5_clust,"MYOEP5")
myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP6_clust,"MYOEP6")
# myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP7_clust,"MYOEP7")
# myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP8_clust,"MYOEP8")
# myoepi_pheno<-replace(myoepi_pheno,myoepi_pheno %in% MYOEP9_clust,"MYOEP9")
# myoepi_pheno<-replace(epi_pheno,epi_pheno %in% MYOEP2_clust,"MYOEP2")
data_fs_clusters_myoep$phenotype<-myoepi_pheno
# # data_immune[mast_rows,]$stroma_pheno<-"mast"
# # data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_myoep, file="200401_DCIScohort_FLOWSOM_Myoep16_phenotype.csv",row.names = FALSE)

#############################################  ASSESS EPI CLUSTER FREQ BY CONDITIONS  #############################################  

# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

##..Import data..##

datamyoep<-read.csv("200401_DCIScohort_FLOWSOM_Myoep16_phenotype.csv")

##..Subset the cell data and sample identifiers..##

cell_data_myoep<-datamyoep %>% select(Point_Num, phenotype)

##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_myoep <- as.data.frame(table(cell_data_myoep$Point_Num, cell_data_myoep$phenotype))
names(Freqs_myoep) <- c("SampleID","phenotype","count")

##..Get overall totals on a per sample basis..##

cell_totals_myoep<-aggregate(Freqs_myoep$count, by=list(Category=Freqs_myoep$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_myoep$SampleID)) {
  frequencies <- Freqs_myoep[Freqs_myoep$SampleID==i,"count"] / cell_totals_myoep[cell_totals_myoep$Category==i,2]
  Freqs_myoep[Freqs_myoep$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_myoep<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_myoep$Point_Num)
# filter annotation data by PointNum
annotation_data_myoep_cohort <- droplevels(annotation_data_myoep[annotation_data_myoep$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_myoep_cohort$PointNumber <- factor(annotation_data_myoep_cohort$PointNumber, levels=pointNums)
annotation_data_cohort_myoep<-annotation_data_myoep_cohort[order(annotation_data_myoep_cohort$PointNumber),]
# cast the annotations to the frequency data
anno_data_myoep<-rep(annotation_data_myoep_cohort,1)
annotated_Freqs_myoep<-cbind(Freqs_myoep, anno_data_myoep)



# cell totals per sample
write.csv(cell_totals_myoep, file="Myoep_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_myoep, file="Myoep_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_myoep, file="Myoep_freqs_per_sample_annotated.csv",row.names = FALSE)
# write out single cell data epi
write.csv(data_fs_clusters_myoep, file="SingleCellData_Myoep_annotated.csv",row.names = FALSE)
# write out single cell data stroma
# write.csv(data_fs_clusters_stroma, file="SingleCellData_T_annotated.csv",row.names = FALSE)

################################ PLOT MYOEP Phenos IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(6)


data_myoep<-read.csv("Myoep_freqs_per_sample_annotated.csv")
plot_data<-data_myoep
#plot_data<-droplevels(data_epi[data_epi$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_myoep[data_myoep$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
#plot_data<-droplevels(data_epi[data_epi$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$phenotype),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$phenotype <- factor(plot_data$phenotype, levels=cluster_order)
plot_data<-plot_data[order(plot_data$phenotype),]

##..Plot all clusters across all samples..##
pdf(file = "200401_myoepFreqsAllbreast", height =8, width = 10) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(phenotype),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(phenotype))) + 
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

pdf(file = "200401_myoepAcrossallbreast", height =12, width = 40) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(phenotype))) + 
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


# status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
#   geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(phenotype))) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~phenotype, scales = "free")
# status

pdf(file = "200401_myoepAcrossCaseCtrl", height =6, width = 7 ) 
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon','normal','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])


pdf(file = "200401_MyoepAcrossCaseCtrl_invasiveCasesOnly", height =6, width = 7) 
status<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('tonsil','lymphnode','placenta','colon','concurrent','contcis','ipscis','continv','ipsinv','normal'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

pdf(file = "200401_myoepAcrossOutcome", height =6, width = 7) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
Tissue
dev.off()



#IMMUNE CLUST##################################################################################################################################################################################################
#########################################################################################################################################################
###########################################################################################################################
######################..FlowSOM Cluster round 5: Immune Cells..############################################
#############################################################################################################################
##################################################################################################################################################################################################
##################################################################################################################################################################################################
library(dplyr)
# stroma<-c("stroma") #can append if wanting to include non-immune clusters that look contaminated
# data_stroma<-data_fs_clusters[data_fs_clusters$lineage %in% stroma,]

immune<-c("immune") #can append if wanting to include non-immune clusters that look contaminated
dataimmune<-data_fs_clusters[data_fs_clusters$lineage %in% immune,]
# immunectrls<-c("tonsil","lymphnode")
# data_immune_trim<-droplevels(data_immune[!data_immune$Tissue %in% immunectrls,])

ff_immune <- flowFrame(exprs = data.matrix(dataimmune[dataimmune$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"), -c(63:73)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_3=c("CD45","Tryptase","HLADRDPDQ","CD3","CD20","CD4","CD8", "CD11c", "CD68", "CD14")

# sc.data.scale.transform.norm[,4:55] <- data.frame(t(t(sc.data[,4:55]) / as.numeric(percentile.vector)))
# #sc.data.scale.transform.norm<-sc.data.scale.transform.norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)
# 
# #TrimBadPoints
# badpoints<-c(2101,2301,6101,5311,5501)
# sc.data.scale.transform.norm.trimbad<-droplevels(sc.data.scale.transform.norm[!sc.data.scale.transform.norm$Point_Num %in% badpoints,]
##..Run FlowSOM random seed for reproducibility..##

set.seed(656)
out_fSOM_immune <- FlowSOM::ReadInput(ff_immune, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_immune <- FlowSOM::BuildSOM(out_fSOM_immune, colsToUse = clusterChannels_3, xdim=5, ydim=5)
out_fSOM_immune <- FlowSOM::BuildMST(out_fSOM_immune)
labels_immune <- out_fSOM_immune$map$mapping[,1]

out_fSOM_immune <- UpdateNodeSize(out_fSOM_immune, reset = TRUE)
#tiff("190203_FS-round2-nodes.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_immune, view = "grid", markers = clusterChannels_3)
#dev.off()
out_fSOM_immune <- UpdateNodeSize(out_fSOM_immune)
#tiff("190203_FS-round2-MST.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_immune, view = "MST", markers = clusterChannels_3)
FlowSOM::PlotStars(out_fSOM_immune, view="tSNE",markers = clusterChannels_3)
#dev.off()

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_immune<-out_fSOM_immune[["map"]][["mapping"]]
fs_clusters_immune<-fs_clusters_immune[,1]
data_fs_clusters_immune <- cbind(dataimmune, fs_clusters_immune)

# go through all clusters and calculate mean for every channel
hm_allclusters_immune <- matrix(, nrow = length(unique(fs_clusters_immune)), ncol = length(clusterChannels_3))
for(i in 1:length(unique(fs_clusters_immune))) {
  temp_mat <- data_fs_clusters_immune[data_fs_clusters_immune[,"fs_clusters_immune"] == i, clusterChannels_3]
  hm_allclusters_immune[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_immune) <- paste("cluster", 1:length(unique(fs_clusters_immune)), sep = "")
colnames(hm_allclusters_immune) <- clusterChannels_3  
hm_allclusters_immune

# plot heatmap of all clusters
#tiff("plots/draft_figs/190606_FS-myeloid-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/190606_FS-myeloid-hmap.eps",width=15, height=11)
pdf(file = "200401_immuneFLOWSOM25", height = 4, width = 5)  
heatmap.2(hm_allclusters_immune, 
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
table(data_fs_clusters_immune$fs_clusters_immune)
write.csv(data_fs_clusters_immune, file="200401_DCIScohort_FLOWSOMimmune20.csv",row.names = FALSE)

cluster <-23
table(data_fs_clusters_immune[data_fs_clusters_immune$fs_clusters_immune==cluster,]$Point_Num)
hist(data_fs_clusters_stroma[data_fs_clusters_stroma$fs_clusters_stroma==cluster,]$VIM)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
MONO_clust<-c(7,6,1)
IMMUNEOTHER_clust<-c(23)
MACS_clust<-c(2,12,4,8)
DC_clust<-c(18)
DCMONO_clust<-c(3)
MACSDC_clust<-c(17,11)
MAST_clust<-c(16,22,21)
CD8T_clust<-c(25,24,20)
CD4T_clust<-c(10,15,14,5,9)
BCELL_clust<-c(13)
TCELL_clust<-c(19)
immune_pheno<-fs_clusters_immune
immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO_clust,"MONO")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACS_clust,"MACS")
immune_pheno<-replace(immune_pheno,immune_pheno %in% DC_clust,"DC")
immune_pheno<-replace(immune_pheno,immune_pheno %in% DCMONO_clust,"DC")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACSDC_clust,"MACSDC")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MAST_clust,"MAST")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD8T_clust,"CD8T")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T_clust,"CD4T")
immune_pheno<-replace(immune_pheno,immune_pheno %in% BCELL_clust,"BCELL")
immune_pheno<-replace(immune_pheno,immune_pheno %in% TCELL_clust,"TCELL")
immune_pheno<-replace(immune_pheno,immune_pheno %in% IMMUNEOTHER_clust,"IMMUNEOTHER")

data_fs_clusters_immune$sublineage<-immune_pheno
# data_immune[mast_rows,]$stroma_pheno<-"mast"
# data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_immune, file="200401_DCIScohort_FLOWSOM_immune25_subLineage.csv",row.names = FALSE)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
CD8T1_clust<-c(25)
CD8T2_clust<-c(24)
CD8T3_clust<-c(20)
CD4T1_clust<-c(10)
CD4T2_clust<-c(15)
CD4T3_clust<-c(14)
CD4T4_clust<-c(5)
CD4T5_clust<-c(9)
MACS1_clust<-c(2)
MACS2_clust<-c(12)
MACS3_clust<-c(4)
MACS4_clust<-c(8)
MAST1_clust<-c(16)
MAST2_clust<-c(22)
MAST3_clust<-c(21)
BCELL1_clust<-c(13)
TCELL1_clust<-c(19)
IMMUNEOTHER1_clust<-c(23)
MONO1_clust<-c(7)
MONO2_clust<-c(6)
MONO3_clust<-c(1)
MACSDC1_clust<-c(17)
MACSDC2_clust<-c(11)
DC1_clust<-c(18)
DCMONO1_clust<-c(3)
# TUMOR8_clust<-c(2)
# TUMOR9_clust<-c(6)
immune_pheno<-fs_clusters_immune
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD8T1_clust,"CD8T1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD8T2_clust,"CD8T2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD8T3_clust,"CD8T3")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T1_clust,"CD4T1")  
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T2_clust,"CD4T2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T3_clust,"CD4T3")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T4_clust,"CD4T4")
immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T5_clust,"CD4T5")
immune_pheno<-replace(immune_pheno,immune_pheno %in% TCELL1_clust,"TCELL")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4T4_clust,"CD4T5")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% CD3T_clust,"CD3T")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% CD4CD8T_clust,"CD4CD8T")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACSDC1_clust,"MACSDC1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACSDC2_clust,"MACSDC2")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% MACSDC3_clust,"MACSDC3")
immune_pheno<-replace(immune_pheno,immune_pheno %in% DC1_clust,"DC1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% DCMONO1_clust,"DCMONO1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACS1_clust,"MACS1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACS3_clust,"MACS3")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACS2_clust,"MACS2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MACS4_clust,"MACS4")
immune_pheno<-replace(immune_pheno,immune_pheno %in% IMMUNEOTHER1_clust,"IMMUNEOTHER1")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% IMMUNEOTHER2_clust,"IMMUNEOTHER2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO1_clust,"MONO1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO2_clust,"MONO2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO3_clust,"MONO3")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO4_clust,"MONO4")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO5_clust,"MONO5")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% MONO6_clust,"MONO6")
immune_pheno<-replace(immune_pheno,immune_pheno %in% BCELL1_clust,"BCELL1")
# immune_pheno<-replace(immune_pheno,immune_pheno %in% BCELL2_clust,"BCELL2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MAST1_clust,"MAST1")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MAST2_clust,"MAST2")
immune_pheno<-replace(immune_pheno,immune_pheno %in% MAST3_clust,"MAST3")
data_fs_clusters_immune$phenotype<-immune_pheno
# data_immune[mast_rows,]$stroma_pheno<-"mast"
# data_immune[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_immune, file="200401_DCIScohort_FLOWSOM_immune25_phenotype.csv",row.names = FALSE)



#############################################  ASSESS immune CLUSTER FREQ BY CONDITIONS  #############################################  

# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# all cell types (out of total) broken down by sample ID. Saves the frequency data as a csv.

##..Import data..##

dataimmune<-read.csv("200401_DCIScohort_FLOWSOM_immune25_phenotype.csv")

##..Subset the cell data and sample identifiers..##

cell_data_immune<-dataimmune %>% select(Point_Num, phenotype)

##..Create a dataframe with the counts of each cell cluster across point number..##

Freqs_immune <- as.data.frame(table(cell_data_immune$Point_Num, cell_data_immune$phenotype))
names(Freqs_immune) <- c("SampleID","phenotype","count")

##..Get overall totals on a per sample basis..##

cell_totals_immune<-aggregate(Freqs_immune$count, by=list(Category=Freqs_immune$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Freqs_immune$SampleID)) {
  frequencies <- Freqs_immune[Freqs_immune$SampleID==i,"count"] / cell_totals_immune[cell_totals_immune$Category==i,2]
  Freqs_immune[Freqs_immune$SampleID==i,"frequency"] <- frequencies
}

##..Add back annotations of the Session, Status, Recurrence, and Tissue..##

# read in data
# annotation_data_immune<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_immune<-info.csv
# get list of pointnNums is frequency data
pointNums<-unique(cell_data_immune$Point_Num)
# filter annotation data by PointNum
annotation_data_immune_cohort <- droplevels(annotation_data_immune[annotation_data_immune$PointNumber %in% pointNums, ])
# ensure order of points matches that of cell frequency data
annotation_data_immune_cohort$PointNumber <- factor(annotation_data_immune_cohort$PointNumber, levels=pointNums)
annotation_data_cohort_immune<-annotation_data_immune_cohort[order(annotation_data_immune_cohort$PointNumber),]
# cast the annotations to the frequency data
anno_data_immune<-rep(annotation_data_immune_cohort,1)
annotated_Freqs_immune<-cbind(Freqs_immune, anno_data_immune)



# cell totals per sample
write.csv(cell_totals_immune, file="immune_phenotype_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Freqs_immune, file="immune_phenotype_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_immune, file="immune_phenotype_freqs_per_sample_annotated.csv",row.names = FALSE)
# write out single cell data immune
#write.csv(data_fs_clusters_immune, file="SingleCellData_immunesublineage_annotated.csv",row.names = FALSE)
# write out single cell data stroma
# write.csv(data_fs_clusters_stroma, file="SingleCellData_T_annotated.csv",row.names = FALSE)

################################ PLOT immune Phenos IN FACET BAR GRAPHS ###############################


# Overview: Reads in cell frequency data. Plots the frequency of all cell types in bulk, 
# per patient, and also for a given variable (ie. sesson, recurrence, etc) builds a facet grid 
# of frequency of all cell types across the variable type.

library(ggplot2)
library(forcats)
library(dplyr)

##..Make a lil custom color palette..##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color<-gg_color_hue(25)


data_immune<-read.csv("immune_phenotype_freqs_per_sample_annotated.csv")
plot_data<-data_immune
#plot_data<-droplevels(data_immune[data_immune$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_immune[data_immune$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
#plot_data<-droplevels(data_immune[data_immune$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$phenotype),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$phenotype <- factor(plot_data$phenotype, levels=cluster_order)
plot_data<-plot_data[order(plot_data$phenotype),]

##..Plot all clusters across all samples..##
pdf(file = "200401_immuneSublinFreqsAllbreast", height =6, width = 10) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(phenotype),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(phenotype))) + 
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

pdf(file = "200401_immuneSublinAcrossallbreast", height =12, width = 40) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(phenotype))) + 
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
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = phenotype)) +
#   scale_fill_manual(values=rev(color)) +1
#   theme +
#   facet_wrap(.~phenotype, scales = "free")
# session

# status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
#   geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(phenotype))) +
#   scale_fill_manual(values=rev(color)) +
#   theme +
#   facet_wrap(.~phenotype, scales = "free")
# status

pdf(file = "200401_immuneSublinAcrossCaseCtrl", height =8, width = 9) 
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon','normal','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

plot_data_invasive<-droplevels(plot_data[plot_data$Reccurence %in% c('ipsinv','continv','na'),])


pdf(file = "200401_immuneSublinAcrossCaseCtrl_invasiveCasesOnly", height =8, width = 9) 
status<-ggplot(plot_data_invasive[!plot_data_invasive$Status %in% c('tonsil','lymphnode','placenta','colon','normal','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
status
dev.off()

pdf(file = "2003023_immuneSublinAcrossOutcome", height =8, width = 9) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(phenotype))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~phenotype, scales = "free")
Tissue
dev.off()



##################################################################################################################################################################################################
##########################################################################################################################
######################..FlowSOM Cluster round 6: Other Cells..############################################
#############################################################################################################################
##################################################################################################################################################################################################
library(dplyr)
# stroma<-c("stroma") #can append if wanting to include non-other clusters that look contaminated
# data_stroma<-data_fs_clusters[data_fs_clusters$lineage %in% stroma,]

other<-c("unknown") #can append if wanting to include non-other clusters that look contaminated
dataother<-data_fs_clusters[data_fs_clusters$lineage %in% other,]
# otherctrls<-c("tonsil","lymphnode")
# data_other_trim<-droplevels(data_other[!data_other$Tissue %in% otherctrls,])

ff_other <- flowFrame(exprs = data.matrix(dataother[dataother$Tissue %in% c("DCIS", "concurrent", "NewCIS", "NewInv", "normal"), -c(63:73)]), desc = list(FIL = 1)) #must remove all nonquantitative columns from the flowframe e.g. TissueType, Status etc
# ff_new <- flowFrame(exprs = data.matrix(sc.data.scale.transform.norm.trimbad[sc.data.scale.transform.norm.trimbad$Tissue %in% c('DCIS','concurrent','NewEvent','normal','tonsil','lymphnode'),-c(1)]), desc = list(FIL = 1)) #exclude tissue type column
clusterChannels_6=c("CK7", "CK5","VIM","P63","PanKRT","ECAD","SMA","CD45","Tryptase","HLADRDPDQ","CD3","CD20","CD4","CD8", "CD11c", "CD68", "CD14","FOXP3","MPO")

# sc.data.scale.transform.norm[,4:55] <- data.frame(t(t(sc.data[,4:55]) / as.numeric(percentile.vector)))
# #sc.data.scale.transform.norm<-sc.data.scale.transform.norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)
# 
# #TrimBadPoints
# badpoints<-c(2101,2301,6101,5311,5501)
# sc.data.scale.transform.norm.trimbad<-droplevels(sc.data.scale.transform.norm[!sc.data.scale.transform.norm$Point_Num %in% badpoints,]
##..Run FlowSOM random seed for reproducibility..##

set.seed(655)
out_fSOM_other <- FlowSOM::ReadInput(ff_other, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_other <- FlowSOM::BuildSOM(out_fSOM_other, colsToUse = clusterChannels_6, xdim=4, ydim=4)
out_fSOM_other <- FlowSOM::BuildMST(out_fSOM_other)
labels_other <- out_fSOM_other$map$mapping[,1]

out_fSOM_other <- UpdateNodeSize(out_fSOM_other, reset = TRUE)
#tiff("190203_FS-round2-nodes.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_other, view = "grid", markers = clusterChannels_6)
#dev.off()
out_fSOM_other <- UpdateNodeSize(out_fSOM_other)
#tiff("190203_FS-round2-MST.tiff", units="in", width=15, height=11, res=300)
FlowSOM::PlotStars(out_fSOM_other, view = "MST", markers = clusterChannels_6)
FlowSOM::PlotStars(out_fSOM_other, view="tSNE",markers = clusterChannels_6)
#dev.off()

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_other<-out_fSOM_other[["map"]][["mapping"]]
fs_clusters_other<-fs_clusters_other[,1]
data_fs_clusters_other <- cbind(dataother, fs_clusters_other)

# go through all clusters and calculate mean for every channel
hm_allclusters_other <- matrix(, nrow = length(unique(fs_clusters_other)), ncol = length(clusterChannels_6))
for(i in 1:length(unique(fs_clusters_other))) {
  temp_mat <- data_fs_clusters_other[data_fs_clusters_other[,"fs_clusters_other"] == i, clusterChannels_6]
  hm_allclusters_other[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_other) <- paste("cluster", 1:length(unique(fs_clusters_other)), sep = "")
colnames(hm_allclusters_other) <- clusterChannels_6  
hm_allclusters_other

# plot heatmap of all clusters
#tiff("plots/draft_figs/190606_FS-myeloid-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/190606_FS-myeloid-hmap.eps",width=15, height=11)
pdf(file = "200401_otherFLOWSOM25", height = 3, width = 4)  
heatmap.2(hm_allclusters_other, 
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
table(data_fs_clusters_other$fs_clusters_other)
write.csv(data_fs_clusters_other, file="200313_DCIScohort_FLOWSOMother16.csv",row.names = FALSE)

cluster <-9
table(data_fs_clusters_other[data_fs_clusters_other$fs_clusters_other==cluster,]$Point_Num)
hist(data_fs_clusters_stroma[data_fs_clusters_stroma$fs_clusters_stroma==cluster,]$VIM)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
MONO_clust<-c(5)
OTHER_clust<-c(16,15,7,10,13,9,11,4,3,8,12,6)
MACS_clust<-c(2)
MAST_clust<-c(14)
CD4T_clust<-c(1)
other_pheno<-fs_clusters_other
other_pheno<-replace(other_pheno,other_pheno %in% MONO_clust,"MONO")
other_pheno<-replace(other_pheno,other_pheno %in% MACS_clust,"MACS")
other_pheno<-replace(other_pheno,other_pheno %in% DC_clust,"DC")
other_pheno<-replace(other_pheno,other_pheno %in% MAST_clust,"MAST")
other_pheno<-replace(other_pheno,other_pheno %in% CD4T_clust,"CD4T")
other_pheno<-replace(other_pheno,other_pheno %in% OTHER_clust,"OTHER")

data_fs_clusters_other$sublineage<-other_pheno
# data_other[mast_rows,]$stroma_pheno<-"mast"
# data_other[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_other, file="200401_DCIScohort_FLOWSOM_other16_subLineage.csv",row.names = FALSE)

#### NAME EACH CLUSTER ################################
##..Add cell type to each event..##
CD4T6_clust<-c(1)
MONO4_clust<-c(5)
OTHER1_nothing_clust<-c(7)
OTHER2_clust<-c(16,15,10,13,9,11,4,3,8,12,6)
MAST4_clust<-c(14)
MACS5_clust<-c(2)
# TUMOR8_clust<-c(2)
# TUMOR9_clust<-c(6)
other_pheno<-fs_clusters_other
other_pheno<-replace(other_pheno,other_pheno %in% CD4T6_clust,"CD4T6")
other_pheno<-replace(other_pheno,other_pheno %in% MONO4_clust,"MONO4")
other_pheno<-replace(other_pheno,other_pheno %in% OTHER1_nothing_clust,"OTHER1")
other_pheno<-replace(other_pheno,other_pheno %in% OTHER2_clust,"OTHER2")
other_pheno<-replace(other_pheno,other_pheno %in% MAST4_clust,"MAST4")
other_pheno<-replace(other_pheno,other_pheno %in% MACS5_clust,"MACS5")

table(data_fs_clusters_other$sublineage)
data_fs_clusters_other$phenotype<-other_pheno
# data_other[mast_rows,]$stroma_pheno<-"mast"
# data_other[neut_rows,]$stroma_pheno<-"neutrophil"
write.csv(data_fs_clusters_other, file="200401_DCIScohort_FLOWSOM_other25_phenotype.csv",row.names = FALSE)

##################################################################################
#####..Concatenate into one matrix with terminal cell-type ID for all cells..#####
##################################################################################
# immunectrls<-c("tonsil","lymphnode")
# data_meta_clusters30_cohort<-droplevels(data_meta_clusters30[!data_meta_clusters30$Tissue %in% immunectrls,])
# 
# 
# data_meta_clusters30_cohort$cell_type<-data_meta_clusters30_cohort$lineage
# data_meta_clusters30_cohort[rownames(data_fs_clusters_immune),]$cell_type<-data_fs_clusters_immune$fs_clusters_immune
# 
# data_meta_clusters30_cohort[rownames(data_fs_clusters_stroma),]$cell_type<-data_fs_clusters_stroma$fs_clusters_stroma
# data_gran$cell_type<-data_gran_norm$cell_type
# 
# #manually annotate gdT cells
# CD3_thresh<-mean(data_gran_norm[data_gran_norm$cell_type=="CD4_T",]$CD3)
# gdT_row_idx<-which(data_gran_norm$gdTCR > 0.5 & data_gran_norm$CD3 >= CD3_thresh)
# gdT_rows<-rownames(data_gran_norm[gdT_row_idx,])
# data_gran[gdT_rows,]$cell_type<-"gdT_cell"
# data_gran_norm[gdT_rows,]$cell_type<-"gdT_cell"
# 
# #update cell lineage in case any gd T cells had been originally assigned incorrectly as immunethlelial, endothelial, fibroblast, or unidentified
# 
# incorrect_rows<-rownames(data_gran_norm[data_gran_norm$cell_type=="gdT_cell" & data_gran_norm$lineage!="immune",])
# data_gran_norm[incorrect_rows,]$lineage="immune"
# data_gran[incorrect_rows,]$lineage="immune"
# 
# #add details for major cell type (lymphocyte, myeloid, granulocyte, fibroblast, endothelial, epithelial)
# myeloid<-c(unique(data_myeloid$myeloid_pheno))
# lymphocyte<-c("CD8_T","CD4_T","B_cell","Treg","gdT_cell")
# granulocyte<-c("mast","neutrophil")
# imm_other<-c("imm_other")
# nonimmune<-c("fibroblast","endothelial","epithelial","unidentified")
# 
# 
# data_gran_norm<-data_gran_norm %>% mutate(cell_lin=case_when(data_gran_norm$cell_type %in% myeloid ~ "myeloid",
#                                                              data_gran_norm$cell_type %in% lymphocyte ~ "lymphocyte",
#                                                              data_gran_norm$cell_type %in% granulocyte ~ "granulocyte",
#                                                              data_gran_norm$cell_type %in% imm_other ~ "other",
#                                                              data_gran_norm$cell_type %in% nonimmune ~ "nonimmune")) 
# data_gran<-data_gran %>% mutate(cell_lin=case_when(data_gran$cell_type %in% myeloid ~ "myeloid",
#                                                    data_gran$cell_type %in% lymphocyte ~ "lymphocyte",
#                                                    data_gran$cell_type %in% granulocyte ~ "granulocyte",
#                                                    data_gran$cell_type %in% imm_other ~ "other",
#                                                    data_gran$cell_type %in% nonimmune ~ "nonimmune")) 
# 
# write.csv(data_gran_norm, file="granA_cellpheno_CS-asinh-norm.csv",row.names = FALSE) #save the annotated asinh-99th percentile scaled data
# write.csv(data_gran, file="granA_cellpheno_CS.csv",row.names = FALSE) #save the annotatetd CS, but untransformed data

####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################################################################################################################
######################## CONCATENATE ################################################################################################################
####################################################################################################################################################################################################################################
##################################################################################################################################################
########################################################################################################################################################

##### IMMMMPOORRRTTANNNNTTTTT ##########
#YOU MUST GO INTO EACH OF THESE CSV FILES AND MAKE THE "fs_clusters_XXXX" read "fs_cluster_lineage
myoep<-read.csv("200401_DCIScohort_FLOWSOM_Myoep6_phenotype.csv")
endo<-read.csv("200401_DCIScohort_FLOWSOM_Endo4_phenotype.csv")
tumor<-read.csv("200401_DCIScohort_FLOWSOM_tumor6_phenotype.csv")
immune<-read.csv("200401_DCIScohort_FLOWSOM_immune25_phenotype.csv")
stroma<-read.csv("200401_DCIScohort_FLOWSOM_Fibro6_phenotype.csv")
other<-read.csv("200401_DCIScohort_FLOWSOM_other25_phenotype.csv")
concat_frame<-rbind(tumor,immune,stroma,myoep,other,endo)
concat_frame<-concat_frame[order(concat_frame$Point_Num),]

write.csv(concat_frame, file="200325_CellTable_TumorImmuneEndoFibroMyoepOther.csv",row.names = FALSE)



# pdf(file = "200315_MPOhist", height =4, width = 4) 
hist(concat_frame$MPO)
# dev.off()



######## SET A THRESHOLD ##########
concatenatedcohort<-read.csv("200325_CellTable_TumorImmuneEndoFibroMyoepOther.csv")

# 
# ER_row_idx<-which(data_tumor$ER > 0.2) #use positivity threshold of 0.5
# ER_rows<-rownames(data_tumor[ER_row_idx,])
# data_tumor[ER_rows,]$lineage<-"12"

##### IMMMMPOORRRTTANNNNTTTTT ##########
#### GENERATE A NEW COLUMN IN YOUR CONCATEDNATED CSV FOR MPOPANKRT WHERE YOU TAKE MPO - PANKRT #########

MPO_row_idx<-which(concatenatedcohort$MPOPANKRT > 0.3) #use positivity threshold of 0.3
MPO_rows<-rownames(concatenatedcohort[MPO_row_idx,])
levels(concatenatedcohort$sublineage)<-c(levels(concatenatedcohort$sublineage),"NEUT" ) 
concatenatedcohort[MPO_rows,]$sublineage<-"NEUT"
write.csv(concatenatedcohort, file="200326_CellTable_TumorImmuneEndoFibroMyoepOtherNeut.csv",row.names = FALSE)

concatenatedcohort2<-read.csv("200326_CellTable_TumorImmuneEndoFibroMyoepOtherNeut.csv")
MPO_row_idx<-which(concatenatedcohort2$MPOPANKRT > 0.3) #use positivity threshold of 0.5
MPO_rows<-rownames(concatenatedcohort2[MPO_row_idx,])
levels(concatenatedcohort2$phenotype)<-c(levels(concatenatedcohort2$phenotype), "NEUT1" ) 
concatenatedcohort2[MPO_rows,]$phenotype<-"NEUT1"
write.csv(concatenatedcohort2, file="200326_CellTable_TumorImmuneEndoFibroMyoepOtherNeut.csv",row.names = FALSE)

######################################################################
######## Add_Cluster_Columns ########################################
######################################################################

##### IMMMMPOORRRTTANNNNTTTTT ##########

# I created a csv where I added extra columns for cluster annotation for future analyses - take the celltable above and add "cluster_code" and "phenotype_code"
# In future refinement or renaming of clusters, just change the code here to include them e.g. if Fibro1 and Fibro2 will be a future "CAF" just make CAF cluster coded for both, keep reference elsewhere
concatenated_cohort_cluster_number<-read.csv("200326_CellTable_TumorImmuneEndoFibroMyoepOtherNeut.csv")
concatenated_cohort_cluster_number<-concatenated_cohort_cluster_number %>% mutate(cluster_code=case_when(concatenated_cohort_cluster_number$sublineage %in% c("TUMOR") ~ 1,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("MYOEP") ~ 2,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("FIBROBLAST") ~ 3,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("ENDO") ~ 4,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("DCMONO") ~ 5,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("MONO") ~ 6,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("MACS") ~ 7,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("DC") ~ 8,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("MACSDC") ~ 9,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("MAST") ~ 10,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("CD8T") ~ 11,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("CD4T") ~ 12,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("BCELL") ~ 13,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("TCELL") ~ 14,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("IMMUNEOTHER") ~ 15,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("OTHER") ~ 16,
                                                                                  concatenated_cohort_cluster_number$sublineage %in% c("NEUT") ~ 17))
                                                                                  
concatenated_cohort_cluster_number<-concatenated_cohort_cluster_number %>% mutate(phenotype_code=case_when(concatenated_cohort_cluster_number$phenotype %in% c("TUMOR1") ~ 1,
                                                   concatenated_cohort_cluster_number$phenotype %in% c("TUMOR2") ~ 2,
                                                   concatenated_cohort_cluster_number$phenotype %in% c("TUMOR3") ~ 3,
                                                   concatenated_cohort_cluster_number$phenotype %in% c("TUMOR4") ~ 4,
                                                  concatenated_cohort_cluster_number$phenotype %in% c("TUMOR5") ~ 5,
                                                  concatenated_cohort_cluster_number$phenotype %in% c("TUMOR6") ~ 6,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MYOEP1") ~ 7,
                                                 concatenated_cohort_cluster_number$phenotype %in% c("MYOEP2") ~ 8,
                                                    concatenated_cohort_cluster_number$phenotype %in% c("MYOEP3") ~ 9,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MYOEP4") ~ 10,
                                                   concatenated_cohort_cluster_number$phenotype %in% c("MYOEP5") ~ 11,
                                                  concatenated_cohort_cluster_number$phenotype %in% c("MYOEP6") ~ 12,
                                                 concatenated_cohort_cluster_number$phenotype %in% c("ENDO1") ~ 13,
                                                 concatenated_cohort_cluster_number$phenotype %in% c("ENDO2") ~ 14,
                                               concatenated_cohort_cluster_number$phenotype %in% c("ENDO3") ~ 15,
                                              concatenated_cohort_cluster_number$phenotype %in% c("ENDO4") ~ 16,
                                             concatenated_cohort_cluster_number$phenotype %in% c("FIBRO1") ~ 17,
                                               # concatenated_cohort_cluster_number$phenotype %in% c("FIBRO2") ~ 18,
                                               concatenated_cohort_cluster_number$phenotype %in% c("CAF1") ~ 19,
                                               concatenated_cohort_cluster_number$phenotype %in% c("CAF2") ~ 20,
                                           concatenated_cohort_cluster_number$phenotype %in% c("CAF3") ~ 21,
                                            concatenated_cohort_cluster_number$phenotype %in% c("NORMFIBRO") ~ 22,
                                         # concatenated_cohort_cluster_number$phenotype %in% c("ENDO1") ~ 22,
                                         #     concatenated_cohort_cluster_number$phenotype %in% c("ENDO2") ~ 23,
                                              concatenated_cohort_cluster_number$phenotype %in% c("CD8T1") ~ 23,
                                               concatenated_cohort_cluster_number$phenotype %in% c("CD8T2") ~ 24,
                                          concatenated_cohort_cluster_number$phenotype %in% c("CD8T3") ~ 25,                                                     
                                          concatenated_cohort_cluster_number$phenotype %in% c("CD4T1") ~ 26,
                                            concatenated_cohort_cluster_number$phenotype %in% c("CD4T2") ~ 27,
                                      concatenated_cohort_cluster_number$phenotype %in% c("CD4T3") ~ 28,
                                           concatenated_cohort_cluster_number$phenotype %in% c("CD4T4") ~ 29,
                                              concatenated_cohort_cluster_number$phenotype %in% c("CD4T5") ~ 30,
                                              concatenated_cohort_cluster_number$phenotype %in% c("CD4T6") ~ 31,                    
                                              concatenated_cohort_cluster_number$phenotype %in% c("MACSDC1") ~ 32,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MACSDC2") ~ 33,
                                              concatenated_cohort_cluster_number$phenotype %in% c("DCMONO1") ~ 34,
                                                concatenated_cohort_cluster_number$phenotype %in% c("DC1") ~ 35,
                                                concatenated_cohort_cluster_number$phenotype %in% c("MACS1") ~ 36,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MACS2") ~ 37,
                                                concatenated_cohort_cluster_number$phenotype %in% c("MACS3") ~ 38,
                                             concatenated_cohort_cluster_number$phenotype %in% c("MACS4") ~ 39,
                                           concatenated_cohort_cluster_number$phenotype %in% c("MACS5") ~ 40,
                                             concatenated_cohort_cluster_number$phenotype %in% c("IMMUNEOTHER1") ~ 41,
                                              concatenated_cohort_cluster_number$phenotype %in% c("TCELL") ~ 42,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MONO1") ~ 43,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MONO2") ~ 44,
                                              concatenated_cohort_cluster_number$phenotype %in% c("MONO3") ~ 45,
                                           concatenated_cohort_cluster_number$phenotype %in% c("MONO4") ~ 46,
                                      concatenated_cohort_cluster_number$phenotype %in% c("MONO5") ~ 47,
                                      concatenated_cohort_cluster_number$phenotype %in% c("MONO6") ~ 48,
                                      concatenated_cohort_cluster_number$phenotype %in% c("MONO7") ~ 49,
                                               concatenated_cohort_cluster_number$phenotype %in% c("BCELL1") ~ 50,
                                                concatenated_cohort_cluster_number$phenotype %in% c("MAST1") ~ 51,
                                      concatenated_cohort_cluster_number$phenotype %in% c("MAST2") ~ 52,
                                      concatenated_cohort_cluster_number$phenotype %in% c("MAST3") ~ 53,
                                      concatenated_cohort_cluster_number$phenotype %in% c("MAST4") ~ 54,
                                      concatenated_cohort_cluster_number$phenotype %in% c("OTHER1") ~ 55,
                                      concatenated_cohort_cluster_number$phenotype %in% c("OTHER2") ~ 56,
                                               concatenated_cohort_cluster_number$phenotype %in% c("NEUT1") ~ 57),)

write.csv(concatenated_cohort_cluster_number, file="200326_CellTable_TumorImmuneStromaMyoepOtherNeut_ClusterCoded.csv",row.names = FALSE)



table(concatenated_cohort_cluster_number$sublineage)
table(concatenated_cohort_cluster_number$phenotype)



################################################### ################################################### ################################################### ################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ###################################################
################################################### GET FREQEUNCIES FROM WHOLE CONCATENATED FILE  ################################################### ####################################
############### ################################################### ################################################### ###################################################

ALLdata<-read.csv("200326_CellTable_TumorImmuneStromaMyoepOtherNeut_ClusterCoded.csv")



cell_data_ALL<-ALLdata %>% select(Point_Num, sublineage, phenotype, cluster_code, phenotype_code)

##..Create a dataframe with the counts of each cell cluster across point number..##

Sublineage_Freqs_ALL <- as.data.frame(table(cell_data_ALL$Point_Num, cell_data_ALL$sublineage))
names(Sublineage_Freqs_ALL) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

cell_totals_ALL<-aggregate(Sublineage_Freqs_ALL$count, by=list(Category=Sublineage_Freqs_ALL$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(Sublineage_Freqs_ALL$SampleID)) {
  frequencies <- Sublineage_Freqs_ALL[Sublineage_Freqs_ALL$SampleID==i,"count"] / cell_totals_ALL[cell_totals_ALL$Category==i,2]
  Sublineage_Freqs_ALL[Sublineage_Freqs_ALL$SampleID==i,"frequency"] <- frequencies
}
#..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_ALL<-info.csv
# get list of pointnNums is frequency data
pointNumsALL<-unique(cell_data_ALL$Point_Num)
# filter annotation data by PointNum
annotation_data_ALL <- droplevels(annotation_data_ALL[annotation_data_ALL$PointNumber %in% pointNumsALL, ])
# ensure order of points matches that of cell frequency data
annotation_data_ALL$PointNumber <- factor(annotation_data_ALL$PointNumber, levels=pointNumsALL)
annotation_data_ALL<-annotation_data_ALL[order(annotation_data_ALL$PointNumber),]
# cast the annotations to the frequency data
anno_data_ALL<-rep(annotation_data_ALL,1)
annotated_Freqs_ALL<-cbind(Sublineage_Freqs_ALL, anno_data_ALL)

# cell totals per sample
write.csv(cell_totals_ALL, file="Sublineage_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(Sublineage_Freqs_ALL, file="Sublineage_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_ALL, file="Sublineage_freqs_per_sample_annotated.csv",row.names = FALSE)


###PLOT####
color<-gg_color_hue(16)


data_ALL<-read.csv("Sublineage_freqs_per_sample_annotated.csv")
plot_data<-data_ALL
#plot_data<-droplevels(data_ALL[data_ALL$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_ALL[data_ALL$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
#plot_data<-droplevels(data_ALL[data_ALL$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$cell_type),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$cell_type <- factor(plot_data$cell_type, levels=cluster_order)
plot_data<-plot_data[order(plot_data$cell_type),]

##..Plot all clusters across all samples..##
pdf(file = "200325_SUBLINEAGE_FreqsAllbreast", height =6, width = 10) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(cell_type),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(cell_type))) + 
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

pdf(file = "200325_SUBLINEAGE_Acrossallbreast", height =12, width = 40) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(cell_type))) + 
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
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = cell_type)) +
#   scale_fill_manual(values=rev(color)) +1
#   theme +
#   facet_wrap(.~cell_type, scales = "free")
# session

status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status

pdf(file = "200325_PHENOTYPE_AcrossCaseCtrl", height =10, width = 15) 
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon','normal','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status
dev.off()

pdf(file = "200325_PHENOTYPE_AcrossOutcome", height =8, width = 9) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
Tissue
dev.off()



######PHENOTYPE######

phenotype_Freqs_ALL <- as.data.frame(table(cell_data_ALL$Point_Num, cell_data_ALL$phenotype))
names(phenotype_Freqs_ALL) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

cell_totals_ALL<-aggregate(phenotype_Freqs_ALL$count, by=list(Category=phenotype_Freqs_ALL$SampleID), FUN=sum)

##..Determine frequecy of each cell type in each sample..##

for(i in unique(phenotype_Freqs_ALL$SampleID)) {
  frequencies <- phenotype_Freqs_ALL[phenotype_Freqs_ALL$SampleID==i,"count"] / cell_totals_ALL[cell_totals_ALL$Category==i,2]
  phenotype_Freqs_ALL[phenotype_Freqs_ALL$SampleID==i,"frequency"] <- frequencies
}
#..Add back annotations of the Session, Status, Recurrence, and Tissue..##
# read in data
# annotation_data_epi<-read.csv("~/Desktop/DCIS/Segmentation/NewWatershed/Metadata/200120_COHORT_METADATA.csv")
annotation_data_ALL<-info.csv
# get list of pointnNums is frequency data
pointNumsALL<-unique(cell_data_ALL$Point_Num)
# filter annotation data by PointNum
annotation_data_ALL <- droplevels(annotation_data_ALL[annotation_data_ALL$PointNumber %in% pointNumsALL, ])
# ensure order of points matches that of cell frequency data
annotation_data_ALL$PointNumber <- factor(annotation_data_ALL$PointNumber, levels=pointNumsALL)
annotation_data_ALL<-annotation_data_ALL[order(annotation_data_ALL$PointNumber),]
# cast the annotations to the frequency data
anno_data_ALL<-rep(annotation_data_ALL,1)
annotated_Freqs_ALL<-cbind(phenotype_Freqs_ALL, anno_data_ALL)

# cell totals per sample
write.csv(cell_totals_ALL, file="phenotype_totals.csv",row.names = FALSE)
# cluster frequencies per sample
write.csv(phenotype_Freqs_ALL, file="phenotype_freqs_per_sample.csv",row.names = FALSE)
# cluster frequencies per sample with annotations
write.csv(annotated_Freqs_ALL, file="phenotype_freqs_per_sample_annotated.csv",row.names = FALSE)

###PLOT####
color<-gg_color_hue(54)


data_ALL<-read.csv("phenotype_freqs_per_sample_annotated.csv")
plot_data<-data_ALL
#plot_data<-droplevels(data_ALL[data_ALL$Tissue %in% c('tonsil','lymphnode','DCIS','normal','NewEvent','concurrent'),])
plot_data<-droplevels(data_ALL[data_ALL$Tissue %in% c('DCIS','normal','NewCIS','NewInv','concurrent'),])
#plot_data<-droplevels(data_ALL[data_ALL$Tissue %in% c('DCIS'),])

# ## If you want to only include a portion of the patients then use the following:
# condition<-c('ctrl','case') 
# # ex. condition<-('tonsil','DCIS')
# # ex. condition<-c(1101, 1102)
# plot_data<-droplevels(data[data$Status %in% condition, ])

# reorder by descending median frequency (clusters)
cluster_order<-levels(fct_reorder(as.factor(plot_data$cell_type),plot_data$frequency,.fun=median,.desc=FALSE))
plot_data$cell_type <- factor(plot_data$cell_type, levels=cluster_order)
plot_data<-plot_data[order(plot_data$cell_type),]

##..Plot all clusters across all samples..##
pdf(file = "200325_PHENOTYPE_FreqsAllbreast", height =6, width = 15) 
cell_box<-ggplot(plot_data, aes(x=fct_reorder(as.factor(cell_type),frequency,.fun=median,.desc=TRUE), y=frequency, fill=as.factor(cell_type))) + 
  geom_boxplot() +
  scale_fill_manual(values = rev(color)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 20)) +
  labs(x="Cluster") +
  labs(y="Frequency of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
cell_box
dev.off()

##..Plot frequency broken down by Point Number..##

pdf(file = "200325_PHENOTYPE_Acrossallbreast", height =12, width = 40) 
cell_bar<-ggplot(plot_data, aes(x=as.factor(SampleID), y=frequency, fill=as.factor(cell_type))) + 
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
#   geom_boxplot(aes(x=as.factor(session), y=frequency, fill = cell_type)) +
#   scale_fill_manual(values=rev(color)) +1
#   theme +
#   facet_wrap(.~cell_type, scales = "free")
# session

status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('normal', 'ctrl', 'case','concurrent','contcis','ipscis','continv','ipsinv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status

pdf(file = "200325_PHENOTYPE_AcrossCaseCtrl", height =10, width = 15) 
status<-ggplot(plot_data[!plot_data$Status %in% c('tonsil','lymphnode','placenta','colon','normal','concurrent','contcis','ipscis','continv','ipsinv'),]) +
  geom_boxplot(aes(x=factor(Status, level = c('ctrl', 'case')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
status
dev.off()

pdf(file = "200325_PHENOTYPE_AcrossOutcome", height =12, width = 15) 
Tissue<-ggplot(plot_data[!plot_data$Tissue %in% c('tonsil','lymphnode','placenta','colon'),]) +
  geom_boxplot(aes(x=factor(Tissue, level = c('normal', 'DCIS','NewCIS','concurrent','NewInv')), y=frequency, fill = as.factor(cell_type))) +
  scale_fill_manual(values=rev(color)) +
  theme +
  facet_wrap(.~cell_type, scales = "free")
Tissue
dev.off()



###Combine Mask Data####################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### ################################################### ################################################### ###################################################
################################################### ################################################### #

data_cohort<-read.csv("200406_CellTable_TumorImmuneStromaMyoepOtherNeut_ClusterCoded.csv")

data_mask<-read.csv("200407_cell_mask_annotations.csv")

data_cohort_masklocations<-cbind(data_cohort, data_mask)










###SUBLINEAGEHEATMAP####

HMChannels_7=c("CK7","CK5","VIM","P63","PanKRT","ECAD","SMA","CD45","Tryptase","HLADRDPDQ","CD3","CD20","CD4","CD8", "CD11c", "CD68", "CD14","FOXP3","MPO")
# go through all clusters and calculate mean for every channel
###SUBLINEAGEHEATMAP####
hm_metaclusters_ALL <- matrix(, nrow = chosen_k, ncol = length(HMChannels_7))
for(i in 1:chosen_k) {
  temp_mat <- concatenated_cohort_cluster_number[concatenated_cohort_cluster_number[,"sublineage"] == i, HMChannels_7]
  concatenated_cohort_cluster_number[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
rownames(concatenated_cohort_cluster_number) <- paste("cluster", 1:chosen_k, sep = "")
colnames(concatenated_cohort_cluster_number) <- clusterChannels_7
concatenated_cohort_cluster_number

# make a metacluster heatmap

#tiff("plots/190514/190514_FS_MC-lymph-hmap.tiff", units="in", width=15, height=11, res=300)
heatmap.2(hm_metaclusters_lymph,
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.6, cexCol = 0.9, margins = c(8,14))
#dev.off()




