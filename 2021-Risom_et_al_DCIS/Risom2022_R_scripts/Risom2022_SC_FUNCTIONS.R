############################### Author: Sam Kimmey ###########################################
#
# Below are general use functions to source into my RMD files - specifically for IONpath WP
#

################################# Libraries ##################################################

#library(flowCore) # for FCS file management and data access
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

################################# Functions ##################################################

med_sd_violin <- function(x) {
  # This function is used to display median +/- SD for violin plots
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

asinTransform <- function(dt, fa=factors, arc.sin.h.val = 5) {
  # asinh transforms
  # Inputs:
  #   dat - data.table
  #   fa - character vector of channels not to transform
  # Outputs:
  #   dt - data.table
  
  if(!is.data.table(dt)){
    dt <- data.table(dt)
  }
  print("Asinh transforming")
  to.transform <- colnames(dt)[!colnames(dt) %in% fa]
  dt[, (to.transform) := asinh(dt[, to.transform, with=F]/arc.sin.h.val)] # replace each column with the ashin h transformed values
  return(dt)
}

printViolinsForParam <- function(dt=csvdt, plot.vars, t=title, y_lab = "Normalized Expression"){
  # this function is to generate and save a violin plot of channels
  # Inputs:
  # dt - data table with single cell expression values
  # fa - vector of str factor channels
  # f - filename to save plot as (.png)
  
  dt = sc.data[sc.data$gate.name == "B.cell",]
  plot.vars <- c("pS6") 
  y_lab = "Norm Exp"
  vio_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  dt.long <- melt(dt, measure.vars = plot.vars)
  
  g <- ggplot(dt.long, aes(x = variable, y = value)) + labs(title = t) +
    geom_violin(scale = "width") + 
    stat_summary(fun.data=med_sd_violin, geom="pointrange") + 
    #geom_jitter(width = 0.3, size = 0.1) + 
    vio_theme + ylab(y_lab) + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "bold"))
  
  print(g)
}

printViolinsForParam.facet <- function(dt=csvdt, vars, t=title, xparam, y_lab = "Transformed expression\n asinh(x/5)"){
  # this function is to generate and save a violin plot of channels
  # Inputs:
  # dt - data table with single cell expression values
  # fa - vector of str factor channels
  # f - filename to save plot as (.png)
  vio_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  
  dt.long <- melt(dt, measure.vars = vars)
  
  g <- ggplot(dt.long, aes(x = dt.long[[xparam]], y = value)) + labs(title = t) +
    geom_violin(scale = "width") + 
    stat_summary(fun.data=med_sd_violin, geom="pointrange") + 
    #geom_jitter(width = 0.3, size = 0.1) + 
    vio_theme + ylab(y_lab) + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "bold"))
  
  print(g)
}

create.uMAP.csv <- function(x=csvTable, p=parameters, config = umap.defaults){
  # function to generate uMAP axis for parameters of interest
  # x - sub sampled and transformed data table
  # p - parameters to use for dimension reduction
  start_time <- Sys.time()
  
  set.seed(403)
  if(!is.data.table(x)){
    x <- data.table(x)
  }
  
  umapMat <- as.matrix(x[,..p])
  
  #print("Parameters for umap:")
  #print(p)
  
  u <- umap(umapMat, config)
  
  x[,umap1:=u$layout[,1]]
  x[,umap2:=u$layout[,2]]
  end_time <- Sys.time()
  print(paste0("Execution begin: ", start_time))
  print(paste0("Execution end: ", end_time))
  return(x)
}

print_umap_by_tumor <- function(umapDF) {
  

  wes <- wes_palette("Darjeeling1", length(unique(umapDF$Tissue_Type)), type = "continuous")
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = Tissue_Type) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 8, byrow = T)) +
    labs(title = "UMAP overlaid with Tissue Site")
  
  
}
print_umap_by_status <- function(umapDF) {
  
  
  wes <- wes_palette("Zissou1", length(unique(umapDF$Status)), type = "continuous")
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = Status) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 8, byrow = T)) +
    labs(title = "UMAP overlaid with Tissue Site")
}

print_umap_by_phenotype <- function(umapDF) {
  
  
  #wes <- wes_palette("Darjeeling1", length(unique(umapDF$phenotype)), type = "continuous")
  #gg_color_hue <- function(n) 
    #hues = seq(15, 375, length = n + 1)
    #hcl(h = hues, l = 65, c = 100)[1:n]
  
  
  wes<-gg_color_hue(30)

  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = phenotype) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 14, byrow = T)) +
    labs(title = "UMAP overlaid with Phenotype")
}

print_umap_by_sublineage <- function(umapDF) {
  
  
  #wes <- wes_palette("Darjeeling1", length(unique(umapDF$phenotype)), type = "continuous")
  #gg_color_hue <- function(n) 
  #hues = seq(15, 375, length = n + 1)
  #hcl(h = hues, l = 65, c = 100)[1:n]
  
  
  # wes<-gg_color_hue(16)
  wes<-c("#63B747","#E4388E","#660066","#CC99FF","#5C311D","#C7003C","#FF7878","#C1C1C1","#0F5428","#3FE0D0","#DFBF9F","#967E69","#054CA3","#FFBF00","#B8AFD5","#7AC4FF")
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = sublineage) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    # guides(color=guide_legend(nrow = 14, byrow = T)) +
    labs(title = "UMAP overlaid with Sublineage")
}

print_umap_by_lineage <- function(umapDF) {
  
  
  #wes <- wes_palette("Darjeeling1", length(unique(umapDF$phenotype)), type = "continuous")
  #gg_color_hue <- function(n) 
  #hues = seq(15, 375, length = n + 1)
  #hcl(h = hues, l = 65, c = 100)[1:n]
  
  
 wes<-gg_color_hue(5)
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = celllineage) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 14, byrow = T)) +
    labs(title = "UMAP overlaid with Lineage")
}

print_umap_by_compartment <- function(umapDF) {
  
  
  wes <- wes_palette("Zissou1", length(unique(umapDF$phenotype)), type = "continuous")
  #gg_color_hue <- function(n) 
  #hues = seq(15, 375, length = n + 1)
  #hcl(h = hues, l = 65, c = 100)[1:n]
  
  
   # wes<-gg_color_hue(2)
  
  # umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = compartment) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    # guides(color=guide_legend(nrow = 14, byrow = T)) +
    labs(title = "UMAP overlaid with Lineage")
}

print_umap_by_annotation <- function(umapDF) {
  
  # umapDF <- umapDF[phenotype!="none",]
  wes <- wes_palette("Zissou1", length(unique(umapDF$phenotype)), type = "continuous")
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = umapDF) +
    aes(x = umap1, y = umap2, color = phenotype) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_color_discrete() +
    #scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 12, ncol = 1, byrow = T)) +
    labs(title = "UMAP overlaid with cluster annotation", color = "Annotation\nGroups")
  
  
}



print_umap_by_annotation_facet <- function(umapDF) {
  
  umapDF <- umapDF[metaClusterAnnotation!="none",]
  wes <- wes_palette("Zissou1", length(unique(umapDF$metaClusterAnnotation)), type = "continuous")
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF.melt <- melt(umapDF, measure.vars = "metaClusterAnnotation")
  
  
  ggplot(data = umapDF.melt) + facet_wrap(~value) +
    aes(x = umap1, y = umap2, color = value) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_color_discrete() +
    #scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 12, ncol = 1, byrow = F)) +
    labs(title = "UMAP overlaid with cluster annotation", color = "Annotation\nGroups")
  
  
}


print_umap_by_tumor_facet <- function(umapDF) {
  
  umapDF <- ss.data
  wes <- wes_palette("Zissou1", length(unique(umapDF$Tissue)), type = "continuous")
  
  #umapDF <- as.data.table(sapply(umapDF, sample, nrow(umapDF)))
  #umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  umapDF.melt <- melt(umapDF, measure.vars = "Tissue")
  
  ggplot(data = umapDF.melt) + facet_wrap(~value) +
    aes(x = umap1, y = umap2, color = value) +
    geom_point(size = 0.5) +
    theme_minimal() + #geom_jitter() +
    scale_colour_manual(values =wes) +
    guides(color=guide_legend(nrow = 8, byrow = T)) +
    labs(title = "UMAP overlaid with Tissue Site") + theme(legend.position = "none")
  
  
}

printUmapByParamList <- function(dt=dataTable, param=parameter){
  # The dataTable must contain FlowSOM metaclusters
  # This function saves a png of the umap of choice, colored by each
  # FLowSOM metacluster
  
  blank_theme <- theme_minimal()+
    theme(
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  for(p in param){
    if(!is.character(dt[[p]])){
      print(paste("printing", p, sep = " "))
      g <- ggplot(dt, aes(x = umap1, y = umap2)) +
        geom_point(aes(color = dt[[p]])) + labs(title = p) + blank_theme +
        xlab("umap1") +
        ylab("umap2") +
        scale_colour_viridis_c(option  = "inferno") +
        labs(color="Expression")
      print(g)
      ggsave(paste0("umap_for_", p, ".pdf"))
    }
  }
}

### FACETING BY MARKER ####

printUmapFacet <- function(dt=dataTable, param=expr.params, legendKey=T, dotSize=0.01){
  # The dataTable must contain FlowSOM metaclusters
  # This function saves a png of the umap of choice, colored by each
  # FLowSOM metacluster
  blank_theme <- theme_minimal()+
    theme(
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  if(legendKey==F){
    blank_theme <- blank_theme + theme(legend.position = "none")
  }
  long.dt <- melt(dt, measure.vars = param)
  long.dt[value>1, value:=1]
  long.dt <- long.dt[sample(nrow(long.dt), nrow(long.dt)),]
  g <- ggplot(long.dt, aes(x = umap1, y = umap2)) +
    geom_point(aes(color = value), size = dotSize) + labs(title = "UMAP for each biomarker") + blank_theme +
    facet_wrap(~variable, nrow=5) +
    xlab("umap1") +
    ylab("umap2") +
    #lims(color = c(0, 1)) + 
    scale_colour_viridis_c(option  = "inferno") + 
    labs(color="Expression")
  print(g)
}
  
printUmapMetaclust_Facet <- function(umapDF=ss.data){
  
  
  umapDF <- as.data.frame(umapDF)
  
  umapDF <- umapDF[sample(nrow(umapDF), nrow(umapDF)),]
  
  
  ggplot(data = ss.data) +
  aes(x = umap1, y = umap2, color = metaCluster) +
  geom_point() +
  geom_point(size = 0.5) +
  scale_colour_viridis_d(option  = "cividis") +
  theme_minimal() +
  guides(color=guide_legend(nrow = 10, byrow = T)) +
  facet_wrap(vars(Tissue))

  
  
}


printUmapBymetaClust <- function(dt=dataTable){
  # The dataTable must contain FlowSOM metaclusters
  # This function saves a png of the umap of choice, colored by each
  # FLowSOM metacluster
  
  blank_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  
  wes <- wes_palette("Zissou1", max(dt$metaCluster), type = "continuous")
  

  g <- ggplot(data = dt) +
    aes(x = umap2, y = umap1, color = as.factor(metaCluster), size = totalEvents) +
    geom_point() + 
    labs(title = "UMAP overlaid with metaCluster ID") + 
    blank_theme +
    xlab("umap1") +
    ylab("umap2") +
    guides(color=guide_legend(nrow = 8, byrow = T)) +
    #scale_color_viridis_d() +
    #scale_color_manual(values =wes) +
    labs(color="Metacluster")

  print(g)
}

printUmapByParam <- function(dt=dataTable, param){
  # The dataTable must contain FlowSOM metaclusters
  # This function saves a png of the umap of choice, colored by each
  # FLowSOM metacluster
  
  blank_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  g <- ggplot(dt, aes(x = umap1, y = umap2)) +
    geom_point(aes(color = dt[[param]])) + labs(title = param) + blank_theme +
    xlab("umap1") +
    ylab("umap2") +
    scale_colour_viridis_c(option  = "inferno", limits = c(0,1)) +
    labs(color="Expression")
  print(g)
}

print_PD1_overlay <- function(dt, xparam, yparam){
  
  gt <- ggplot(dt, aes(x = dt[[xparam]], y = dt[[yparam]], color = PD.1)) + 
    geom_point() +
    xlab(xparam) +
    ylab(yparam)
  
  print(gt)
}

print_biaxial_facet <- function(dt, xparam, yparam, zparam=F){
  
  dt <- melt(dt, measure.vars = "Tissue")
  
  if(isFALSE(zparam)){
    gt <- ggplot(dt, aes(
      x = dt[[xparam]], 
      y = dt[[yparam]] 
      #color = PD.1
    )) + facet_wrap(~value) +
      geom_point() +
      xlab(xparam) +
      ylab(yparam) + theme_minimal()
  } else {
    gt <- ggplot(dt, aes(
      x = dt[[xparam]], 
      y = dt[[yparam]],
      color = dt[[zparam]]
    )) + facet_wrap(~value) +
      geom_point() +
      xlab(xparam) +
      scale_colour_viridis_c(option  = "inferno") + 
      labs(color=paste0(zparam, "\nExpression")) +
      ylab(yparam) + theme_minimal()
  }
  
  
  
  print(gt)
}


print_cell_event_number <- function(dt){
  ev_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  
  
  ggplot(data = dt) +
    aes(x = tumorType) +
    geom_bar(fill = '#737373') +
    labs(title = 'Total cell events per tumor type') +
    theme_minimal() +
    coord_flip() +
    ylim(0,12500) +
    xlab("") + ev_theme
}

print_cell_event_number_tissue <- function(dt){
  ev_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  
  
  ggplot(data = dt) +
    aes(x = reorder(Tissue, dt$TotalTissueEv)) +
    geom_bar(fill = '#737373') +
    labs(title = 'Segmented single cell data') +
    theme_minimal() +
    coord_flip() +
    ylim(0,8000) +
    xlab("") + ev_theme + ylab("Total Cell Events")
}

roundUpNice <- function(x, nice=seq(1, 10, 0.2)) {
  # use this to round up to a nice number
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

print_cell_event_number_run <- function(dt){
  ev_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
    )
  
  max.ev <- max(table(dt[,runNumb]))
  max.ev <- roundUpNice(max.ev)
  
  ggplot(data = dt) +
    aes(x = Run) +
    geom_bar(fill = '#737373') +
    labs(title = 'Total cell events per MIBI run') +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(name = "Total events",
                       labels = comma,
                       limits = c(0,max.ev)
                       #trans = log_trans()
    ) +
    xlab("") + ev_theme
}

print_cell_event_number_annotation <- function(dt){
  
  dt <- dt[metaClusterAnnotation!="none",]
  
  ev_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5)
      #axis.text.x = element_text(angle = 90)
    )
  
  max.ev <- max(table(dt[,metaClusterAnnotation]))
  max.ev <- roundUpNice(max.ev)
  
  ggplot(data = dt) +
    aes(x = metaClusterAnnotation) +
    geom_bar(fill = '#737373') +
    labs(title = 'Total cell events per annotated metacluster') +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(name = "Total events",
                       labels = comma,
                       limits = c(0,max.ev)
                       #trans = log_trans()
    ) +
    xlab("") + ev_theme
}


get_SOM_median_df <- function(dt, fa=non.exp.param){
  # variable
  # dt - data table with somCluster col from clustering, any 
  # output
  # df with median values for clusters, and the imageID and 
  # tumorID with the most cell events in that cluster
  fa=non.exp.param
  dt <- ss.data
  colnames(ss.data)
  dt <- dt[,3:ncol(dt)]
  expr.p <- colnames(dt)[!colnames(dt) %in% fa]

  med_names <- c(expr.p, "topTissue", "totalEvents", "somCluster")
  #add all the median values for the given som cluster to new data frame
  
  som_median_df <- data.frame(matrix(nrow = 0, ncol = length(med_names)))
  colnames(som_median_df) <- med_names
  current.df <- dt[somCluster==1,]
  
  current.df.expr <- current.df[,..expr.p]
  current.med <- colMedians(as.matrix(current.df.expr))
  som_median_df[1,1:length(current.med)] <- current.med
  
  #now add the top imageID, tissueSite, cancerDiagnosis, NULL cell ID (to be later deleted) and the somCluster ID
  tissues <- sort(table(current.df[,"Tissue"]))
  top.tissue <- tissues[length(tissues)]

  current.med <- c(names(top.tissue), nrow(current.df),  1)
  som_median_df[1,(length(expr.p)+1):ncol(som_median_df)] <- current.med
  
  

  for (i in 2:max(dt$somCluster)) {
    # here need to get the median value for each SOM cluster CSV
    
    current.df <- dt[somCluster==i,]
    current.df.expr <- current.df[,..expr.p]
    current.med <- colMedians(as.matrix(current.df.expr))
    som_median_df[i,1:length(current.med)] <- current.med
    
    #now add the top imageID, tissueSite, cancerDiagnosis, NULL cell ID (to be later deleted) and the somCluster ID
    tissues <- sort(table(current.df[,"Tissue"]))
    top.tissue <- tissues[length(tissues)]
    
    current.med <- c(names(top.tissue), nrow(current.df),  i)
    som_median_df[i,(length(expr.p)+1):ncol(som_median_df)] <- current.med
    
    }
  som_median_df$totalEvents <- as.numeric(as.character(som_median_df$totalEvents))
  som_median_df$somCluster <- as.numeric(as.character(som_median_df$somCluster))
  return(som_median_df)
}


somMedian_to_df <- function(somMedianMatrix, expressionParam){
  somMedianMatrix.num <- as.matrix(somMedianMatrix[,expressionParam])
  mode(somMedianMatrix.num) <- "numeric"
  
  somMedian.df <- as.data.frame(somMedianMatrix.num)
  somMedian.df$imageID <- somMedianMatrix$imageID
  somMedian.df$tumorType <- somMedianMatrix$tumorType
  somMedian.df$somCluster <- somMedianMatrix$somCluster
  return(somMedian.df)
}

print_SOM_median_df_superheat <- function(somMedianMatrix, col_to_plot, t = "SOM cluster heatplot", heatmap_scale = NA){

  somMedian.matrix <- as.matrix(somMedianMatrix[,col_to_plot])
  # print(summary(somMedian.matrix))
  if(is.na(heatmap_scale)){
    superheat(
      somMedian.matrix, 
      scale = F,
      row.dendrogram = T,
      col.dendrogram = T,
      left.label.size = 0.2,
      left.label.text.size = 2,
      bottom.label.size = 0.2,
      bottom.label.text.size = 5,
      bottom.label.text.angle = 90,
      left.label = "none",
      heat.pal.values = c(0, 0.1, 0.2, 1),
      title = t
      )
  }else{
    superheat(
      somMedian.matrix, 
      scale = F,
      row.dendrogram = T,
      col.dendrogram = T,
      left.label.size = 0.2,
      left.label.text.size = 2,
      bottom.label.size = 0.2,
      bottom.label.text.size = 4,
      bottom.label.text.angle = 90,
      left.label = "none",
      heat.lim = heatmap_scale,
      heat.pal.values = c(0, 0.1, 0.2, 1),
      title = t
    )
  }
}

print_SOM_median_df_superheat_metaC <- function(somMedianMatrix, 
                                                col_to_plot, 
                                                t = "SOM cluster heatplot", 
                                                heatmap_scale = NA
                                                ){
  
  somMedian.matrix <- as.matrix(somMedianMatrix[,..col_to_plot])
  # print(summary(somMedian.matrix))
  if(is.na(heatmap_scale)){
    superheat(
      somMedian.matrix, 
      scale = F,
      row.dendrogram = F,
      col.dendrogram = T,
      left.label.size = 0.2,
      left.label.text.size = 2,
      bottom.label.size = 0.2,
      bottom.label.text.size = 5,
      bottom.label.text.angle = 90,
      left.label = "none",
      heat.pal.values = c(0, 0.1, 0.2, 1),
      title = t,
      membership.rows = somMedianMatrix$metaCluster,
      smooth.heat = TRUE,
      
      grid.hline = T
    )
  }else{
    superheat(
      somMedian.matrix, 
      scale = F,
      row.dendrogram = T,
      col.dendrogram = T,
      left.label.size = 0.2,
      left.label.text.size = 2,
      bottom.label.size = 0.2,
      bottom.label.text.size = 4,
      bottom.label.text.angle = 90,
      left.label = "none",
      heat.lim = heatmap_scale,
      heat.pal.values = c(0, 0.1, 0.2, 1),
      title = t,
      membership.rows = somMedianMatrix$metaCluster
    )
  }
}


print_density <- function(dt, xparam, yparam){
  
  gt <- ggplot(dt, aes(x = dt[[xparam]], y = dt[[yparam]])) + 
    geom_density_2d() +
    xlab(xparam) +
    ylab(yparam)
  
  print(gt)
}

run_esquisser_gui <- function(){
  esquisse::esquisser()
}

somCluster <- function(dt, som_ch, sq_length = 10) {
  # Clusters with a SOM
  # Inputs:
  #   dt - data.table with all of the data to be clustered
  #   channels - vector of channel names
  # Outputs:
  #   dt - with added column
  dt.slim <- dt[,..som_ch]
  dt.slim <- as.matrix(dt.slim)
  set.seed(666)
  som.out <- SOM(dt.slim, xdim=sq_length, ydim=sq_length, distf=2) # 100 clusters using euclidean distance
  SOM
  return(som.out$mapping[,1])
}

metaCluster_consensus <- function(dt, som_parameters = cluster.markers, k_val){
  s = 403 # set seed
  #dt <- scale.somMedians
  cluster.mat <- as.matrix(dt[,som_parameters])
  dt$metaCluster <- metaClustering_consensus(cluster.mat, k = k_val, seed = s)
  
  return(dt)
}

get_metacluster_num <- function(dt, som_median_tabl){
  # dt <- ss.data
  # som_median_tabl <- somMedians.consensus
  dt[,metaCluster:=0]
  for(i in 1:max(dt$somCluster)){
    c <- som_median_tabl[som_median_tabl$somCluster==i,]
    dt[somCluster==i,metaCluster:=c$metaCluster]
  }
  dt$somCluster <- as.factor(as.numeric(dt$somCluster))
  dt$metaCluster <- as.factor(as.numeric(dt$metaCluster))
  return(dt)
}

get_image_total_ev <- function(dt, som_median_tabl){
  dt <- master_dt.score.t.som.scale
  som_median_tabl <- scale.somMedians
  dt[,totalEvents:=0]
  for(i in unique(dt$imageID)){
    print(i)
    c <- dt[dt$imageID==i,]
    dt[imageID==i,totalEvents:=nrow(c)]
  }
  return(dt)
}

random_sample <- function(dt, n){
  randDt <- dplyr::sample_n(dt, n)
  return(data.table(randDt))
}

scaleData <- function(dt, fa) {
  # scales each run individually
  # Inputs:
  #   dt - data.table
  #   fa - character vector of channels not to scale
  # Outputs:
  #   dt - data.table
  print("Scaling data")
  channels <- colnames(dt)[!colnames(dt) %in% fa]
  print(channels)
  refs <- mapply(quantile, x=dt[, channels, with=F], MARGIN=2, probs=c(0.999), na.rm=T)
  
  print(refs)
  for (i in seq(channels)) {
    if(refs[i] < 0.05){
      print(refs[i])
    }
    dt[, channels[i]:=dt[, channels[i], with=F] / refs[i]]
  }
  return(dt)
}

data_summary <- function(x) {
  # This function is used to display median +/- SD for violin plots
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


printParamViolinListPerMetaClustReOrder <- function(dataTable, paramList, totalClusterNumber = total.metacluster){
  # This function reorders the parameters (highest to lowest from top to bottom)
  # input 
  # - dataTable - data table of single-cell events
  # - paramList - vector of str containing the parameters to plot
  # - totalClusterNumber - int value of the total number of FlowSOM metaclusters
  # return
  # - prints ggplot object
  
  blank_theme <- theme_minimal()+
    theme(
      #panel.grid=element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title=element_text(size=14, face="bold", hjust = 0.5),
      legend.position = "none"
    )
  
  dtl <- melt(dataTable, measure.vars = paramList)
  ev.table <- data.frame(table(dataTable[,metaCluster]))
  
  for(i in 1:totalClusterNumber){
    # print a plot of params in paramList for each metacluster
    g <- ggplot(dtl[metaCluster == i,], aes(x = reorder(variable, value), y = value, group = variable, color = NULL)) + 
      blank_theme +
      geom_violin(scale = "width") + 
      stat_summary(fun.data=data_summary, geom="pointrange") + # plots a dot for the median and line for +/-SD
      labs(title = paste0("FlowSOM Cluter Parameter Expression\nCluster: ", i, ", Total Cells: ", ev.table[i,2])) +
      ylab("Expression asinh (x/5)") +
      xlab("") +
      coord_flip()
    print(g)
    ggsave(filename = paste0("cluster_", i, "_violin.png"), height = 5, width = 6, units = "in")
    #ggsave_both_PSandPNG(Plot = g, fileName = paste0("cluster_", i, "_violin"))
  }
  print("done.")
}

ggsave_both_PSandPNG <- function(Plot, fileName, height = 5, width = 5, units = "in"){
  ggsave(plot, filename = paste0(fileName, ".ps"), height = height, width = width, units = units)
  ggsave(plot, filename = paste0(fileName, ".png"), height = height, width = width, units = units)
}
