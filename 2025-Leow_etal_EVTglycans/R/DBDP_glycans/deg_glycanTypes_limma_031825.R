# Get pval,fc by differential expression of glycans 
# Using edgeR to input data and limma for differential expression analysis
# Author: Hadeesha - adapted by Ke
# Reference: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# Date: 03/17/25

#--------------------------------
# Load packages/functions
#--------------------------------
library(limma)
library(edgeR)
library(tidyverse)
library(data.table)

#--------------------------------
# Load data
#--------------------------------
data <- read_csv('data/DBDP_transition/MALDI/library_matched/DBDP_transition_glyRatio_glycanTypes_031825.csv')

#--------------------------------
# Format data for DGElist (edgeR object)
#--------------------------------
data_filter <- data %>% 
  filter(Region == "DB" | Region == "DP") 
data_int <- data_filter %>% select(R_fucosylated:R_bisecting) 
data_int <- data_int*1000 #scaled data cause limma doesn't work well for intensities between 0 and 1
data_pheno <- data_filter[,c("mask","Region")] 

counts <- as.data.frame(data_int) %>% t()
colnames(counts) <- data_pheno$mask

#--------------------------------
# Load into DGElist
#--------------------------------
### the variable names might not make much sense since I just used them out of their package lol
d <- DGEList(counts, group = data_pheno$Region)
group <- d$samples$group
# d <- calcNormFactors(d)

# snames <-colnames(counts)
# cultivar <- gsub("_.*$","",snames)
# group <- interaction(cultivar)

plotMDS(d, col = as.numeric(group)) ## check how the groups separate by MDS

# library(factoextra)
# library(FactoMineR)
# 
# pca.raw.d <- log2(d$counts+0.5)
# pca.d <- PCA(t(pca.raw.d),graph = F)
# fviz_pca_ind(pca.d, col.ind = group) ## check how the groups separate by MDS

#--------------------------------
# Run differential expression analysis
#--------------------------------
# Specify the model to be fitted
mm <- model.matrix(~ 0 + Region, data=data_pheno)

# 3. Voom transformation and calculation of variance weights
voom.y.d <- voom(d, mm, plot = T)

# 4. Fitting linear models in limma
# lmFit fits a linear model using weighted least squares for each gene
fit <- lmFit(voom.y.d, mm)
coef.fit <- fit$coefficients
head(coef(fit))


#initialize final df
# de <- data.frame()

# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
# Specify which groups to compare:
contr <- makeContrasts(RegionDB - RegionDP,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("DB_DP",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)] %>% 
  mutate(composition = gsub("R_", "", Glycan))

# de <- top.table %>% 
#   inner_join(.,coef.fit %>% select(Glycan, enzymeControl,enzymeLacNAcase)) %>% 
#   select(Glycan, composition, enzymeControl, enzymeLacNAcase, everything())

#--------------------------------
# write de to csv
#--------------------------------
# gly_tab <- read_csv("data/matched_peaks.csv") %>% 
#   na.omit()
# de_export <- left_join(de, gly_tab[,c("composition","lib_mz", "peak", "mass_error")])
  
write.csv(top.table, "data/DBDP_transition/MALDI/library_matched/DE_pval_glycanTypes_DBDP.csv", row.names = F)
# write.csv(coef.fit, "data/limmafit_glycans_popeAttaching_sePEvsPTL.csv")
