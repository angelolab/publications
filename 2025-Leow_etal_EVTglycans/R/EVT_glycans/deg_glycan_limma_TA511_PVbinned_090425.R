# Get pval,fc by differential expression of glycans in evts
# Using edgeR to input data and limma for differential expression analysis
# Author: Hadeesha - adapted by Ke
# Reference: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# Date: 11/14/22

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
#glycan data - not scaled
data <- read_csv('data/MALDI_IF_EVT_glycans/IFmaldi_TA511_mean_gly_permask_evtVctFv_wMeta.csv')

data_update <- data %>% 
  mutate(
    cell_anno_update = case_when(
      cell_anno %in% c("fv", "vct") ~ "PV",
      TRUE ~ cell_anno
    )
  ) %>%
  filter(cell_anno_update %in% c("PV", "EVT_A", "EVT_I", "EVT_E")) %>% 
  select(cell_anno_update, everything())

write.csv(data_update, "data/MALDI_IF_EVT_glycans/IFmaldi_TA511_mean_gly_permask_PVbinned.csv", row.names = F)

#--------------------------------
# Format data for DGElist (edgeR object)
#--------------------------------
data_int <- data_update %>% select_at(vars(contains("H"))) 
data_int <- data_int*1000 #scaled data cause limma doesn't work well for intensities between 0 and 1

data_pheno <- data_update %>% select(cell_anno_update, ga, patient_id, maskname) 

counts <- as.data.frame(data_int) %>% t()
colnames(counts) <- data_pheno$maskname

#--------------------------------
# Load into DGElist
#--------------------------------
### the variable names might not make much sense since I just used them out of their package lol
d <- DGEList(counts, group = data_pheno$cell_anno_update)
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
mm <- model.matrix(~ 0 + group)

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
contr <- makeContrasts(groupEVT_A - groupEVT_E,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("groupEVTA_groupEVTE",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)]

de <- top.table

#--------------------------------
# Repeat DE analysis for all pairs of groups
#--------------------------------
#groupEVT_A - groupEVT_I
# Specify which groups to compare:
contr <- makeContrasts(groupEVT_A - groupEVT_I,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("groupEVTA_groupEVTI",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)]

de <- inner_join(de, top.table, by = "Glycan")

################################
#groupEVT_A - groupPV
# Specify which groups to compare:
contr <- makeContrasts(groupEVT_A - groupPV,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("groupEVTA_groupPV",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)]

de <- inner_join(de, top.table, by = "Glycan")

################################
#groupEVT_E - groupEVT_I
# Specify which groups to compare:
contr <- makeContrasts(groupEVT_E - groupEVT_I,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("groupEVTE_groupEVTI",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)]

de <- inner_join(de, top.table, by = "Glycan")


################################
#groupEVT_E - groupPV
# Specify which groups to compare:
contr <- makeContrasts(groupEVT_E - groupPV,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("groupEVTE_groupPV",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)]

de <- inner_join(de, top.table, by = "Glycan")


################################
#groupEVT_I - groupPV
# Specify which groups to compare:
contr <- makeContrasts(groupEVT_I - groupPV,levels = colnames(coef(fit))) # the two groups you would like to contrast

# Estimate contrast for each gene and Empirical Bayes smoothing of standard errors
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5)

#save results to final df
newcol <- paste("groupEVTI_groupPV",colnames(top.table), sep = "_")
colnames(top.table) <- newcol

top.table$Glycan <- rownames(top.table)
top.table <- top.table[,c(7,1,4,5)]

de <- inner_join(de, top.table, by = "Glycan")



#--------------------------------
# write de to csv
#--------------------------------
write.csv(de, "data/MALDI_IF_EVT_glycans/de_glycans_TA511_PVbinned.csv", row.names = F)

# dat <- read_csv("data/de_glycans_TA511.csv")
