# Test reproducibility of Seurat example
# From Seurat tutorial website: satijalab.org/seurat/articles/pbmc3k_tutorial.html
# Author: Candace Liu
# Date: 8/15/22

library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(FlowSOM)
library(ConsensusClusterPlus)

pbmc.data = Read10X(data.dir="filtered_gene_bc_matrices/hg19")
pbmc = CreateSeuratObject(counts=pbmc.data, project="pbmc3k", min.cells=3, min.features=200)

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern="^MT-")
VlnPlot(pbmc, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)

plot1 = FeatureScatter(pbmc, feature1="nCount_RNA", feature2="percent.mt")
plot2 = FeatureScatter(pbmc, feature1="nCount_RNA", feature2="nFeature_RNA")
plot1 + plot2

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
pbmc = NormalizeData(pbmc, normalization.method="LogNormalize", scale.factor=10000)

# Identification of highly variable features
pbmc = FindVariableFeatures(pbmc, selection.method="vst", nfeatures=2000)
top10 = head(VariableFeatures(pbmc), 10)

plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot=plot1, points=top10, repel=TRUE)
plot1 + plot2

# Scale the data
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features=all.genes)

# Perform linear dimensional reduction
pbmc = RunPCA(pbmc, features=VariableFeatures(object=pbmc))
print(pbmc[["pca"]], dims=1:5, nfeatures=5)
VizDimLoadings(pbmc, dims=1:2, reduction="pca")
DimPlot(pbmc, reduction="pca")
DimHeatmap(pbmc, dims=1, cells=500, balanced=TRUE)

# Determine the dimensionality of the dataset
pbmc = JackStraw(pbmc, num.replicate=100)
pbmc = ScoreJackStraw(pbmc, dims=1:20)
JackStrawPlot(pbmc, dims=1:15)
ElbowPlot(pbmc)


## Test reproducibilty of clustering
all_seeds = c(2022,32,115,6,22)
k = 15 #number of clusters

one_seed <- function(seed) {
  cell_dat = data.table(pbmc$pca@cell.embeddings[,1:10])
  set.seed(seed)
  output = SOM(as.matrix(cell_dat))
  cell_dat$cluster = output$mapping[,1]
  mean_dat = cell_dat[,lapply(.SD,mean), by=cluster]
  mat_dat = data.frame(mean_dat[,1:10])
  # Z-score the columns
  mat_dat = scale(mat_dat)
  # Get hierarchical clusters
  consensus = ConsensusClusterPlus(t(mat_dat), maxK=k, seed=seed)
  clust_to_hclust = consensus[[k]]$consensusClass
  names(clust_to_hclust) = mean_dat$cluster
  cell_dat$hCluster = clust_to_hclust[as.character(cell_dat$cluster)]
  return(cell_dat$hCluster)
}
all_seed_output = lapply(all_seeds, one_seed)
dt = data.table(rep1=as.integer(all_seed_output[[1]]),
                rep2=as.integer(all_seed_output[[2]]),
                rep3=as.integer(all_seed_output[[3]]),
                rep4=as.integer(all_seed_output[[4]]),
                rep5=as.integer(all_seed_output[[5]]))
fwrite(dt,"seurat_flowsom.csv")


