### 2019-2-20

# Project SCS 9 patients Experiment Detox

# pseudotime for SMCs

#-----------------------------------------------------------------------------------------
### load data 
# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 9 patients/Experiment Detox")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.1

### load packages
library(Seurat)        # V2.3.4
library(plyr)          # v1.8.4
library(dplyr)         # v0.7.6
library(tidyr)         # v0.8.2
library(slingshot)     # v1.1.1


# object
seuset <- readRDS("SMCs subset 9 patients.RDS")

#-----------------------------------------------------------------------------------------
# ### get SMCs
# # original data
# seuset <- readRDS("seuset 9 patients.RDS")
# 
# # subset
# seubset <- SubsetData(seuset,
#                       ident.use = 10)
# 
# seubset <- FindVariableGenes(seubset, do.plot = F)
# 
# seubset <- RunPCA(seubset, do.print = F)
# 
# seubset <- FindClusters(object = seubset,
#                         reduction.type = "pca",
#                         resolution = 1,
#                         dims.use = 1:10,
#                         force.recalc = T,
#                         random.seed = 91 )
# 
# seubset <- RunTSNE(seubset, dims.use = 1:10)
# 
# # save subset
# saveRDS(seubset, "SMCs subset 9 patients.RDS")
# 
# # remove from enviornment
# rm(seuset)
# # change
# seuset <- seubset
# # remove
# rm(seubset)

#-----------------------------------------------------------------------------------------
### Run slingshot 
## http://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/slingshot.html

# The minimal input to slingshot is a matrix representing the cells in a reduced-dimensional space and a vector of cluster labels.
# seuset@raw.data / seuset@data
# seuset@meta$res1 ?

## Dataset
library(SingleCellExperiment)
library(splatter)
# get counts
counts <- as.matrix(seuset@data)
# estimate simulation parameters
params <- splatEstimate(counts)
# simulate
sim <- splatSimulate(params, method = "paths")

## Upstream analysis
# Gene filtering
# filter genes down to potential cell-type markers
# at least M reads in at least N cells
geneFilter <- apply(assays(sim)$counts,1,function(x){
  sum(x >= 10) >= 4
})
sim <- sim[geneFilter, ]

# Normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)

# Dimensionality Reduction
library(destiny)

pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = topo.colors(100)[sim$Step], pch=16, asp = 1)

dm <- DiffusionMap(t(log1p(assays(sim)$norm)))

rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)

plot(rd2, col = topo.colors(100)[sim$Step], pch=16, asp = 1)
reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)

## clustering cells
library(mclust, quietly = TRUE)

cl1 <- Mclust(rd1)$classification
colData(sim)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

## using slingshot
# To run slingshot with the dimensionality reduction produced by PCA and cluster labels identified by Gaussian mixutre modeling, we would do the following:
sce <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')

summary(sce$slingPseudotime_1)

# Below, we visuzalize the inferred lineage for the splatter data with points colored by pseudotime.
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

# We can also see how the lineage structure was intially estimated by the cluster-based minimum spanning tree by using the type argument.
# plot lineage
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages')

## Downstream Analysis
# Identifying temporally expressed genes
# After running slingshot, an interesting next step may be to find genes that change their expression over the course of development. W
# We demonstrate one possible method for this type of analysis on the 1,000 most variable genes.
require(gam)

t <- sce$slingPseudotime_1

# for time, only look at the 1,000 most variable genes
Y <- log1p(assays(sim)$norm)
var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
Y <- Y[var1K,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# We can then pick out the top genes based on p-value and visualize their expression over developmental time with a heatmap.
require(clusterExperiment)

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(sim)$norm[rownames(assays(sim)$norm) %in% topgenes, 
                             order(t, na.last = NA)]
heatclus <- sce$GMM[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",visualizeData = 'transformed')

## Detailed Slingshot Functionality
# Identifying global lineage structure
# now need dimensionally reduced data and vector of cluster numbers
# omit starting cluster variable
lin1 <- getLineages(rd, cl)
