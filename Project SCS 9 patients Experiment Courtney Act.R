### 2019-2-19

# Project SCS 9 patients Experiment Courtney Act

#-----------------------------------------------------------------------------------------
### load data 

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 9 patients/Experiment Courtney Act")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.1

### load packages
library(Seurat)        # V2.3.4
library(dplyr)         # v0.7.6
library(DESeq2)

# object
seuset <- readRDS("seuset 9 patients.RDS")
fibrous <- readRDS("Fibrous patients only.RDS")
atheromatous <- readRDS("Atheromatous patients only2.RDS")

#-----------------------------------------------------------------------------------------
# subset based on phenotype
meta.data <- seuset@meta.data
fibrous.cells <- rownames(meta.data[meta.data$Phenotype == "Fibrous",])

fibrous <- SubsetData(object = seuset,
                      cells.use = fibrous.cells)

atheromatous.cells <- rownames(meta.data[meta.data$Phenotype == "Atheromatous" | meta.data$Phenotype == "Fibro-atheromatous",])

atheromatous <- SubsetData(object = seuset,
                           cells.use = atheromatous.cells)

dim(fibrous@meta.data)
dim(atheromatous@meta.data)

#pdf(paste(result.folder, "/", Sys.Date(), " tSNE plots.pdf", sep = ""))
TSNEPlot(fibrous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "all")
TSNEPlot(fibrous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "fibrous")
TSNEPlot(atheromatous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "atheromatous")
dev.off()

# downstream Functions 
downstream.functions <- function(object) {
  object <- FindVariableGenes(object, do.plot = F)
  print(head(object@var.genes, n = 20))
  
  # Run PCA
  object <- RunPCA(object = object,
                   do.print = F)
  
  print(PCElbowPlot(object, num.pc = 20))
  
  return(object)
}


# Find clusters
find.clusters <- function(object) { 
  Random.Seed <- set.seed(53)
  object <- FindClusters(object = object, 
                             reduction.type = "pca", 
                             resolution = 1.2, 
                             dims.use = 1:10, 
                             force.recalc = T,
                             random.seed = 53 )
  
  
  object <- RunTSNE(object,
                        reduction.use = "pca",
                        dims.use = 1:10,
                        seed.use = Random.Seed)
  
  #return(object)
}

#--------------------------------------------------------------------------------------
### Fibrous
## PCA etc
fibrous <- downstream.functions(fibrous)

## find clusters
# res 1.2
fibrous <- find.clusters(fibrous)

to.plot <- c("MYH11", "CD34", "CD14", "CD68","CD3E","CD79A", "KIT")

pdf(paste(result.folder, "/", Sys.Date(), " tSNE plots Fibrous.pdf", sep = ""))
TSNEPlot(fibrous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")
for (i in to.plot) {
  FeaturePlot(object = fibrous,
              features.plot = i)
}
dev.off()

table(fibrous@meta.data$res.1.2)

# save
saveRDS(fibrous, file = "Fibrous patients only.RDS")

#--------------------------------------------------------------------------------------
### atheromatous
## PCA etc
atheromatous <- downstream.functions(atheromatous)

## find clusters
# res 1.2
atheromatous <- find.clusters(atheromatous)

to.plot <- c("MYH11", "CD34", "CD14", "CD68","CD3E","CD79A", "KIT")

pdf(paste(result.folder, "/", Sys.Date(), " tSNE plots Atheromatous.pdf", sep = ""))
TSNEPlot(atheromatous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")
for (i in to.plot) {
  FeaturePlot(object = atheromatous,
              features.plot = i)
}
dev.off()

table(atheromatous@meta.data$res.1.2)

# save
saveRDS(atheromatous, file = "Atheromatous patients only2.RDS")

#--------------------------------------------------------------------------------------
Lianne.genes <- c("EDN1", "EDNRA", "EDNRB")

pdf(paste(result.folder, "/", Sys.Date(), " feature plots Atheromatous.pdf", sep = ""))
TSNEPlot(atheromatous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "Atheromatous")
for (i in Lianne.genes) {
  FeaturePlot(object = atheromatous,
              features.plot = i)
}
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), " feature plots Fibrous.pdf", sep = ""))
TSNEPlot(fibrous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "Fibrous")
for (i in Lianne.genes) {
  FeaturePlot(object = fibrous,
              features.plot = i)
}
dev.off()


### split athermomatous in two sexes
female.cells <- rownames(meta.data[meta.data$AE == 4459 | meta.data$AE == 4470,])
head(female.cells)

male.cells <- rownames(meta.data[meta.data$AE != 4459 | meta.data$AE != 4470,])

female <- SubsetData(object = atheromatous,
                     cells.use = female.cells)

male <- SubsetData(object = atheromatous,
                   cells.use = male.cells)

TSNEPlot(female, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "female atheromatous")
TSNEPlot(male, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "male atheromatous")

# not enough power
rm(female, male, female.cells, male.cells)


#--------------------------------------------------------------------------------------
# fibrous vs atheromatous epithelium and smcs
# select cells of interest
get.cells <- function(object, ident.nr) {
  # transform object@raw.data into matrix
  current.data <- as.matrix(object@raw.data)
  # get all colums where the colnames correspond to the cells with chosen ident number
  current.matrix <- current.data[ ,rownames(object@meta.data[object@meta.data$res.1.2 == ident.nr,])]
  return(current.matrix)
  # delete from environment
  rm(current.data, current.matrix)
}

# create dataset per celltype per phenotype
# since there are no repliates, we cannot use DESeq
fib.vs.athero <- data.frame(endo.fibrous = rowSums(get.cells(fibrous, 7)),
                            endo.atheromatous = rowSums(get.cells(atheromatous, 8)),
                            smc.fibrous = rowSums(get.cells(fibrous, 9)),
                            smc.atheromatous = rowSums(get.cells(atheromatous, 7)))

# transform
fib.vs.athero <- log2(fib.vs.athero + 0.1)

## check voor genes with highest variance
plot(fib.vs.athero[,1:2])
plot(fib.vs.athero[,3:4])
cor(fib.vs.athero)


current.violin.input <- reshape2::melt(fib.vs.athero[,1:2], id.vars = NULL)
ggplot(data = current.violin.input , aes(x = variable, y = value)) +
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5, colour = "darkgrey") + 
  geom_boxplot(width=0.08, colour = "red")

fib.vs.athero.pseudobulk <- fib.vs.athero

#--------------------------------------------------------------------------------------
## DESeq between fibrous patients and atheromatous patients

# prepare
input.meta.data <- read.table('Project SCS 9 patients meta data.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
names(input.meta.data)[1] <- "AE"
input.meta.data$ID <- paste(input.meta.data$AE, ".P",input.meta.data$Plate, sep = "")

patient.phenotype <- distinct(input.meta.data[,c(1,3)])
patient.phenotype$Phenotype2 <- c("Atheromatous", rep("Fibrous",2),"Atheromatous", "Fibrous", rep("Atheromatous",4))

# function to get the correct celltype per patient
get.cells <- function(object, AE.nr) {
 # select correct celltype 
  if (deparse(substitute(object)) == "fibrous") { 
    ident.nr <- 7
  } else  {
    ident.nr <- 8
  }
  
  # transform object@raw.data into matrix & get (raw) data
  current.matrix <- as.matrix(object@raw.data)
  #current.matrix <- as.matrix(object@data)
  # select cells that correspond to the correct ident number and AE number
  current.cells <- object@meta.data[object@meta.data$res.1.2 == ident.nr & object@meta.data$AE == AE.nr,]
  current.cells <- rownames(current.cells)
  
  # print length current cells
  print(length(current.cells))
  
  # get the correct cells and sum 
  current.matrix <- rowSums(current.matrix[,current.cells, drop = F])
 
  return(current.matrix)
  
}

fib.vs.athero <- data.frame(matrix(nrow = dim(seuset@raw.data)[1], ncol = 0))
colnames.fib.vs.athero <- c()
for (AE in patient.phenotype$AE) {
  # get phenotype
  if (patient.phenotype[patient.phenotype$AE == AE,2] == "Fibrous") {
    # get the cells
    correct.cells <- get.cells(fibrous, AE)
    name <- "fib"
  } else {
    correct.cells <- get.cells(atheromatous, AE)
    name <- "athero"
  }
  
  # bind to matrix
  fib.vs.athero <- cbind(fib.vs.athero, correct.cells)
  # save name
  colnames.fib.vs.athero <- c(colnames.fib.vs.athero, paste(AE,".",name, sep = ""))
  
  #remove from environment
  rm(correct.cells,AE, name)
}

# give colnames for clarity
names(fib.vs.athero) <- colnames.fib.vs.athero
head(fib.vs.athero)
# check
colSums(fib.vs.athero)
dim(fib.vs.athero)

### DESeq2
# choose columns to proceed with
DESeq.fib.vs.athero <- fib.vs.athero[,c(2:7)]
head(DESeq.fib.vs.athero, n =1)
coldata <- patient.phenotype[c(2:7),]

# round data
DESeq.fib.vs.athero <- round(DESeq.fib.vs.athero)

dds <- DESeqDataSetFromMatrix(countData =  DESeq.fib.vs.athero,
                              colData = coldata,
                              design = ~ Phenotype2)


dds <- DESeq(dds)                               
pre.results <- results(dds)
head(pre.results)
plotMA(pre.results)

results <- data.frame(pre.results)
results <- results[!is.na(results$pvalue),]

# results <- results[order(results$padj),]
# head(results)

results <- results[order(results$pvalue),]
head(results)
sum(results$pvalue <= 0.05)
write.table(results[results$pvalue <= 0.05,], "SMCs Fib vs Ath from raw.data.txt", sep = "\t", col.names = T, row.names = T)


genes.to.plot <- c("THBS1")

ggplot(data = fib.vs.athero.pseudobulk[,1:2], aes(x = fib.vs.athero.pseudobulk[,1], y = fib.vs.athero.pseudobulk[,2])) +
  geom_point() + 
  geom_point(data = fib.vs.athero.pseudobulk[genes.to.plot,1:2], aes(x = fib.vs.athero.pseudobulk[genes.to.plot,1], y = fib.vs.athero.pseudobulk[genes.to.plot,2]), col = "red", size = 4) +
  ggtitle("Fibrous vs Atheromatous") + xlab("Fibrous") + ylab("Atheromatous") +
  geom_text(data = fib.vs.athero.pseudobulk[genes.to.plot, 1:2], aes(x = fib.vs.athero.pseudobulk[genes.to.plot,1], y = fib.vs.athero.pseudobulk[genes.to.plot,2]),label = genes.to.plot, hjust=0, vjust=0)


