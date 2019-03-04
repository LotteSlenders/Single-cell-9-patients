### 2019-1-29

# Project SCS 9 patients Experiment Alyssa Edwards

#-----------------------------------------------------------------------------------------
### load data 

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 9 patients/Experiment Alyssa Edwards")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.1

### load packages
library(Seurat)        # V2.3.4
library(dplyr)         # v0.7.6
library(org.Hs.eg.db)

#-----------------------------------------------------------------------------------------
### prep counts
# load meta data
meta.data <- read.table('Project SCS 9 patients meta data.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
names(meta.data)[1] <- "AE"
meta.data$ID <- paste(meta.data$AE, ".P",meta.data$Plate, sep = "")
sum(duplicated(meta.data$ID))
#View(meta.data)

# check files
all.files <- list.files(pattern="*.tsv")
all.files
all.files %in% meta.data$File
length(all.files) == dim(meta.data)[1]

### read in files and QC
## update symbols
#Set number of cells / plate.
n = 384
update.symbols <- function(my.counts.table = my.counts.table, n = n){
  #Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
  #Retrieve the mapping
  gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(my.counts.table), columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
  
  #Keep only one mapping per alias
  d <- duplicated(gene.info$ALIAS)
  gene.info <- gene.info[!d,]
  
  #Add mappings to the table
  my.counts.table$SYMBOL <- gene.info$SYMBOL
  
  #Remove non-mappings (old LOCs and stuff that are not availble anymore)
  na <- is.na(my.counts.table$SYMBOL)
  my.counts.table <- my.counts.table[!na,]
  
  #Keep only highest expressed SYMBOL if multiple old aliases are now merged
  #Get maximum number of dups per gene
  num_dups <- max(table(my.counts.table$SYMBOL))
  
  #Loop until only one entry per gene left
  while(num_dups > 1){
    #First retrieve the list of duplicated symbols
    o <- order(my.counts.table$SYMBOL)
    my.counts.table <- my.counts.table[o,]
    d <- duplicated(my.counts.table$SYMBOL)
    
    #Then keep only the highest expressed one (this part will fail if we happen upon 2 rows with equal expression, but for now that hasn't happened yet so meh ;)
    h <- rowSums(my.counts.table[d,1:n]) > rowSums(my.counts.table[which(d==T)-1,1:n])
    h <- c(h, rowSums(my.counts.table[which(d==T)-1,1:n]) > rowSums(my.counts.table[d,1:n]))
    h <- h[which(!h)]
    my.counts.table <- my.counts.table[!row.names(my.counts.table) %in% names(h),]
    
    num_dups <- num_dups - 1 #One dup removed
  }
  
  #Overwrite old symbols in row names and clean up
  row.names(my.counts.table) <- my.counts.table$SYMBOL
  my.counts.table$SYMBOL <- NULL
  
  return(my.counts.table)
}

## read in and process
seurat.object.list <- list(c)
unwanted.genes <- c("^ERCC", "^MT-","^RPL", "^RPS","^UGDH-AS1$","^PGM2P2$", "^LOC100131257$", "^KCNQ1OT1$", 
                    "^MALAT1$", "^PGM5P2$", "^MAB21L3$","^EEF1A1$","^MALAT1$")

for (df in seq_along(all.files)) {
  current.df <- read.delim(all.files[df])
      
  # fix gene names and rownames
  current.df$GENEID <- sub("__.*", "" ,  current.df$GENEID)
  rownames(current.df) <- current.df$GENEID
  current.df <- current.df[,-1]
  # remove unwanted genes
  current.df <- current.df[grep(paste(unwanted.genes, collapse = "|"),rownames(current.df),invert=TRUE),]
  
  # match AE and plate number to file
  current.meta.data <- meta.data[meta.data$File == all.files[df],]
  colnames(current.df) <- paste(current.meta.data$ID, ".", 1:length(current.df), sep = "")
  
  # update gene names
  current.df <- update.symbols(current.df, n)
  
  ## create object
  current.object <- CreateSeuratObject(current.df, project = current.meta.data$ID)
  # add corresponding meta data
  current.object@meta.data <- cbind(current.object@meta.data, current.meta.data)
  # add to list
  seurat.object.list[[df]] <- current.object
  
  if (df == round(0.5*length(all.files))) {
    print("Halfway done")
  } 
  
  if (df == length(all.files)){
    print(":)")
  }

}

# somehow creating the objects prior occurs with a loss of cells for some, I do not know the cause
View(head(current.df))
head(current.object@meta.data)
tail(current.object@meta.data)
dim(current.object@meta.data)
dim(current.object@raw.data)
dim(current.object@data)

#-----------------------------------------------------------------------------------------
### QC
# merge suerat objects from list
seuset <- seurat.object.list[[1]]
for (i in 2:length(x = seurat.object.list)) {
  # only filter genes with too low expression during last merge
  if (i == length(seurat.object.list)){
  seuset <- MergeSeurat(object1 = seuset, object2 = seurat.object.list[[i]], do.normalize = F,
                        min.cells = 10)
  } else {
    seuset <- MergeSeurat(object1 = seuset, object2 = seurat.object.list[[i]], do.normalize = F)
  }
}
# remove object list to save memory
rm(seurat.object.list)

dim(seuset@data)

# Filter cells with too low and too high gene counts
min(seuset@meta.data$nGene)
max(seuset@meta.data$nGene)
mean(seuset@meta.data$nGene)
median(seuset@meta.data$nGene)

seuset <- FilterCells(object = seuset,
                      subset.names = c("nGene"),
                      low.thresholds = 200,
                      high.thresholds = 10000
)
dim(seuset@data)

# Normalize data
seuset <- NormalizeData(object = seuset,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000
)

# scale data (regress out UMI)
seuset <- ScaleData(object = seuset, 
                    vars.to.regress = c("nUMI"))

# Find variable genes 
seuset <- FindVariableGenes(seuset, do.plot = F)
head(seuset@var.genes, n = 20)

# Run ICA
random.seed <- set.seed(99)
seuset <- RunICA(object = seuset,
                 ics.compute = 50,
                 genes.print = 1,
                 seed.use = random.seed)

# Run PCA
seuset <- RunPCA(object = seuset,
                 do.print = F)

PCElbowPlot(seuset, num.pc = 20)

#-----------------------------------------------------------------------------------------
### Clustering ICA
# seuset.ICA <- FindClusters(object = seuset, 
#                            reduction.type = "ica", 
#                            resolution = 1.7, 
#                            dims.use = 1:15, 
#                            force.recalc = T)
# 
# seuset.ICA <- RunTSNE(seuset.PCA,
#                       reduction.use = "ica",
#                       dims.use = 1:15
# )
# 
# TSNEPlot(seuset.ICA, do.label = T, pt.size = 1, label.size = 5, group.by = "ident")

#-----------------------------------------------------------------------------------------
### clustering PCA
random.seed <- set.seed(91)
seuset.PCA <- FindClusters(object = seuset, 
                          reduction.type = "pca", 
                          resolution = 1.7, 
                          dims.use = 1:10, # based on elbow, 7. But 10 is still okay.
                          force.recalc = T,
                          random.seed = 91 )


seuset.PCA <- RunTSNE(seuset.PCA,
                      reduction.use = "pca",
                      dims.use = 1:10,
                      seed.use = random.seed)

## check stuff
# how many cells per type
sort(table(seuset.PCA@meta.data$res.1.7))
# how many cells per plate
sort(table(seuset.PCA@meta.data$ID))
# how many cells per type per plate
table(seuset.PCA@meta.data$res.1.7, seuset.PCA@meta.data$ID)
# how many cells per patient
sort(table(seuset.PCA@meta.data$AE))

barplot(prop.table(x = table(seuset.PCA@ident, seuset.PCA@meta.data$AE)))
barplot(prop.table(x = table(seuset.PCA@ident, seuset.PCA@meta.data$ID)))

## Check tSNE
#pdf(paste(result.folder, "/", Sys.Date(), "Preliminary plots.pdf", sep = ""))
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")

# from prior knowledge
# endothelial cells
FeaturePlot(seuset.PCA, features.plot = c("EDN1"))
FeaturePlot(seuset.PCA, features.plot = c("EDNRA", "EDNRB"))
FeaturePlot(seuset.PCA, features.plot = c("CDH5", "PECAM1"))
FeaturePlot(seuset.PCA, features.plot = c("CD34"))
# SMC
FeaturePlot(seuset.PCA, features.plot = c("MYH11"))
FeaturePlot(seuset.PCA, features.plot = c("LGALS3", "ACTA2"))
# macrophages
FeaturePlot(seuset.PCA, features.plot = c("CD14", "CD68"))
FeaturePlot(seuset.PCA, features.plot = c("CD36"))
# t-cells
FeaturePlot(seuset.PCA, features.plot = c("CD3E"))
# b-cells
FeaturePlot(seuset.PCA, features.plot = c("CD79A"))
# mast cells
FeaturePlot(seuset.PCA, features.plot = c("KIT"))

FeaturePlot(seuset.PCA, features.plot = c("ACKR1"))

TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "Calcification", plot.title = "Calcification")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "AE", plot.title = "AE number")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "ID", plot.title = "AE and plate ID")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "Phenotype", plot.title = "Plaque phenotype")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "SR.score", plot.title = "SR score")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "CD68.score", plot.title = "CD68 score")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "CD3.score", plot.title = "CD3 score")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "CD34.score", plot.title = "CD34 score")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "alpha.SMA.score", plot.title = "alpha-SMA score")
TSNEPlot(seuset.PCA, do.label = T, pt.size = 1, label.size = 5, group.by = "Glyc.c.score", plot.title = "Glyc-c score")

dev.off()

# check PCA
PCAPlot(seuset.PCA,
        group.by = "Phenotype")
PCAPlot(seuset.PCA,
        group.by = "AE")
PCAPlot(seuset.PCA,
        group.by = "Calcification")
PCAPlot(seuset.PCA,
        group.by = "C.H")

#-----------------------------------------------------------------------------------------

saveRDS(seuset.PCA, file = "seuset 9 patients.RDS")
