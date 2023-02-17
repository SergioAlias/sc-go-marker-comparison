#!/usr/bin/env Rscript

# Sergio Alías, 20230116
# Last modified 20230120

###########################
###   scRare            ###
###   preprocessing.R   ###
###########################

# Basic preprocessing based on OSCA Basics book -> http://bioconductor.org/books/3.16/OSCA.basic/


### Imports ###

library(SingleCellExperiment)

###############


## For data from GEO:

sample_id <- "GSM3535279"
dir <- paste0("/home/user/", sample_id) # change it according to your OS and file structure
sce <- DropletUtils::read10xCounts(dir)

# if your data comes from another source - convert it to a SingleCellExperiment object

name <- "axln4"


## QC

sce <- scuttle::addPerCellQCMetrics(sce) # adds QC to colData
reasons <- scuttle::perCellQCFilters(colData(sce))
# colSums(as.matrix(reasons)) # for visualization
sce <- sce[, !reasons$discard]


## Normalization

set.seed(100)
clust <- scran::quickCluster(sce)
# table(clust) # for visualization
# deconv.sf <- scran::calculateSumFactors(sce, cluster = clust) # for visualization
# summary(deconv.sf) # for visualization
sce <- scran::computeSumFactors(sce, cluster = clust)
sce <- scuttle::logNormCounts(sce)


## Feature selection

dec <- scran::modelGeneVar(sce)
sce.original <- sce
chosen <- scran::getTopHVGs(dec, prop=0.1)
sce <- sce[chosen,]
altExp(sce, "original") <- sce.original
rm(sce.original)


## Dimensionality reduction (PCA)

set.seed(100)
sce <- scran::fixedPCA(sce, subset.row = NULL)


## Clustering

nn.clusters <- scran::clusterCells(sce, use.dimred = "PCA")
colLabels(sce) <- nn.clusters


## Saving sce as RDS

saveRDS(sce, paste0(name, ".sce.RDS"))
