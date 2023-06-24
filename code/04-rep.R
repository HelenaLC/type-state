# args <- list(
#     "data/01-fil/t0,s0,b0.rds",
#     "outs/sel-t0,s0,b0,Entropy_Fstat.rds",
#     "data/02-rep/0,s0,b0,Entropy_Fstat.rds")

suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(SingleCellExperiment)
})

sce <- readRDS(args[[1]])
sel <- readRDS(args[[2]])

rowData(sce)$sel_val <- sel$sel_val
sce <- runPCA(sce, subset_row = sel$sel_val)
sce$cluster_re <- clusterCells(sce, use.dimred = "PCA")

saveRDS(sce, args[[3]])
