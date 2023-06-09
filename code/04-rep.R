# args <- list(
#     "data/01-fil/t0,s0,b0.rds",
#     "outs/sel-t0,s0,b0,entropy,top40p.rds",
#     "data/02-rep/t0,s0,b0.rds")

suppressPackageStartupMessages({
    library(scran)
    library(scater)
})

sce <- readRDS(args[[1]])
sel <- readRDS(args[[2]])

sce <- runPCA(sce, subset_row = sel$sel_val)
sce$cluster_id <- clusterCells(sce, use.dimred = "PCA")

saveRDS(sce, args[[3]])
