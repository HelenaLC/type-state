# wcs <- list(sim="t100,s0,b0", sel="DUBStepR")
# args <- list(
#     sprintf("data/01-fil/%s.rds", wcs$sim),
#     sprintf("outs/sel-sim-%s,%s.rds", wcs$sim, wcs$sel),
#     sprintf("data/02-rep/%s,%s.rds", wcs$sim, wcs$sel))

suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(SingleCellExperiment)
})

sce <- readRDS(args[[1]])
res <- readRDS(args[[2]])

idx <- match(rownames(sce), res$gene_id)
rowData(sce)$sel_val <- sel <- res$sel_val[idx]

sce <- runPCA(sce, subset_row=sel)
sce <- runUMAP(sce, dimred="PCA")
sce$cluster_re <- clusterCells(sce, use.dimred="PCA")

dr <- data.frame(as.list(reducedDims(sce)))
cd <- data.frame(row.names=NULL, wcs, dr, colData(sce))
if (!is.null(wcs$sim)) cd <- data.frame(cd, metadata(sce))

saveRDS(sce, args[[3]])
saveRDS(cd, args[[4]])