suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(SingleCellExperiment)
    library(bluster)
})

sce <- readRDS(args[[1]])
res <- readRDS(args[[2]])

idx <- match(rownames(sce), res$gene_id)
rowData(sce)$sel_val <- sel <- res$sel_val[idx]

sce <- runPCA(sce, subset_row=sel)
if (!is.null(wcs$dat)) {
    sce$cluster_re <- clusterCells(sce,
        use.dimred = "PCA",
        BLUSPARAM = NNGraphParam(cluster.fun = "louvain",
            cluster.args = list(resolution = 0.4)))
    sce <- runUMAP(sce, dimred="PCA")
} else {
    sce$cluster_re <- clusterCells(sce, use.dimred="PCA")
}

dr <- data.frame(unname(as.list(reducedDims(sce))))
cd <- data.frame(row.names=NULL, wcs, dr, colData(sce))
if (!is.null(wcs$sim)) cd <- data.frame(cd, metadata(sce))

saveRDS(sce, args[[3]])
saveRDS(cd, args[[4]])