# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(igraph)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
res <- readRDS(args[[2]])

# get selected features
idx <- match(rownames(sce), res$gene_id)
rowData(sce)$sel_val <- sel <- res$sel_val[idx]

# clustering & reduction
sce <- runPCA(sce, subset_row=sel)
g <- buildSNNGraph(sce, use.dimred="PCA", type="jaccard")
k <- cluster_louvain(g, resolution=1)$membership
sce$cluster_re <- factor(k, seq_along(unique(k)))
if (!is.null(wcs$dat)) sce <- runUMAP(sce, dimred="PCA")

# wrangling
cd <- data.frame(row.names=NULL, wcs, colData(sce))
dr <- data.frame(unname(as.list(reducedDims(sce))))
cd <- cd[!names(cd) %in% names(dr)]; cd <- cbind(cd, dr)
if (!is.null(wcs$sim)) cd <- data.frame(cd, metadata(sce))

# saving
saveRDS(sce, args[[3]])
saveRDS(cd, args[[4]])