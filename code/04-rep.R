# args <- list(
#     "data/01-fil/t0,s0,b0.rds",
#     "outs/sel-t0,s0,b0,Entropy_Fstat.rds",
#     "data/02-rep/0,s0,b0,Entropy_Fstat.rds")

suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(SingleCellExperiment)
    library(harmony)
    library(stringr)
})

sce <- readRDS(args[[1]])
sel <- readRDS(args[[2]])

rowData(sce)$sel_val <- sel$sel_val
rowData(sce)$sel <- sel$sel

sce <- runPCA(sce, subset_row = sel$sel_val)
sce <- runUMAP(sce, subset_row = sel$sel_val)
sce$cluster_re <- clusterCells(sce, use.dimred = "PCA")

pca <- reducedDim(sce, "PCA")
umap <- reducedDim(sce, "UMAP")
if (str_detect(args[[1]], "fil")) {
    cd <- data.frame(colData(sce), metadata(sce), wcs,
        pca, umap, row.names = NULL)
} else {
    cd <- data.frame(colData(sce), wcs,
        pca, umap, row.names = NULL)
}

#rd <- data.frame(rowData(sce), wcs)

saveRDS(sce, args[[3]])
saveRDS(cd, args[[4]])