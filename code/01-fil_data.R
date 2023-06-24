suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(Matrix)
    library(SingleCellExperiment)
    library(igraph)
})

# load data
x <- readRDS(args$sim)

# standardize cell metadata
colData(x) <- DataFrame(
    row.names = NULL,
    cluster_id = x$Group,
    sample_id = x$Sample,
    group_id = x$Condition)

# filtering;
# keep genes with count > 1 in at least 10 cells,
# and cells with at least 10 detected genes
y <- counts(x)
gs <- rowSums(y > 1) >= 10
cs <- colSums(y > 0) >= 10
x <- x[gs, cs]

# preprocessing;
# log-library size normalization, 
# feature selection, dimension reduction
x <- logNormCounts(x)
tbl <- modelGeneVar(x, block = x$sample_id)
rowData(x)$hvg <- hvg <- tbl$bio > 0
rowData(x)$bio <- tbl$bio 
#x <- runPCA(x, subset_row = hvg, ncomponents = 10)
x <- runPCA(x, subset_row = hvg, ncomponents = 10)

# table of gene/cell metadata
# and simulation parameters
dr <- reducedDim(x, "PCA")

# Clustering at high resolution
g <- buildSNNGraph(x, use.dimred = "PCA")
x$cluster_hi <- cluster_louvain(g, resolution = 2)$membership

# Clustering at low resolution
cluster_lo <- cluster_louvain(g, resolution = 0)$membership
gc <- table(cluster_lo, x$group_id)

if (sum(rowSums(gc == 0) == 1) == nrow(gc)) {
    x$cluster_lo <- 1
} else {
    x$cluster_lo <- cluster_lo
}

md <- data.frame(metadata(x), wcs)
cd <- cbind(md, dr, colData(x))
rd <- cbind(md, rowData(x), md)

# save data
saveRDS(x, args$fil)
saveRDS(rd, args$rd)
saveRDS(cd, args$cd)