suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(Matrix)
    library(SingleCellExperiment)
    library(igraph)
    library(leiden)
    library(harmony)
})

# load data
x <- readRDS(args$sim)

# filtering:
# keep genes with count > 1 in at least 10 cells,
# and cells with at least 10 detected genes
y <- counts(x)
gs <- rowSums(y > 1) >= 10
cs <- colSums(y > 0) >= 10
x <- x[gs, cs]

# preprocessing:
# log-library size normalization, 
# feature selection, dimension reduction
x <- logNormCounts(x)
tbl <- modelGeneVar(x, block=x$sample_id)
rowData(x)$hvg <- hvg <- tbl$bio > 0
rowData(x)$bio <- tbl$bio 
x <- runPCA(x, subset_row=hvg, ncomponents=10)

# high- & low-resolution clustering
g <- buildSNNGraph(x, use.dimred="PCA")
x$cluster_hi <- cluster_louvain(g, resolution=2)$membership
x$cluster_lo <- cluster_louvain(g, resolution=0)$membership

# set cluster identifier if low-resolution
# clusters aren't represented in both groups
ns <- table(x$cluster_lo, x$group_id)
na <- rowSums(ns == 0)
if (sum(na == 1) == nrow(ns)) 
    x$cluster_lo <- 1

# factorize cluster identifiers
ks_lo <- sort(unique(x$cluster_lo))
ks_hi <- sort(unique(x$cluster_hi))
x$cluster_lo <- factor(x$cluster_lo, ks_lo)
x$cluster_hi <- factor(x$cluster_hi, ks_hi)

# table of gene/cell metadata
# and simulation parameters
dr <- reducedDim(x, "PCA")
md <- data.frame(metadata(x), wcs)
rd <- cbind(md, rowData(x))
cd <- cbind(md, colData(x), dr)

# save data
saveRDS(x, args$fil)
saveRDS(rd, args$rd)
saveRDS(cd, args$cd)