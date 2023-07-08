suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(Matrix)
    library(SingleCellExperiment)
    library(igraph)
})

x <- readRDS(args$dat)
# remove outlier
qc <- perCellQCMetrics(x)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
x <- x[, !ol]

# remove lowly expressed genes
x <- x[rowSums(counts(x) > 1) >= 10, ]

# compute sum-factors & normalize
x <- computeLibraryFactors(x)
x <- logNormCounts(x)

# sample
gid <- sample(seq_len(nrow(x)), 3000)
cid <- sample(seq_len(ncol(x)), 5000)
x <- x[gid, cid]

# hvg
tbl <- modelGeneVar(x, block = x$sample_id)
rowData(x)$hvg <- hvg <- tbl$bio > 0
rowData(x)$bio <- tbl$bio 
x <- runPCA(x, subset_row = hvg, ncomponents = 10)
x <- runUMAP(x)
# table of gene/cell metadata
# and simulation parameters
dr <- reducedDim(x, "PCA")

# high- & low-resolution clustering
g <- buildSNNGraph(x, use.dimred = "PCA")
x$cluster_hi <- cluster_louvain(g, resolution = 2)$membership
x$cluster_lo <- cluster_louvain(g, resolution = 0.1)$membership

for (. in names(colData(x)))
    x[[.]] <- factor(x[[.]])
cd <- cbind(dr, colData(x))
rd <- rowData(x)


# save data
saveRDS(x, args$pro)
saveRDS(rd, args$rd)
saveRDS(cd, args$cd)

