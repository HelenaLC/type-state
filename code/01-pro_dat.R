# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(Matrix)
    library(igraph)
    library(harmony)
})

# loading
x <- readRDS(args[[1]])

# remove lowly expressed genes
x <- x[rowSums(counts(x) > 1) >= 20, ]

# normalization
x <- computeLibraryFactors(x)
x <- logNormCounts(x)

# selection
tbl <- modelGeneVar(x, block=x$sample_id)
table(hvg <- (rowData(x)$bio <- tbl$bio) > 0)

# reduction
x <- runPCA(x, subset_row=hvg)
x <- runUMAP(x, dimred="PCA")

# integration
y <- HarmonyMatrix(
    data_mat=reducedDim(x, "PCA"),
    meta_data=colData(x),
    vars_use="sample_id")
names(y) <- paste0("PC", names(y))
reducedDim(x, "HARMONY") <- y

# high- & low-resolution clustering
g <- buildSNNGraph(x, use.dimred = "PCA")
x$cluster_hi <- cluster_louvain(g, resolution=2)$membership
x$cluster_lo <- cluster_louvain(g, resolution=0.1)$membership

# table of gene/cell metadata
# and simulation parameters
dr <- reducedDim(x, "PCA")
for (. in names(colData(x)))
    x[[.]] <- factor(x[[.]])
rd <- data.frame(wcs, rowData(x))
cd <- data.frame(wcs, colData(x), dr)

# saving
saveRDS(x, args[[2]])
saveRDS(rd, args[[3]])
saveRDS(cd, args[[4]])
