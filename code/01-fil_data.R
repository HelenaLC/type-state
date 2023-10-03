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

# pca <- HarmonyMatrix(data_mat = assay(x, "counts"), 
#     meta_data = colData(x), 
#     vars_use = "sample_id")
# reducedDim(x, "PCA") <- pca


# table of gene/cell metadata
# and simulation parameters
dr <- reducedDim(x, "PCA")

# high- & low-resolution clustering
g <- buildSNNGraph(x, use.dimred = "PCA")
x$cluster_hi <- cluster_louvain(g, resolution = 2)$membership
x$cluster_lo <- cluster_louvain(g, resolution = 0)$membership
# x$cluster_hi <- leiden(g, resolution_parameter = 2)
# x$cluster_lo <- leiden(g, resolution_parameter = 0.1)


    
# mock cluster identifier if low-resolution
# clusters aren't represented in both groups
ns <- table(x$cluster_lo, x$group_id)
na <- rowSums(ns == 0)
if (sum(na == 1) == nrow(ns)) 
    x$cluster_lo <- 1

# factorize cluster identifiers
x$cluster_lo <- factor(x$cluster_lo, sort(unique(x$cluster_lo)))
x$cluster_hi <- factor(x$cluster_hi, sort(unique(x$cluster_hi)))

md <- data.frame(metadata(x), wcs)
cd <- cbind(md, dr, colData(x))
rd <- cbind(md, rowData(x), md)

# save data
saveRDS(x, args$fil)
saveRDS(rd, args$rd)
saveRDS(cd, args$cd)