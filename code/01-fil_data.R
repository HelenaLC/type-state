suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(Matrix)
    library(SingleCellExperiment)
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
x <- runPCA(x, subset_row = hvg, ncomponents = 10)

# table of gene/cell metadata
# and simulation parameters
dr <- reducedDim(x, "PCA")
md <- data.frame(metadata(x), wcs)
cd <- cbind(md, dr, colData(x))
rd <- cbind(md, rowData(x), md)

# save data
saveRDS(x, args$fil)
saveRDS(rd, args$rd)
saveRDS(cd, args$cd)