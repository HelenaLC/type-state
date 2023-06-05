suppressPackageStartupMessages({
    library(Matrix)
    library(scater)
    library(scran)
    library(SingleCellExperiment)
})
# read in raw data
x <- readRDS(args[[1]])

colData(x) <- DataFrame(
    sample_id = x$Sample,
    cluster_id = x$Group,
    condition = x$Condition,
    row.names = NULL)

by <- c("cluster_id", "sample_id", "condition")
by <- intersect(by, names(colData(x)))
# keep genes with count > 1 in at least 10 cells,
# and cells with at least 10 detected genes
x <- x[
  rowSums(counts(x) > 1) >= 10,
  colSums(counts(x) > 0) >= 10]

# normalization
x <- logNormCounts(x)

# feature selection
tbl <- modelGeneVar(x, block = x$sample_id)
hvg <- tbl$bio > 0

# dimension reduction
x <- runPCA(x, ncomponents = 10, subset_row = hvg)
#x$cluster_id <- quickCluster(x, method = "igraph", graph.fun = "louvain", k = 5)

# drop missing factor levels
for (. in by) x[[.]] <- droplevels(factor(x[[.]]))

saveRDS(x, args[[2]])