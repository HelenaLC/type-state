suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

# read in raw data
x <- readRDS(args[[1]])

by <- c("cluster_id", "sample_id", "condition")
by <- intersect(by, names(colData(x)))
# keep genes with count > 1 in at least 10 cells,
# and cells with at least 10 detected genes
x <- x[
  rowSums(counts(x) > 1) >= 10,
  colSums(counts(x) > 0) >= 10]

# drop missing factor levels
for (. in by) x[[.]] <- droplevels(factor(x[[.]]))

saveRDS(x, args[[2]])