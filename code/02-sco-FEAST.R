suppressPackageStartupMessages({
    library(FEAST)
    library(SingleCellExperiment)
})

fun <- \(x) {
    y <- as.matrix(assay(x, "logcounts"))
    k <- length(unique(x$cluster_id))
    ids <- Consensus(y, k=k)$cluster
    res <- cal_F2(y, ids)$F_scores
    setNames(res, rownames(x))
}
