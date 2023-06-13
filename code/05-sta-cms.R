suppressPackageStartupMessages({
    library(CellMixS)
    library(SingleCellExperiment)
})

fun <- \(x) {
    # clusters should be well separated
    # within samples (low value)
    knn <- max(table(x$sample_id))/5
    idx <- split(seq(ncol(x)), x$sample_id)
    res_by_s <- vapply(idx, \(.) {
        res <- cms(x[, .], knn, "cluster_id",
            n_dim = ncol(reducedDim(x, "PCA")))
        mean(res$cms)
    }, numeric(1))
    # samples should be well mixed 
    # within clusters (high value)
    knn <- max(table(x$cluster_id))/5
    idx <- split(seq(ncol(x)), x$cluster_id)
    res_by_k <- vapply(idx, \(.) {
        res <- cms(x[, .], knn, "sample_id",
            n_dim = ncol(reducedDim(x, "PCA")))
        mean(res$cms)
    }, numeric(1))
    # average across samples/clusters & groupings
    res_by_s <- mean(res_by_s, na.rm = TRUE)
    res_by_k <- mean(res_by_k, na.rm = TRUE)
    res <- ( (1 - res_by_s) + res_by_k ) / 2
    data.frame(sta_val = res, row.names = NULL)
}
