suppressPackageStartupMessages({
    library(scater)
    library(bluster)
    library(SingleCellExperiment)
})

fun <- \(x) {
    y <- reducedDim(x, "PCA")
    # clusters should be well separated
    # within samples (high value)
    idx <- split(seq(ncol(x)), x$sample_id)
    res_by_s <- vapply(idx, \(.) {
        ids <- x$cluster_re[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- approxSilhouette(y[., ], ids)
        mean(res$width)
    }, numeric(1))
    # samples should be well mixed
    # within clusters (low value)
    idx <- split(seq(ncol(x)), x$cluster_re)
    res_by_k <- vapply(idx, \(.) {
        ids <- x$sample_id[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- approxSilhouette(y[., ], ids)
        mean(res$width)
    }, numeric(1))
    # average across samples/clusters & groupings
    res_by_s <- mean(res_by_s, na.rm = TRUE)
    res_by_k <- mean(res_by_k, na.rm = TRUE)
    res <- ( res_by_s + (1 - res_by_k) ) / 2
    data.frame(sta_val = res, row.names = NULL)
}
