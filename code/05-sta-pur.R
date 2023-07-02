suppressPackageStartupMessages({
    library(scater)
    library(bluster)
    library(SingleCellExperiment)
})

fun <- \(x) {
    y <- reducedDim(x, "PCA")
    # cluster purity should be high
    # within samples (high value)
    idx <- split(seq(ncol(x)), x$sample_id)
    res_by_s <- vapply(idx, \(.) {
        ids <- x$cluster_re[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- neighborPurity(y[., ], ids)
        mean(res$purity)
    }, numeric(1))
    # sample purity should be low
    # within clusters (low value)
    idx <- split(seq(ncol(x)), x$cluster_re)
    res_by_k <- vapply(idx, \(.) {
        ids <- x$sample_id[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- neighborPurity(y[., ], ids)
        mean(res$purity)
    }, numeric(1))
    # average across samples/clusters
    res_s <- mean(res_by_s, na.rm = TRUE)
    res_k <- mean(res_by_k, na.rm = TRUE)
    # return as separate statistics
    df_s <- data.frame(
        sta = "pur_s",
        sta_val = res_s,
        row.names = NULL)
    df_k <- data.frame(
        sta = "pur_k",
        sta_val = res_k,
        row.names = NULL)
    df <- rbind(df_s, df_k)
    return(df)
}