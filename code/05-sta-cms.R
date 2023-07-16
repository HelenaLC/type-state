suppressPackageStartupMessages({
    library(CellMixS)
    library(SingleCellExperiment)
})

fun <- \(x) {
    npc <- ncol(reducedDim(x, "PCA"))
    # clusters should be well separated
    # within samples (low value)
    knn <- max(table(x$sample_id))/5
    idx <- split(seq(ncol(x)), x$sample_id)
    res_by_k <- vapply(idx, \(.) {
        res <- cms(x[, .], knn, "cluster_re", n_dim = npc)
        mean(res$cms) # average across cells
    }, numeric(1))
    # samples should be well mixed 
    # within clusters (high value)
    knn <- max(table(x$cluster_re))/5
    idx <- split(seq(ncol(x)), x$cluster_re)
    res_by_s <- vapply(idx, \(.) {
        res <- cms(x[, .], knn, "sample_id", n_dim = npc)
        mean(res$cms) # average across cells
    }, numeric(1))
    # average across samples/clusters
    res_s <- mean(res_by_s, na.rm = TRUE)
    res_k <- mean(res_by_k, na.rm = TRUE)
    # return as separate statistics
    df_s <- data.frame(
        sta = "cms_s",
        sta_val = res_s,
        row.names = NULL)
    df_k <- data.frame(
        sta = "cms_k",
        sta_val = res_k,
        row.names = NULL)
    df <- rbind(df_s, df_k)
    return(df)
}
