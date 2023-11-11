suppressPackageStartupMessages({
    library(cluster)
    library(PCAtools)
    library(SingleCellExperiment)
})

fun <- \(x) {
    # use 'elbow' method to determine number of PCs
    # to consider in cell-cell distance computation
    i <- seq_len(ncol(x))
    y <- reducedDim(x, "PCA")
    pve <- attr(y, "varExplained")
    npc <- findElbowPoint(pve)
    y <- y[, seq_len(npc), drop=FALSE]
    # utility function to compute average
    # silhouette across 'id's for each 'by'
    f <- \(x, id, by) {
        vapply(split(i, x[[by]]), \(.) {
            ids <- as.integer(factor(x[[id]][.]))
            if (length(unique(ids)) == 1) return(NA)
            res <- silhouette(ids, dist(y[., ]))
            # average across cells
            mean(res[, "sil_width"]) 
        }, numeric(1))
    }
    # groups should be mixed
    # within clusters (low value)
    res_by_g <- f(x, "group_id", "cluster_re")
    # clusters should be separated
    # within groups (high value)
    res_by_k <- f(x, "cluster_re", "group_id")
    # average across samples/clusters
    res_g <- mean(res_by_g, na.rm=TRUE)
    res_k <- mean(res_by_k, na.rm=TRUE)
    # return as separate statistics
    df_g <- data.frame(
        sta="sil_g",
        sta_val=res_g,
        row.names=NULL)
    df_k <- data.frame(
        sta="sil_k",
        sta_val=res_k,
        row.names=NULL)
    df <- rbind(df_g, df_k)
    return(df)
}
