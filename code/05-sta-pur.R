suppressPackageStartupMessages({
    library(scater)
    library(bluster)
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
    # cluster purity should be high
    # within groups (high value)
    j <- split(i, x$group_id)
    res_by_k <- vapply(j, \(.) {
        ids <- x$cluster_re[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- neighborPurity(y[., ], ids)
        mean(res$purity)
    }, numeric(1))
    # group purity should be low
    # within clusters (low value)
    j <- split(i, x$cluster_re)
    res_by_g <- vapply(j, \(.) {
        ids <- x$group_id[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- neighborPurity(y[., ], ids)
        mean(res$purity)
    }, numeric(1))
    # average across samples/clusters
    res_g <- mean(res_by_g, na.rm=TRUE) 
    res_k <- mean(res_by_k, na.rm=TRUE) 
    # return as separate statistics
    df_g <- data.frame(
        sta="pur_g",
        sta_val=res_g,
        row.names=NULL)
    df_k <- data.frame(
        sta="pur_k",
        sta_val=res_k,
        row.names=NULL)
    df <- rbind(df_g, df_k)
    return(df)
}