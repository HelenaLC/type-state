suppressPackageStartupMessages({
    library(PCAtools)
    library(CellMixS)
    library(SingleCellExperiment)
})

fun <- \(x) {
    # use 'elbow' method to determine number of PCs
    # to consider in cell-cell distance computation
    i <- seq_len(ncol(x))
    y <- reducedDim(x, "PCA")
    pve <- attr(y, "varExplained")
    npc <- findElbowPoint(pve)
    # utility function to compute
    # mixing across 'id's for each 'by'
    f <- \(x, id, by) {
        vapply(split(i, x[[by]]), \(.) {
            res <- cms(x[, .], k=30, id, n_dim=npc)$cms
            ks.test(res[!duplicated(res)], "punif")$p.value
        }, numeric(1))
    }
    # clusters should be well separated
    # within samples (CMS flat > low KS)
    res_by_k <- f(x, "cluster_re", "sample_id")
    # samples should be well mixed
    # within clusters (CMS skewed > high KS)
    res_by_s <- f(x, "sample_id", "cluster_re")
    # average across samples/clusters
    res_s <- mean(res_by_s, na.rm=TRUE)
    res_k <- mean(res_by_k, na.rm=TRUE)
    # return as separate statistics
    df_s <- data.frame(
        sta="cms_s",
        sta_val=res_s,
        row.names=NULL)
    df_k <- data.frame(
        sta="cms_k",
        sta_val=res_k,
        row.names=NULL)
    df <- rbind(df_s, df_k)
    return(df)
}
