suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(MKmisc)
})

fun <- \(x, 
    fun = mean, 
    fun_pb = "mean",
    penalty = FALSE, 
    assay_to_use = "logcounts", 
    cluster_to_use = "cluster_id") {
    # one vs rest t test
    y <- assay(x, assay_to_use)
    ids <- unique(x[[cluster_to_use]])
    res <- sapply(ids, \(k) {
        tmp <- x
        id <- tmp[[cluster_to_use]]
        j <- !(i <- id == k)
        ij <- c(which(i), which(j))
        df <- data.frame(row.names = ij)
        df[i, "group"] <- "Group1"
        df[j, "group"] <- "Group2"
        df$group <- factor(df$group)
        z <- y[, as.numeric(rownames(df)), drop = FALSE]
        mod.t.test(as.matrix(z), group = df[, "group"])$adj.p.value
    })
    rownames(res) <- rownames(x)
    res <- apply(res, 1, FUN = fun, na.rm = TRUE)

    if (!penalty) return(-log(res))
    #x[[cluster_to_use]] <- cluster_ids(x, clustering_to_use)
    pb <- muscat::aggregateData(x,
         by = cluster_to_use,
         assay = assay_to_use,
         fun = fun_pb)
    minmax_norm <- \(x) (x - min(x)) / (max(x) - min(x))
    norm_mat <- apply(assay(pb), 2, minmax_norm)
    max_mat <- apply(norm_mat, 1, max)
    max_mat*(-log(res))
}
