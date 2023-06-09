suppressPackageStartupMessages({
    library(SummarizedExperiment)
})

fun <- \(x, ...,
    fun = mean, 
    fun_pb = "mean",
    penalty = FALSE, 
    assay_to_use = "logcounts", 
    cluster_to_use = "cluster_id") {
    # one vs rest t test
    ng <- seq_len(nrow(x))
    y <- assay(x, assay_to_use)
    ids <- unique(x[[cluster_to_use]])
    res <- sapply(ids, \(k) {
        tmp <- x
        ids <- tmp[[cluster_to_use]]
        j <- !(i <- ids == k)
        ij <- c(which(i), which(j))
        df <- data.frame(row.names = ij)
        df[i, "group"] <- "Group1"
        df[j, "group"] <- "Group2"
        df$group <- factor(df$group)
        z <- y[, as.numeric(rownames(df)), drop = FALSE]
        sapply(ng, function(.) {
            if (length(unique(df$group)) == 2) {
                t.test(z[., ] ~ df$group)$p.value
            } else {
                t.test(z[., ] ~ 1)$p.value
            }
        })
    })
    rownames(res) <- rownames(x)
    res <- apply(res, 1, FUN = fun, na.rm=TRUE)
    
    if (!penalty) return(-log(res))
    x[[cluster_to_use]] <- cluster_ids(x, clustering_to_use)
    pb <- muscat::aggregateData(x, 
        by = cluster_to_use, 
        assay = assay_to_use, 
        fun = fun_pb)
    minmax_norm <- \(x) (x - min(x)) / (max(x) - min(x))
    norm_mat <- apply(assay(pb), 2, minmax_norm)
    max_mat <- apply(norm_mat, 1, max)
    max_mat*(-log(res))
}
