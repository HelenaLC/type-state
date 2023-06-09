suppressPackageStartupMessages({
    library(muscat)
    library(scater)
})

fun <- \(x, 
    cluster = "cluster_id", 
    sample = "sample_id",    
    assay = "counts", 
    fun = "sum",
    base = 2) {
    # compute pseudo-bulks by cluster-sample
    y <- aggregateAcrossCells(x, 
        colData(x)[c(cluster, sample)], 
        coldata.merge = FALSE,
        use.assay.type = assay, 
        statistics = fun)
    # split cell indices by samples
    idx <- split(seq(ncol(y)), y[[sample]])
    # compute proportions across clusters
    # (matrix of dim. features x clusters)
    p <- lapply(idx, \(.) {
        z <- assay(y)[, .]
        z <- as.matrix(z)
        z[z < 0] <- 0
        prop.table(z, 1) # equivalent to x/sum(x)
    })
    # compute entropy...
    h <- \(p) {
        p <- p[p > 0]
        n <- log(length(p), base = base)
        -sum(p*log(p, base = base))/n
    }
    # ...across clusters for every sample
    hs <- sapply(p, apply, 1, h)
    z <- aggregateAcrossCells(x, 
        x[[sample]], coldata.merge = FALSE,
        use.assay.type = assay, statistics = fun)
    res <- data.frame(
        row.names = NULL, entropy = c(hs),
        marker_id = rep(rownames(hs), ncol(hs)),
        sample_id = rep(colnames(hs), each = nrow(hs)))
    # when cluster size is small,
    # a gene is likely to only express in one cluster
    res[is.na(res)] <- 0
    c(by(res, res$marker_id, \(.) mean(1-.$entropy)))
}
