suppressPackageStartupMessages({
    library(muscat)
    library(scater)
})

fun <- \(x) {
    # compute pseudo-bulks by cluster-sample
    x <- SingleCellExperiment(assays = list(counts = assay(x, "counts")),
        colData = DataFrame(sample_id = x$sample_id,
            condition = x$group_id,
            cluster_id = x$cluster_hi))
    
    y <- aggregateAcrossCells(x, 
        colData(x)[c("cluster_id", "sample_id")], 
        coldata.merge = FALSE,
        use.assay.type = "counts", 
        statistics = "sum")
    
    # split cell indices by samples
    idx <- split(seq(ncol(y)), y$sample_id)
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
        n <- log(length(p), base = 2)
        -sum(p*log(p, base = 2))/n
    }
    # ...across clusters for every sample
    hs <- sapply(p, apply, 1, h)
    z <- aggregateAcrossCells(x, 
        x$sample_id, coldata.merge = FALSE,
        use.assay.type = "counts", statistics = "sum")
    res <- data.frame(
        row.names = NULL, entropy = c(hs),
        marker_id = rep(rownames(hs), ncol(hs)),
        sample_id = rep(colnames(hs), each = nrow(hs)))
    # when cluster size is small,
    # a gene is likely to only express in one cluster
    res[is.na(res)] <- 0
    c(by(res, res$marker_id, \(.) mean(1-.$entropy)))
}
