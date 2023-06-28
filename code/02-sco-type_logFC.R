suppressPackageStartupMessages({
    library(scran)
})

fun <- \(x){
    pb <- aggregateAcrossCells(x,
        ids = x$cluster_hi,
        use.assay.type = "logcounts", statistics = "mean")
    y <- assay(pb, "logcounts")
    lfc <- sapply(seq_len(ncol(y)), \(i){
        not_i <- setdiff(seq_len(ncol(y)), i)
        mgk <- sapply(not_i, \(j) {
            log(y[,i]/y[,j], base = 2)
        })
        rowMeans(mgk)
    })
    apply(abs(lfc), 1, max)
}
