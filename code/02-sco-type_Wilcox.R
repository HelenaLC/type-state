suppressPackageStartupMessages({
    library(scuttle)
    library(SummarizedExperiment)
})

fun <- \(x, ...,
    fun = "mean",
    penalty = FALSE,
    assay_to_use = "logcounts",
    cluster_to_use = "cluster_id") {

    ng <- nrow(x)
    y <- logcounts(x)
    idx <- split(seq(ncol(x)), x$cluster_id)
    res <- vapply(idx, \(i) {
        grp <- c(1, 2)[1 + (seq(ncol(x)) %in% i)]
        vapply(seq(ng), \(g) {
            wilcox.test(y[g, ] ~ grp)$p.value
        }, numeric(1))
    }, numeric(ng))
    
    res <- -log(rowMeans(res))
    if (!penalty) return(res)
    
    pbs <- aggregateAcrossCells(x, 
        ids = x$cluster_id, 
        statistics = "mean", 
        use.assay.type = "logcounts")
    
    .f <- \(.) (.-min(.)) / diff(range(.))
    pbs <- apply(assay(pbs), 2, .f)
    return(rowMaxs(pbs)*res)
}
