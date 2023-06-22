suppressPackageStartupMessages({
    library(scuttle)
    library(SummarizedExperiment)
})

fun <- \(x) {

    ng <- nrow(x)
    y <- logcounts(x)
    idx <- split(seq(ncol(x)), x$cluster_hi)
    res <- vapply(idx, \(i) {
        grp <- c(1, 2)[1 + (seq(ncol(x)) %in% i)]
        vapply(seq(ng), \(g) {
            wilcox.test(y[g, ] ~ grp)$p.value
        }, numeric(1))
    }, numeric(ng))
    
    res <- -log(res)
    apply(res, 1, mean)
}
