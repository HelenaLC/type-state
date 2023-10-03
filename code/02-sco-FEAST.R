suppressPackageStartupMessages({
    library(FEAST)
})

fun <- \(x) {
    y <- assay(x, "counts")
    y <- as.matrix(y)
    con <- Consensus(Y, k = length(unique(x$cluster_id)))
    res <- cal_F2(y, con$cluster)
    res$F_scores
}