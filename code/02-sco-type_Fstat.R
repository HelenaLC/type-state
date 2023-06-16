suppressPackageStartupMessages({
    library(FEAST)
    library(SummarizedExperiment)
})

fun <- \(x) {
    
    ids <- unique(x$sample_id)
    res <- sapply(ids, \(i) {
        y <- x[, which(x$sample_id == i)]
        mtx <- as.matrix(assay(y, "logcounts"))
        cal_F2(mtx, y$cluster_hi)$F_scores
    })
    rownames(res) <- rownames(x)
    rowMeans(res, na.rm = TRUE)
}
