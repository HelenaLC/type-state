suppressPackageStartupMessages({
    library(scater)
    library(bluster)
    library(SingleCellExperiment)
})

fun <- \(x) {
    mtx <- reducedDim(x, "PCA")
    idx <- split(seq(ncol(x)), x$sample_id)
    res <- vapply(idx, \(.) {
        ids <- x$cluster_id[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- neighborPurity(mtx[., ], ids)
        mean(res$purity)
    }, numeric(1))
    data.frame(
        sample_id = names(res),
        sta_val = res, row.names = NULL)
}