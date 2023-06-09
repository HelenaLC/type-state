suppressPackageStartupMessages({
    library(CellMixS)
    library(SingleCellExperiment)
})

fun <- \(x) {
    idx <- split(seq(ncol(x)), x$cluster_id)
    res <- vapply(idx, \(.) {
        res <- cms(x[, .], 15, "sample_id",
            n_dim = ncol(reducedDim(x, "PCA")))
        mean(res$cms)
    }, numeric(1))
    data.frame(
        sample_id = names(res),
        sta_val = res, row.names = NULL)
}
