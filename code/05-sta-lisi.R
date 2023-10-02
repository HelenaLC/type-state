suppressPackageStartupMessages({
    library(lisi)
    library(SingleCellExperiment)
})

fun <- \(x){
    X <- reducedDim(x, "PCA")
    meta <- colData(x)
    lisi <- compute_lisi(X, meta, "group_id")
    data.frame(sta_val = mean(lisi$group_id) - 1, row.names = NULL)
}
