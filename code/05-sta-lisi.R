suppressPackageStartupMessages({
    library(lisi)
    library(SingleCellExperiment)
})

fun <- \(x){
    X <- reducedDim(x, "PCA")
    meta <- colData(x)
    lisi_g <- compute_lisi(X, meta, "group_id")
    lisi_k <- compute_lisi(X, meta, "cluster_id")
    rbind(
        data.frame(sta="lisi_g", sta_val=mean(lisi_g$group_id) - 1),
        data.frame(sta="lisi_k", sta_val=mean(lisi_k$cluster_id) - 2))
}
