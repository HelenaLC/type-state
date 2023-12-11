suppressPackageStartupMessages({
    library(matrixStats)
    library(variancePartition)
    library(SingleCellExperiment)
})

fun <- \(x) {
    # subset selected features
    y <- assay(x, "logcounts")
    y <- y[rowData(x)$sel_val, ]
    # fit LMM to estimate fraction of variance
    # attributable to cluster, sample, group
    cd <- data.frame(colData(x))
    mod <- ~(1|group_id)+(1|cluster_id)
    res <- fitExtractVarPartModel(y, mod, cd, quiet=TRUE)
    # drop residuals & re-scale proportions
    tmp <- as.matrix(res <- res[, -ncol(res)])
    res <- sweep(tmp, 1, rowSums(tmp), `/`)
    res[rowAlls(is.na(res)), ] <- 0
    res <- data.frame(res)
    # return variance fractions
    # separately for each variable
    rbind(
        data.frame(sta="PVE_g", sta_val=res$group_id),
        data.frame(sta="PVE_k", sta_val=res$cluster_id))
}