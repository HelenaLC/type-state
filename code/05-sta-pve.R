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
    cd <- cd[c("cluster_id", "sample_id", "group_id")]
    mod <- ~(1|cluster_id)+(1|sample_id)+(1|group_id)
    res <- fitExtractVarPartModel(y, mod, cd, quiet = TRUE)
    # drop residuals & re-scale proportions
    tmp <- as.matrix(res <- res[, -4])
    res[, -4] <- sweep(tmp, 1, rowSums(tmp), `/`)
    res[rowAlls(is.na(res)), ] <- 0
    # return variance fractions
    # separately for each variable
    rbind(
        data.frame(sta = "pve_k", sta_val = res$cluster_id),
        data.frame(sta = "pve_s", sta_val = res$sample_id),
        data.frame(sta = "pve_g", sta_val = res$group_id))
}