suppressPackageStartupMessages({
    library(variancePartition)
})

fun <- \(x) {
    # subset selected features
    y <- assay(x, "logcounts")
    cd <- data.frame(colData(x))
    cd <- cd[c("cluster_hi", "sample_id", "group_id")]
    mod <- ~(1|cluster_hi)+(1|sample_id)+(1|group_id)
    res <- fitExtractVarPartModel(y, mod, cd, quiet = TRUE)
    # drop residuals & re-scale proportions
    # tmp <- as.matrix(res <- res[, -3])
    # res[, -3] <- sweep(tmp, 1, rowSums(tmp), `/`)
    res[rowAlls(is.na(res)), ] <- 0
    # return variance fractions
    # separately for each variable
    setNames(res$group_id, rownames(x))
}