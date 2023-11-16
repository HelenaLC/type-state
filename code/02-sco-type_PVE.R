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
    res[rowAlls(is.na(res)), ] <- 0
    # return variance fraction attributable to cluster
    setNames(res$cluster_hi, rownames(res))
}