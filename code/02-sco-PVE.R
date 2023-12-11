suppressPackageStartupMessages({
    library(variancePartition)
})

fun <- \(x) {
    y <- assay(x, "logcounts")
    cd <- data.frame(colData(x))
    cd <- cd[c("cluster_hi", "sample_id", "group_id")]
    mod <- ~(1|cluster_hi)+(1|sample_id)+(1|group_id)
    res <- fitExtractVarPartModel(y, mod, cd, quiet = TRUE)
    res[rowAlls(is.na(res)), ] <- 0
    lapply(list(
        tPVE=res$cluster_hi, 
        sPVE=res$group_id), 
        setNames, rownames(res))
}