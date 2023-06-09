suppressPackageStartupMessages({
    library(FEAST)
    library(SummarizedExperiment)
})

fun <- \(x, 
    assay_to_use = "logcounts", 
    cluster_to_use = "cluster_id",
    sample_to_use = "sample_id") {
    
    ids <- unique(x[[sample_to_use]])
    res <- sapply(ids, \(i) {
        y <- x[, which(x[[sample_to_use]] == i)]
        mtx <- as.matrix(assay(y, assay_to_use))
        cal_F2(mtx, y[[cluster_to_use]])$F_scores
    })
    rownames(res) <- rownames(x)
    rowMeans(res, na.rm = TRUE)
}
