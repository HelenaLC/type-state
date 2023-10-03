suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(IEntropy)
    library(rsvd)
})

fun <- \(x){
    y <- assay(x, "logcounts")
    ent <- Get_entropy(y, K = 2)
    idx <- match(rownames(x), ent$Gene)
    res <- ent$entropy_int[idx]
    res[is.na(res)] <- 0
    return(res)
}