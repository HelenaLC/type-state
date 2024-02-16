suppressPackageStartupMessages({
    library(lisi)
    library(SingleCellExperiment)
})

fun <- \(x){
    y <- reducedDim(x, "PCA")
    cd <- colData(x)
    by <- c(
        LISI_g="group_id", 
        LISI_k="cluster_id")
    res <- sapply(by, \(.) {
        res <- compute_lisi(y, cd, .)
        n <- length(unique(cd[[.]]))
        #mean(res[[1]])-(n-1)
        (mean(res[[1]])-1)/(n-1)
    })
    data.frame(
        row.names=NULL, 
        sta=names(res), 
        sta_val=res)
}
