suppressPackageStartupMessages({
    library(MLVSBM)
    library(SummarizedExperiment)
})

fun <- \(x){
    res <- ARI(x$cluster_id, x$cluster_re)
    data.frame(sta_val = res, row.names = NULL)
}