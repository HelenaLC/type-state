suppressPackageStartupMessages({
    library('variancePartition')
    library(SingleCellExperiment)
    library(dplyr)
})


fun <- \(x){
    info <- data.frame(colData(x)) %>% 
        select(cluster_re, sample_id, group_id)
    form <- ~ (1|cluster_re) + (1|sample_id) + (1|group_id)
    res <- fitExtractVarPartModel(assay(x, "logcounts"),
        form, info)
    val <- (colMeans(res)[1] - colMeans(res)[2])/(colMeans(res)[1]+colMeans(res)[2])
    data.frame(sta_val = val, row.names = NULL)
}