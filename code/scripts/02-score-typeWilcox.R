suppressPackageStartupMessages(
    {
        source("code/scripts/utils.R")
        library(SummarizedExperiment)
    }
)


fun <- \(x, 
         sample = "sample_id", 
         assay_to_use = "logcounts", 
        ...){
    
    #res <- sapply(unique(colData(x)[, sample]), \(i){
    #    temp <- x[, which(colData(x)[,sample] == i)]
    res <- .one_vs_rest_wilcox_score(sce, assay_to_use = assay_to_use, ...)
    #})
    res
}

