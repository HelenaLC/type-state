suppressPackageStartupMessages({
    library(singleCellHaystack)
    library(SummarizedExperiment)
})


fun <- \(x){
    res <- haystack(x = reducedDim(x, "PCA"), 
        expression = assay(x, "logcounts"))
    res$results$D_KL
}