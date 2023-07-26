suppressPackageStartupMessages({
    library(DUBStepR)
})

fun <- \(x) {
    y <- assay(x, "logcounts")
    colnames(y) <- seq_len(ncol(y))
    #y <- ScaleData(y)
    dub <- DUBStepR(y, optimise.features = FALSE)
    res <- rep(0, nrow(x))
    idx <- match(dub$corr.info$feature.genes, rownames(x))
    res[idx] <- dub$corr.info$corr.range
    return(res)
}