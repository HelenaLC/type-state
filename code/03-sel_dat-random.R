suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    o <- order(y$sco_val, decreasing=TRUE)
    y$gene_id[o[seq_len(1e3)]]
}
