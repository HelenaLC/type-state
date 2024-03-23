suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    n <- 2e3
    y <- x$random
    o <- order(y$sco_val, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}
