suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    n <- 2e3
    t <- rank((y <- x$tF)$sco_val)
    o <- order(t, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}