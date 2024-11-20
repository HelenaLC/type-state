suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    n <- 2e3
    t <- rank((y <- x$tF)$sco_val)
    s <- rank(x$sPBDS$sco_val)
    o <- order(t-s, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}