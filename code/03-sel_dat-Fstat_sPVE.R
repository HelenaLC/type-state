suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    n <- 2e3
    s <- rank(x$sPVE$sco_val)
    t <- rank((y <- x$Fstat)$sco_val)
    o <- order(t-s, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}