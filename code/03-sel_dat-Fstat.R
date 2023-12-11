suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    t <- rank((y <- x$Fstat)$sco_val)
    o <- order(t, decreasing=TRUE)
    y$gene_id[o[seq_len(1e3)]]
}