suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    t <- rank((y <- x$type_Fstat)$sco_val)
    s <- rank(x$state_edgeR$sco_val)
    o <- order(t-s, decreasing=TRUE)
    y$gene_id[o[seq_len(1e3)]]
}