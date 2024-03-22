fun <- \(x) {
    n <- 2e3
    z <- (y <- x$HVG)$sco_val
    o <- order(z, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}