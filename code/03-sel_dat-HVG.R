fun <- \(x) {
    z <- (y <- x$HVG)$sco_val
    o <- order(z, decreasing=TRUE)
    y$gene_id[o[seq_len(1e3)]]
}