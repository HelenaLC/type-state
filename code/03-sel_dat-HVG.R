fun <- \(x) {
    #n <- sum(sco$ratio$sco_val>1)
    n <- 2e3
    z <- (y <- x$HVG)$sco_val
    o <- order(z, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}