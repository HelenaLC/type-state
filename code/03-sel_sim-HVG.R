fun <- \(x) {
    y <- x$HVG
    z <- y$sco_val > 0
    y$gene_id[z]
}