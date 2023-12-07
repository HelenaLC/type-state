fun <- \(x) {
    y <- x$HVG
    z <- y$sco_val
    y$gene_id[z > 0]
}