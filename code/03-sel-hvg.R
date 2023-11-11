fun <- \(x) {
    y <- x$hvg
    z <- y$sco_val > 0
    y$gene_id[z]
}