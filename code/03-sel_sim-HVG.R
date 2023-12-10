fun <- \(x) {
    y <- x$hvg
    z <- y$bio > 0
    y$gene_id[z]
}