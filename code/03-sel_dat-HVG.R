fun <- \(x) {
    y <- x$hvg
    z <- y$bio
    y$gene_id[z > 0]
}