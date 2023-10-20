fun <- \(x) {
    y <- x$type_Fstat
    y <- y[!is.na(y$sco_val), ]
    y$gene_id[y$sco_val >= 0.95]
}
