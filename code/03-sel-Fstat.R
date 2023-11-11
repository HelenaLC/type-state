fun <- \(x) {
    y <- x$type_Fstat
    y <- y[!is.na(y$sco_val), ]
    z <- order(y$sco_val, decreasing=TRUE)[seq_len(round(nrow(y))*0.25)]
    y$gene_id[z]
}
