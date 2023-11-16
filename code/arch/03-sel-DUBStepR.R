fun <- \(x) {
    y <- x$DUBStepR
    y <- y[!is.na(y$sco_val), ]
    y$gene_id[y$sel_val]
}
