fun <- \(x) {
    y <- x$FEAST
    z <- order(y$sco_val, decreasing=TRUE)[seq_len(round(nrow(y))*0.25)]
    y$gene_id[z]
}
