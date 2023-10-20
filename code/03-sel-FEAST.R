fun <- \(x) {
    y <- x$FEAST
    z <- order(y$sco_val, decreasing=TRUE)
    y$gene_id[z <= nrow(y)/4]
}
