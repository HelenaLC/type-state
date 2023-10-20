fun <- \(x) {
    y <- x$random
    z <- order(y$sco_val, decreasing=FALSE)
    y$gene_id[z <= nrow(y)/4]
}
