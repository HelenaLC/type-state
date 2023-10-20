fun <- \(x) {
    t <- x$type_Fstat
    s <- x$state_edgeR
    z <- rank(t$sco_val)-rank(s$sco_val)
    z <- order(z, decreasing=TRUE)
    t$gene_id[z <= nrow(t)/4]
}