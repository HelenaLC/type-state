fun <- \(x) {
    t <- x$type_Fstat
    s <- x$state_edgeR
    z <- rank(t$sco_val)-rank(s$sco_val)
    idx <- order(z, decreasing=TRUE)[seq_len(round(nrow(t))*0.25)]
    t$gene_id[idx]
}