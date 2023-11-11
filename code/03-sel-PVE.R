fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    f <- x[["type_PVE"]]
    e <- x[["state_PVE"]]
    o <- rank(f$sco_val)
    l <- rank(e$sco_val)
    s <- order(o - l, decreasing = TRUE)[seq_len(round(nrow(e)*0.25))]
    f$gene_id[s]
    #s <= round(nrow(e)*0.3)
    
}