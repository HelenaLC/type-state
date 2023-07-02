fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    f <- x[["type_Fstat"]]
    e <- x[["state_edgeR"]]
    o <- rank(f$sco_val)
    l <- rank(e$sco_val)
    s <- order(o - l, decreasing = TRUE)[seq_len(round(nrow(e)*0.25))]
    res <- rep(FALSE, nrow(f))
    res[s] <- TRUE
    return(res)
    #s <= round(nrow(e)*0.3)
    
}