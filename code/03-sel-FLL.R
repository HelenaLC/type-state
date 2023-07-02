fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    f <- x[["type_Fstat"]]
    e <- x[["type_logFC"]]
    a <- x[["state_limma"]]
    o <- rank(f$sco_val)
    l <- rank(e$sco_val)
    m <- rank(a$sco_val)
    s <- order(o + l - m, decreasing = TRUE)[seq_len(round(nrow(e)*0.25))]
    res <- rep(FALSE, nrow(f))
    res[s] <- TRUE
    return(res)
    #s <= round(nrow(e)*0.3)

}