fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    f <- x$type_Fstat
    e <- x$type_entropy
    o <- rank(f$sco_val)
    l <- rank(e$sco_val)
    s <- order(o + l, decreasing = FALSE)
    s <= round(nrow(f)*0.25)
}