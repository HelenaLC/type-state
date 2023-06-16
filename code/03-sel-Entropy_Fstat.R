


fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    f <- x[["type_Fstat"]]
    e <- x[["type_entropy"]]
    o <- order(f$sco_val, decreasing = FALSE)
    o <= round(length(f)*0.5)
    
    l <- order(e$sco_val, decreasing = FALSE)
    l <= round(length(e)*0.5)
    intersect(o, l)
}