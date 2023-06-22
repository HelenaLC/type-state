fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    y <- x[["type_Fstat"]]
    o <- order(y$sco_val, decreasing = TRUE)
    o <= round(length(y)*0.25)
}