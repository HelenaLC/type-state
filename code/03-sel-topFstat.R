fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    y <- x[["type_Fstat"]]
    o <- order(y$sco_val, decreasing = TRUE)
    o <= round(nrow(y)*0.25)
}