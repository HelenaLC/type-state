fun <- \(x) {
    y <- x$type_Fstat
    o <- order(y$sco_val, decreasing = TRUE)
    return(o <= round(nrow(y)*0.25))
}