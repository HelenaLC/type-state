fun <- \(x) {
    y <- x$type_entropy
    o <- order(y$sco_val, decreasing = TRUE)
    return(o <= round(nrow(y)*0.3))
}
