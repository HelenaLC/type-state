fun <- \(x) {
    y <- x$state_limma
    o <- order(y$sco_val, decreasing = FALSE)
    o <= round(nrow(y)*0.25)
}
