fun <- \(x) {
    y <- x$random
    o <- order(y$sco_val, decreasing = FALSE)
    o <= round(nrow(y)*0.25)
}