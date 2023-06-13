fun <- \(x) {
    y <- x$sco_val[x$sco == "entropy"]
    o <- order(y, decreasing = FALSE)
    o <= round(length(y)*0.25)
}
