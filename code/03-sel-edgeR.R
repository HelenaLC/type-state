fun <- \(x) {
    y <- x$state_edgeR
    o <- order(y$sco_val, decreasing = TRUE)[seq_len(round(nrow(y)*0.25))]
    res <- rep(FALSE, nrow(y))
    res[o] <- TRUE
    return(res)
}
