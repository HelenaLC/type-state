fun <- \(x) {
    y <- x$random
    o <- order(y$sco_val, decreasing = FALSE)[seq_len(round(nrow(y)*0.3))]
    res <- rep(FALSE, nrow(y))
    res[o] <- TRUE
    return(res)
}
