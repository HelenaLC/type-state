fun <- \(x) {
    f <- x$type_Fstat
    l <- x$state_limma
    #e <- x$type_entropy
    #t <- x$type_LR
    #c <- x$type_logFC
    f$rank <- rank(f$sco_val)
    #e$rank <- rank(e$sco_val)
    l$rank <- rank(l$sco_val)
    #t$rank <- rank(t$sco_val)
    #c$rank <- rank(c$sco_val)
    r <- f$rank - l$rank
    idx <- order(r, decreasing = TRUE)[seq_len(round(nrow(f)*0.25))]
    res <- rep(FALSE, nrow(f))
    res[idx] <- TRUE
    return(res)
}