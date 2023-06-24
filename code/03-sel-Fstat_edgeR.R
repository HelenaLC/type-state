fun <- \(x) {
    f <- x$type_Fstat
    e <- x$state_edgeR
    f$rank <- rank(f$sco_val)
    e$rank <- rank(e$sco_val)
    r <- f$rank - e$rank
    o <- order(r, decreasing = TRUE)
    return(o <= round(length(r)*0.25))
}