fun <- \(x) {
    #y <- x$sco_val[x$sco == "entropy"]
    f <- x[["type_Fstat"]]
    e <- x[["state_edgeR"]]
    f$rank <- rank(f$sco_val)
    e$rank <- rank(e$sco_val)
    f$fe <- f$rank - e$rank
    
    o <- order(f$fe, decreasing = TRUE)
    o <= round(length(f)*0.25)

}