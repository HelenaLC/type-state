fun <- \(x) {
    y <- runif(nrow(x))
    names(y) <- rownames(x)
    return(y)
}
