fun <- \(x) {
    df <- x[[1]]
    gs <- seq_len(nrow(df))
    ng <- round(0.1*length(gs))
    sel <- sample(gs, ng)
    gs %in% sel
}