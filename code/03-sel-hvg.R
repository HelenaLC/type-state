fun <- \(x) {

    y <- x[[1]]
    o <- order(y$bio, decreasing = TRUE)[seq_len(round(nrow(y)*0.25))]
    res <- rep(FALSE, nrow(y))
    res[o] <- TRUE
    return(res)
    # idx <- order(x[[1]]$bio, decreasing = FALSE)[seq_len(nrow(x)*0.25)]
    # idx <- match(hvg, rownames(sce))
    # rowData(sce)$hvg_sel <- FALSE
    # rowData(sce)$hvg_sel[idx] <- TRUE
    #return(x[[1]]$hvg)

}