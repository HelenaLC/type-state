suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

fun <- \(x) {
    y <- table(x$sample_id, x$cluster_re)
    y <- t(as.matrix(unclass(y)))
    
    ids <- colnames(y)
    idx <- match(ids, x$sample_id)
    df <- data.frame(
        row.names = ids,
        sample_id = ids,
        group_id = x$group_id[idx])
    mm <- model.matrix(~ group_id, df)
    
    y <- tryCatch(
        voom(DGEList(y), mm, plot = FALSE),
        error = function(e) NULL)
    df <- if (!is.null(y)) {
        fit <- lmFit(y, mm)
        fit <- contrasts.fit(fit, c(0, 1))
        fit <- eBayes(fit, trend = FALSE)
        
        tbl <- topTable(fit, number = Inf, sort.by = "none")
        idx <- match(c("P.Value", "adj.P.Val"), names(tbl))
        names(tbl)[idx] <- c("p_val", "p_adj")
        
        cell_n <- tabulate(x$cluster_re)
        cell_i <- split(seq(ncol(x)), x$cluster_re)
        DataFrame(
            row.names = NULL, tbl,
            cluster_re = rownames(y), 
            cell_n, cell_i = I(cell_i))
    }
    return(df)
}
