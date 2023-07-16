suppressPackageStartupMessages({
    library(edgeR)
    library(scuttle)
    library(SingleCellExperiment)
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
    dge <- DGEList(y, samples = df)
    dge <- estimateDisp(dge, mm, trend = "none")
    tbl <- tryCatch({
        fit <- glmQLFit(dge, mm, robust = TRUE, abundance.trend = FALSE)
        glmQLFTest(fit, coef = ncol(mm))},
        error = function(e) NULL)
    res <- if (!is.null(tbl)) {
        tbl <- data.frame(topTags(tbl, sort.by = "none"))
        idx <- match(c("PValue", "FDR"), colnames(tbl))
        colnames(tbl)[idx] <- c("p_val", "p_adj")
        cell_n <- tabulate(x$cluster_re)
        cell_i <- split(seq(ncol(x)), x$cluster_re)
        DataFrame(
            row.names = NULL, tbl,
            cluster_re = rownames(tbl), 
            cell_n[rownames(tbl)], 
            cell_i = I(cell_i)[rownames(tbl)])
    }
    
    return(res)
}