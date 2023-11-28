suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    # rank by type-state score
    t <- rank((y <- x$type_Fstat)$sco_val)
    s <- rank(x$state_PVE$sco_val)
    o <- order(t-s, decreasing=TRUE)
    n <- if (!is.null(y$dat)) {
        # if ground truth unavailable, select 2,000
        2e3
    } else {
        # select number of genes that are truly DE but not DS
        de <- grep("^GroupDE", names(y))
        ds <- grep("^ConditionDE", names(y))
        sum(rowAnys(y[de] != 1) & rowAlls(y[ds] == 1))
    }
    y$gene_id[o[seq_len(n)]]
}