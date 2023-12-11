suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    s <- rank(x$sPVE$sco_val)
    t <- rank((y <- x$tPVE)$sco_val)
    o <- order(t-s, decreasing=TRUE)
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    n <- sum(rowAnys(y[de] != 1) & rowAlls(y[ds] == 1))
    y$gene_id[o[seq_len(n)]]
}