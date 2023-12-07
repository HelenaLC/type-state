suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    t <- rank((y <- x$type_Fstat)$sco_val)
    o <- order(t, decreasing=TRUE)
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    n <- sum(rowAnys(y[de] != 1) & rowAlls(y[ds] == 1))
    y$gene_id[o[seq_len(n)]]
}