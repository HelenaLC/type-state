suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$HVG
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    n <- sum(rowAnys(y[de] != 1) & rowAlls(y[ds] == 1))
    o <- order(y$sco_val, decreasing=TRUE)
    y$gene_id[o[seq_len(n)]]
}