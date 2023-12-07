suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    o <- order(y$sco_val, decreasing=TRUE)
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    n <- sum(rowAnys(y[de] != 1) & rowAlls(y[ds] == 1))
    y$gene_id[o[seq_len(n)]]
}
