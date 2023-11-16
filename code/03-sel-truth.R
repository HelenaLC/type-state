suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    # select genes that are truly DE but not DS
    y <- x$random
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    y$gene_id[rowAnys(y[de] != 1) & rowAlls(y[ds] == 1)]
}