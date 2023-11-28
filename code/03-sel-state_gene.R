suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    # select genes that are truly DE but not DS
    y <- x$random
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    y$gene_id[rowAnys(y[ds] != 1)]
}