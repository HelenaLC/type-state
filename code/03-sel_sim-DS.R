suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    ds <- grep("^ConditionDE", names(y))
    ds <- rowAnys(y[ds] != 1)
    y$gene_id[ds]
}