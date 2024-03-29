suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    de <- rowAnys(y[de] != 1)
    ds <- rowAnys(y[ds] != 1)
    y$gene_id[de & !ds]
}