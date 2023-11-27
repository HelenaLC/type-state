suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    n <- if (!is.null(y$dat)) {
        # if ground truth unavailable, select 2,000
        2e3
    } else {
        # select number of genes that are truly DE but not DS
        de <- grep("^GroupDE", names(y))
        ds <- grep("^ConditionDE", names(y))
        sum(rowAnys(y[de] != 1) & rowAlls(y[ds] == 1))
    }
    o <- order(y$sco_val, decreasing=FALSE)
    y$gene_id[o <= n]
}
