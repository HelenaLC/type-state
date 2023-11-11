suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    z <- y[grep("GroupDE", names(y))]
    z <- !rowAlls(as.matrix(z) == 1)
    y$gene_id[z]
}