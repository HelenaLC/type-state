suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    de <- grep("^GroupDE", names(y))
    de <- rowAnys(y[de] != 1)
    y$gene_id[de]
}