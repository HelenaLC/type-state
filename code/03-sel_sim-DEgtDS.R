suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    y <- x$random
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    de <- sapply(de, \(i) {
        j <- setdiff(de, i)
        lfc <- log2(y[j]/y[[i]])
        rowMeans(lfc) }) |> rowMeans()
    ds <- sapply(ds, \(i) {
        j <- setdiff(ds, i)
        lfc <- log2(y[j]/y[[i]])
        rowMeans(lfc) }) |> rowMeans()
    y$gene_id[de > ds]
}