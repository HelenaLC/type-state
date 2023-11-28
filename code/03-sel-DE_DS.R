suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    # select genes that are truly DE but not DS
    y <- x$random
    de <- grep("^GroupDE", names(y))
    ds <- grep("^ConditionDE", names(y))
    # type effect
    y$de <- sapply(de, \(i) {
        j <- setdiff(de, i)
        lfc <- log2(y[[i]]/y[j])
        abs(rowMeans(lfc))
    }) |> rowMeans()
    
    # state effect
    y$ds <- sapply(ds, \(i) {
        j <- setdiff(ds, i)
        lfc <- log2(y[[i]]/y[j])
        abs(rowMeans(lfc))
    }) |> rowMeans()
    
    y$de_ds <- y$de - y$ds
    y$gene_id[which(y$de_ds > 0)] 
}