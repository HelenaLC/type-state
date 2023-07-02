suppressPackageStartupMessages({
    library(matrixStats)
})

fun <- \(x) {
    #rep(TRUE, nrow(x[[1]]))
    #rd <- data.frame(rowData(x))
    rd <- x[[1]]
    de <- rd[grep("GroupDE", names(rd))]
    idx <- !rowAlls(as.matrix(de) == 1)
    res <- rep(FALSE, nrow(x[[1]]))
    res[idx] <- TRUE
    return(res)
}