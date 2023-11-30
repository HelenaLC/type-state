suppressPackageStartupMessages({
    library(caret)
    library(SummarizedExperiment)
    library(matrixStats)
})

fun <- \(x) {
    rd <- data.frame(rowData(x))
    de <- grep("^GroupDE", names(rd))
    ds <- grep("^ConditionDE", names(rd))
    idx <- which(rowAnys(rd[de] != 1) & rowAlls(rd[ds] == 1))
    # de <- grep("^GroupDE", names(rd))
    # ds <- grep("^ConditionDE", names(rd))
    # # type effect
    # rd$de <- sapply(de, \(i) {
    #     j <- setdiff(de, i)
    #     lfc <- log2(rd[[i]]/rd[j])
    #     abs(rowMeans(lfc))
    # }) |> rowMeans()
    # 
    # # state effect
    # rd$ds <- sapply(ds, \(i) {
    #     j <- setdiff(ds, i)
    #     lfc <- log2(rd[[i]]/rd[j])
    #     abs(rowMeans(lfc))
    # }) |> rowMeans()
    # 
    # rd$de_ds <- rd$de - rd$ds
    # idx <- which(rd$de_ds > 0)
    rd$truth <- FALSE
    rd$truth[idx] <- TRUE
    cm <- confusionMatrix(factor(rd$sel_val), factor(rd$truth),
        positive = "TRUE")
    f1 <- cm$byClass["F1"]
    data.frame(sta_val = f1, row.names = NULL)
}