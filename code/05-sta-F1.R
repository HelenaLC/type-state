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
    rd$truth <- FALSE
    rd$truth[idx] <- TRUE
    cm <- confusionMatrix(factor(rd$sel_val), factor(rd$truth))
    f1 <- cm$byClass["F1"]
    data.frame(sta_val = f1, row.names = NULL)
}