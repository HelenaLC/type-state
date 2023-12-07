suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(matrixStats)
    library(caret)
})

fun <- \(x) {
    rd <- data.frame(rowData(x))
    de <- grep("^GroupDE", names(rd))
    ds <- grep("^ConditionDE", names(rd))
    de <- rowAnys(rd[de] != 1)
    ds <- rowAnys(rd[ds] != 1)
    cm <- confusionMatrix(
        factor(rd$sel_val), 
        factor(ds & !de), 
        positive="TRUE")
    data.frame(
        row.names=NULL, 
        sta_val=cm$byClass["F1"])
}