suppressPackageStartupMessages({
    library(caret)
    library(SingleCellExperiment)
})

fun <- \(x) {
    rd <- rowData(x)
    rd <- rd[grep("GroupDE", names(rd))]
    de <- !rowAlls(as.matrix(rd) == 1)
    tf <- c("TRUE", "FALSE")
    cm <- confusionMatrix(
        factor(rowData(x)$sel, tf), 
        factor(de, tf))
    data.frame(
        TPR = cm$byClass["Recall"],
        FDR = 1 - cm$byClass["Specificity"])
}
