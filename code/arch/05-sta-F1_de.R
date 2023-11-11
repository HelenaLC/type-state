suppressPackageStartupMessages({
    library(caret)
    library(SingleCellExperiment)
})

fun <- \(x) {
    rd <- data.frame(rowData(x))
    de <- rd[grep("^GroupDE", names(rd))]
    de <- !rowAlls(as.matrix(de) == 1)
    tf <- c(TRUE, FALSE)
    cm <- confusionMatrix(
        data=factor(rd$sel_val, tf),
        reference=factor(de, tf))
    data.frame(sta_val=cm$byClass["F1"])
}
