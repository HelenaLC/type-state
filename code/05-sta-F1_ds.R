suppressPackageStartupMessages({
    library(caret)
    library(SingleCellExperiment)
})


fun <- \(x){
    rd <- data.frame(rowData(x))
    cde <- rd[grep("ConditionDE", names(rd))]
    idx <- !rowAlls(as.matrix(cde) == 1)
    rd$true <- FALSE
    rd$true[idx] <- TRUE
    res <- confusionMatrix(
        data = factor(rd$sel_val, levels = c(TRUE, FALSE)), 
        reference = factor(rd$true, levels = c(TRUE, FALSE))
    )
    data.frame(sta_val = res$byClass["F1"], row.names = NULL)
}
