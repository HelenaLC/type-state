suppressPackageStartupMessages({
    library(caret)
    library(SingleCellExperiment)
})

fun <- \(x){
    rd <- data.frame(rowData(x))
    de <- rd[grep("GroupDE", names(rd))]
    idx <- !rowAlls(as.matrix(de) == 1)
    rd$true <- FALSE
    rd$true[idx] <- TRUE
    res <- confusionMatrix(
        data = factor(rd$sel_val, levels = c(TRUE, FALSE)), 
        reference = factor(rd$true, levels = c(TRUE, FALSE))
    )
    data.frame(sta_val = res$byClass["F1"], row.names = NULL)
}

