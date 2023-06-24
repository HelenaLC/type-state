suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(SingleCellExperiment)
    source("code/utils.R")
})


fun <- \(x) {
    rd <- data.frame(rowData(x))
    idx <- .lfc_markers(x)
    rd$true <- FALSE
    rd$true[idx] <- TRUE
    res <- confusionMatrix(
        data = factor(rd$sel_val, levels = c(TRUE, FALSE)), 
        reference = factor(rd$true, levels = c(TRUE, FALSE))
    )
    data.frame(sta_val = res$byClass["F1"], row.names = NULL)

}