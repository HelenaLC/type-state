suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(SingleCellExperiment)
    source("code/utils.R")
})


fun <- \(x) {
    idx <- .lfc_markers(x)
    rd <- data.frame(rowData(x))
    rd$true <- FALSE
    rd$true[idx] <- TRUE
    tf <- c(TRUE, FALSE)
    dat <- factor(rd$sel_val, tf)
    ref <- factor(rd$true, tf)
    res <- confusionMatrix(data = dat, reference = ref)
    data.frame(sta_val = res$byClass["F1"], row.names = NULL)
}