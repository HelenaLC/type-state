suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(SingleCellExperiment)
    source("code/utils.R")
})

fun <- \(x) {
    #idx <- .lfc_markers(x, cutoff = 1)
    de <- x[grep("GroupDE", names(x))]
    idx <- !rowAlls(as.matrix(de) == 1)
    x$true <- FALSE
    x$true[idx] <- TRUE
    n <- seq(0,1,0.05)
    
    nf <- lapply(n, \(th){
        x$predicted <- FALSE
        ids <- order(x$sco_val, decreasing = TRUE)[seq_len(round(nrow(x)*th))]
        x$predicted[ids] <- TRUE
        res <- confusionMatrix(
            data = factor(x$predicted, levels = c(TRUE, FALSE)), 
            reference = factor(x$true, levels = c(TRUE, FALSE))
        )
        list(F1 = res$byClass["F1"],
            t = x$t[1],
            b = x$b[1],
            s = x$s[1],
            sco = x$sco[1],
            n = th)
    }) %>% bind_rows()
    nf[is.na(nf)] <- 0
    return(nf)
}