suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(SingleCellExperiment)
})

fun <- \(x) {
    de <- x[grep("^GroupDE", names(x))]
    de <- !rowAlls(as.matrix(de) == 1)
    x$ref <- factor(de, tf <- c("TRUE", "FALSE"))
    names(th) <- th <- seq(0, 1, 0.05)
    df <- lapply(th, \(.) {
        o <- order(x$sco_val, decreasing=TRUE)
        x$dat <- factor(o <= round(nrow(x)*.), tf)
        cm <- confusionMatrix(x$dat, x$ref)
        md <- x[1, c("t", "s", "b")]
        f1 <- cm$byClass["F1"]
        data.frame(md, eva_val=f1)
    }) %>% bind_rows(.id="th")
    df[is.na(df)] <- 0
    return(df)
}
