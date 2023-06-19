suppressPackageStartupMessages({
    library(dplyr)
    library(miloR)
    library(SingleCellExperiment)
})

fun <- \(x) {
    cd <- data.frame(colData(x))
    df <- distinct(cd[c("sample_id", "group_id")])
    rownames(df) <- df$sample_id
    
    y <- makeNhoods(buildGraph(Milo(x)))
    y <- countCells(y, "sample_id", cd)
    y <- calcNhoodDistance(y, d = 30)
    res <- testNhoods(y, ~ group_id, df)
    
    idx <- match(c("PValue", "FDR"), names(res))
    names(res)[idx] <- c("p_val", "p_adj")
    return(res)
}
