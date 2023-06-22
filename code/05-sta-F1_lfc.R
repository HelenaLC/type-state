suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(SingleCellExperiment)
})


fun <- \(x) {
    rd <- data.frame(rowData(x))
    groupDE <- rd %>%
        select(contains("GroupDE")) 
    
    mg <- sapply(seq_len(ncol(groupDE)), \(i){
        not_i <- setdiff(seq_len(ncol(groupDE)), i)
        mgk <- sapply(not_i, \(j) {
            log(groupDE[,i]/groupDE[,j], base = 2)
        })
        rowMeans(mgk)
    })
    
    rownames(mg) <- rownames(groupDE)
    # use abs because of cares about down-regulated markers
    true <- apply(abs(mg), 1, max)
    idx <- which(true > 1)
    rd$true <- FALSE
    rd$true[idx] <- TRUE
    res <- confusionMatrix(
        data = factor(rd$sel_val, levels = c(TRUE, FALSE)), 
        reference = factor(rd$true, levels = c(TRUE, FALSE))
    )
    data.frame(sta_val = res$byClass["F1"], row.names = NULL)
    

}