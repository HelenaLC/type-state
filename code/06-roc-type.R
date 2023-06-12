suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(SingleCellExperiment)
})



fun <- \(x){
    #n <- sum(x$sel_val)
    groupDE <- data.frame(x) %>%
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
    x$true <- FALSE
    x$true[idx] <- TRUE
    
    thrs <- seq(min(x$sco_val), 
        max(x$sco_val), 
        diff(range(x$sco_val))/20)
    
    roc <- lapply(thrs, \(th){
        x$predicted <- FALSE
        x$predicted[which(x$sco_val > th)] <- TRUE
        res <- confusionMatrix(
            data = factor(x$predicted, levels = c(TRUE, FALSE)), 
            reference = factor(x$true, levels = c(TRUE, FALSE))
            )
        list(TPR = res$byClass["Recall"], 
            FPR = 1 - res$byClass["Specificity"],
            t = x$t[1],
            b = x$b[1],
            s = x$s[1],
            sco = x$sco[1],
            thr = th)
    }) %>% bind_rows()
    
    return(roc)
}


