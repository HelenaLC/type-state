suppressPackageStartupMessages({
    library(limma)
    library(SummarizedExperiment)
})


fun <- \(x) {
    
    x$cluster_hi <- factor(x$cluster_hi)
    x$cluster_hi <- droplevels(x$cluster_hi)
    ids <- unique(x$cluster_hi)
    y <- assay(x, "logcounts")
    res <- sapply(ids, \(k) {
        tmp <- x
        id <- tmp$cluster_hi
        j <- !(i <- id == k)
        ij <- c(which(i), which(j))
        df <- data.frame(row.names = ij)
        df[i, "group"] <- "Group1"
        df[j, "group"] <- "Group2"
        df$group <- factor(df$group)
        z <- y[, as.numeric(rownames(df)), drop = FALSE]
        group <- df[,"group"]
        design <- model.matrix(~ group)
        fit <- lmFit(z, design = design)
        fit <- contrasts.fit(fit, c(0,1))
        fit <- eBayes(fit, trend = TRUE)
        fit$F
    })
    rowMeans(res)
}




