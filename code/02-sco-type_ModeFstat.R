suppressPackageStartupMessages({
    library(limma)
    library(SummarizedExperiment)
})


fun <- \(x,
    cluster_to_use = "cluster_id",
    assay_to_use = "logcounts"){
    x[[cluster_to_use]] <- factor(x[[cluster_to_use]])
    x[[cluster_to_use]] <- droplevels(x[[cluster_to_use]])
    ids <- unique(x[[cluster_to_use]])
    y <- assay(x, assay_to_use)
    res <- sapply(ids, \(k) {
        tmp <- x
        id <- tmp[[cluster_to_use]]
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




