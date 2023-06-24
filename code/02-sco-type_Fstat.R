suppressPackageStartupMessages({
    library(limma)
    library(SummarizedExperiment)
})


fun <- \(x) {
    x$cluster_hi <- factor(x$cluster_hi)
    x$cluster_hi <- droplevels(x$cluster_hi)
    y <- assay(x, "logcounts")
    
    z <- x
    z$cluster_hi <- factor(x$cluster_hi, sample(unique(x$cluster_hi)))
    cd <- data.frame(colData(z))
    mm <- model.matrix(~ sample_id + cluster_hi, data = cd)
    
    sapply(rownames(y), \(g) {
        z <- y[g,]
        fit <- lmFit(z, mm)
        fit <- eBayes(fit, trend = FALSE)
        cs <- colnames(fit$cov.coefficients)
        idx <- !(colnames(mm) %in% colnames(fit$cov.coefficients))
        if (sum(idx) == 0) {
            cs <- grep("cluster", cs)
            topTable(fit, coef = cs, sort.by = "none")$F
        } else {
            mm <- mm[,-which(idx)]
            fit <- lmFit(z, mm)
            fit <- eBayes(fit, trend = FALSE)
            cs <- colnames(fit$cov.coefficients)
            cs <- grep("cluster", cs)
            topTable(fit, coef = cs, sort.by = "none")$F
        }
    })

    
    
}


