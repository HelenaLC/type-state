suppressPackageStartupMessages({
    library(limma)
    library(SummarizedExperiment)
})

fun <- \(x) {
    y <- assay(x, "logcounts")
    cd <- data.frame(colData(x))
    f <- ~ sample_id + cluster_hi
    mm <- model.matrix(f, data=cd)
    rownames(mm) <- colnames(x)
    apply(y, 1, \(g) {
        fit <- lmFit(g, mm)
        fit <- eBayes(fit, trend=FALSE)
        cs <- colnames(fit$cov.coefficients)
        nan <- !colnames(mm) %in% cs
        if (any(nan)) {
            fit <- lmFit(g, mm[, !nan])
            fit <- eBayes(fit, trend=FALSE)
        }
        cs <- grep("cluster", cs)
        tt <- topTable(fit, coef=cs, sort.by="none")
        as.numeric(1-tt["adj.P.Val"])
    }) 
}
