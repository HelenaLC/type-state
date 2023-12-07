suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

fun <- \(x) {
    f <- \(x, id) {
        # regress PCs on 'id's
        y <- reducedDim(x, "PCA")
        ve <- attr(y, "varExplained")
        fit <- summary(lm(y~x[[id]]))
        # sum across coefficients of determination weighted
        # by variance explained & scale by total variance
        r2 <- vapply(fit, \(.) .$adj.r.squared, numeric(1))
        r2 <- r2[i <- r2 > 0]; ve <- ve[i]
        sum(ve*r2)/sum(ve)
    }
    ids <- c(
        pcr_g="group_id",
        pcr_k="cluster_id")
    lapply(names(ids), \(sta) {
        sta_val <- f(x, ids[sta])
        data.frame(sta, sta_val)
    }) |> do.call(what=rbind)
}
