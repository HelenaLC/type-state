suppressPackageStartupMessages({
    library(edgeR)
    library(scuttle)
    library(SingleCellExperiment)
    #library(caret)
})

fun <- \(x) {
    ids <- colData(x)[c("sample_id", "cluster_re")]
    y <- aggregateAcrossCells(x, ids)
    idx <- split(seq(ncol(y)), y$cluster_re)
    res <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        gs <- unique(z$group_id)
        if (length(gs) == 2 & length(idx[[k]]) > 2) {
            mm <- model.matrix(~ z$group_id)
            z <- DGEList(assay(z))
            z <- calcNormFactors(z)
            z <- estimateDisp(z, mm)
            fit <- glmQLFit(z, mm)
            lrt <- glmQLFTest(fit)
            tbl <- topTags(lrt, n = Inf, sort.by = "none")
            data.frame(
                row.names = NULL,
                gene = rownames(tbl), 
                cluster_re = k, tbl)
        }
    })
    res <- do.call(rbind, res)
    if (!is.null(res)) {
        idx <- match(c("PValue", "FDR"), names(res))
        names(res)[idx] <- c("p_val", "p_adj")
    }
    return(res)
}