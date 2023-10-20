suppressPackageStartupMessages({
    library(edgeR)
    library(poolr)
    library(scuttle)
})

fun <- \(x) {
    ids <- colData(x)[c("sample_id", "cluster_lo")]
    y <- aggregateAcrossCells(x, ids)
    idx <- split(seq(ncol(y)), y$cluster_lo)
    res <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        gs <- unique(z$group_id)
        if (length(gs) == 2) {
            mm <- model.matrix(~z$group_id)
            z <- DGEList(assay(z))
            z <- calcNormFactors(z)
            z <- estimateDisp(z, mm)
            fit <- glmQLFit(z, mm)
            lrt <- glmQLFTest(fit, 2)
            tbl <- topTags(lrt, n=Inf, sort.by="none")$table
            old <- c("logFC", "PValue", "FDR")
            new <- c("lfc", "p_val", "p_adj")
            names(tbl)[match(old, names(tbl))] <- new
            data.frame(
                row.names=NULL,
                gene=rownames(tbl), 
                cluster_id=k, tbl)
        }
    })
    res <- do.call(rbind, res)
    if (!is.null(res)) {
        gene <- rownames(x)
        out <- vapply(gene, \(g) {
            p <- res[res$gene == g, "p_adj"]
            fisher(p)$statistic
        }, numeric(1))
        return(out)
    }
}
