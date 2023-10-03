suppressPackageStartupMessages({
    library(dplyr)
    library(edgeR)
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
            tbl <- topTags(lrt, n = Inf, sort.by = "none")$table
            tbl <- rename(tbl, p_val = "PValue", p_adj = "FDR")
            data.frame(
                row.names = NULL,
                gene = rownames(tbl), 
                cluster_id = k, tbl)
        }
    })
    res <- do.call(rbind, res)

    # if (!is.null(res)) {
    #     # average across clusters
    #     res <- group_by(res, gene)
    #     res <- summarize(res, mean(-log(p_adj)))
    #     out <- numeric(nrow(x))
    #     names(out) <- rownames(x)
    #     out[res$gene] <- res[[2]]
    #     return(out)
    # }
    if (!is.null(res)) {
        # average across clusters
        # res <- group_by(res, gene)
        # res <- summarize(res, mean(p_adj))
        # res <- summarize(res, poolr::fisher(p_adj))
        # out <- numeric(nrow(x))
        # names(out) <- rownames(x)
        # out[res$gene] <- res[[2]]
        gene <- rownames(x)
        out <- sapply(gene, \(g) {
            p <- res[res$gene == g, "p_adj"]
            poolr::fisher(p)$statistic
        })

        return(out)
    }
}
