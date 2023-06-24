suppressPackageStartupMessages({
    library(dplyr)
    library(limma)
    library(scuttle)
})

fun <- \(x) {
    # pseudo-bulks by cluster-sample
    y <- aggregateAcrossCells(x,
        ids = colData(x)[c("cluster_lo", "sample_id")],
        use.assay.type = "logcounts", statistics = "mean")
    # test for DS by cluster
    idx <- split(seq(ncol(y)), y$cluster_lo)
    tbl <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        ws <- table(z$sample_id)
        mm <- model.matrix(~ z$group_id)
        fit <- lmFit(assay(z), design = mm, weights = ws)
        fit <- contrasts.fit(fit, c(0, 1))
        fit <- eBayes(fit, trend = FALSE)
        tbl <- topTable(fit, sort = "none", n = Inf)
        rnm <- match(c("P.Value", "adj.P.Val"), names(tbl))
        names(tbl)[rnm] <- c("p_val", "p_adj")
        data.frame(
            cluster_id = k,
            gene = rownames(tbl), 
            tbl, row.names = NULL)
    })
    tbl <- do.call(rbind, tbl)
    # average across clusters
    avg <- tbl %>% 
        group_by(gene) %>%
        summarize(mean(p_adj))
    res <- numeric(nrow(x))
    names(res) <- rownames(x)
    res[avg$gene] <- avg[[2]]
    return(res)
}
