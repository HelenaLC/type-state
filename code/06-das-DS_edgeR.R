suppressPackageStartupMessages({
    library(edgeR)
    library(scuttle)
    library(SingleCellExperiment)
})

fun <- \(x) {
    # pseudo-bulking by sample & cluster
    ids <- colData(x)[c("sample_id", "cluster_re")]
    y <- aggregateAcrossCells(x, ids)
    idx <- split(seq(ncol(y)), y$cluster_re)
    res <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        gs <- unique(z$group_id)
        if (length(gs) == 2 & length(idx[[k]]) > 2) {
            # differential testing
            mm <- model.matrix(~ z$group_id)
            z <- DGEList(assay(z))
            z <- calcNormFactors(z)
            z <- estimateDisp(z, mm)
            fit <- glmQLFit(z, mm)
            lrt <- glmQLFTest(fit)
            tbl <- topTags(lrt, n=Inf, sort.by="none")$table
            # output standardization
            old <- c("logFC", "PValue", "FDR")
            new <- c("lfc", "p_val", "p_adj")
            names(tbl)[match(old, names(tbl))] <- new
            pbCell_id <- I(rep(list(idx[[k]]), nrow(tbl)))
            data.frame(
                row.names=NULL,
                gene_id=rownames(tbl), 
                pbCell_id, cluster_re=k, tbl)
        }
    }) |> do.call(what=rbind)
    if (!is.null(res)) {
      #cd <- colData(x)
      res$cell_id <- sapply(res$cluster_re, 
        \(i) list(which(colData(x)$cluster_re==i)))
      res$prop <- sapply(res$cell_id, 
        \(i) max(prop.table(table(colData(x)[i, "group_id"]))))
      res$nCells <- sapply(res$cell_id, length)
      res$nSubpop <- length(unique(res$cluster_re))
      res$cell_id <- NULL
    }
    return(res)
}