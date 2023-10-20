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
    lapply(names(idx), \(k) {
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
            cell_id <- I(rep(list(idx[[k]]), nrow(tbl)))
            data.frame(
                row.names=NULL,
                gene_id=rownames(tbl), 
                cell_id, cluster_re=k, tbl)
        }
    }) |> do.call(what=rbind)
}