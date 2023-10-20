# args <- list(
#     "code/02-sco-random.R",
#     "data/01-fil/t100,s20,b0.rds")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

source(args[[1]])
sce <- readRDS(args[[2]])

df <- if (!is.null(sce)) {
    res <- fun(sce)
    # standardize output
    df <- data.frame(row.names=NULL, rowData(sce))
    if (is.data.frame(res)) {
        idx <- match(df$gene_id, res$gene_id)
        df <- data.frame(df, res[idx, ])
    } else {
        idx <- match(df$gene_id, names(res))
        df$sco_val <- res[idx]
    }
    if (!is.null(wcs$sim)) {
        # stash metadata
        md <- metadata(sce)
        md <- md[setdiff(names(md), names(df))]
        df <- data.frame(df, md)
    }
    # stash wildcards
    wcs <- wcs[setdiff(names(wcs), names(df))]
    if (length(wcs)) df <- data.frame(df, wcs)
}

saveRDS(df, args[[3]])
