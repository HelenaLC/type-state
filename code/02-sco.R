# dependencies
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

# loading
source(args[[1]])
sce <- readRDS(args[[2]])

df <- if (!is.null(sce)) {
    res <- fun(sce)
    # standardize output
    df <- data.frame(row.names=NULL, rowData(sce))
    if (is.list(res)) {
        df <- lapply(names(res), \(.) {
            i <- match(df$gene_id, names(res[[.]]))
            df$sco <- .; df$sco_val <- res[[.]][i]; df
        }) |> do.call(what=rbind)
    } else {
        i <- match(df$gene_id, names(res))
        df$sco_val <- res[i]
    }
    df$sco_val[is.na(df$sco_val)] <- 0
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


# saving
saveRDS(df, args[[3]])
