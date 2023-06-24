# args <- list(
#     "code/02-sco-type_Fstat.R",
#     "data/01-fil/t0,s0,b0.rds")

suppressPackageStartupMessages({
    #library(FEAST)
    library(scran)
    library(igraph)
    library(data.table)
})

source(args[[1]])
sce <- readRDS(args[[2]])

df <- if (!is.null(sce)) {
    res <- fun(sce)
    data.frame(wcs,
        metadata(sce), rowData(sce), 
        sco_val = res, row.names = NULL)
}

saveRDS(df, args[[3]])
