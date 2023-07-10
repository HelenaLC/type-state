# args <- list(
#     fun = "code/05-sta-silhouette.R",
#     sce = "data/02-rep/t40,s80,b0,Fstat,top40p.rds",
#     res = "outs/sta-t40,s80,b0,entropy,top40p,silhouette.rds")
suppressPackageStartupMessages({
    library(stringr)
})
source(args[[1]])
sce <- readRDS(args[[2]])

res <- fun(sce)

if (str_detect(args[[2]], "sim")) {
    md <- metadata(sce)
    ex <- c(names(res), names(md))
    wcs <- wcs[setdiff(names(wcs), ex)]
    df <- data.frame(md, wcs, res)
} else {
    wcs <- wcs[setdiff(names(wcs), colnames(res))]
    df <- data.frame(wcs, res)
}


saveRDS(df, args[[3]])