# loading
source(args[[1]])
sce <- readRDS(args[[2]])
res <- fun(sce)

# wrangling
ex <- names(df <- res)
if (!is.null(wcs$sim)) {
    md <- metadata(sce)
    ex <- c(ex, names(md))
    df <- data.frame(md, df)
}
wcs <- wcs[setdiff(names(wcs), ex)]
df <- data.frame(wcs, df)

# saving
saveRDS(df, args[[3]])