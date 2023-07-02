# args <- list(
#     fun = "code/05-sta-silhouette.R",
#     sce = "data/02-rep/t40,s80,b0,Fstat,top40p.rds",
#     res = "outs/sta-t40,s80,b0,entropy,top40p,silhouette.rds")

source(args[[1]])
sce <- readRDS(args[[2]])

res <- fun(sce)

md <- metadata(sce)
ex <- c(names(res), names(md))
wcs <- wcs[setdiff(names(wcs), ex)]
df <- data.frame(md, wcs, res)

saveRDS(df, args[[3]])