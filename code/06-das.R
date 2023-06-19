# args <- list(
#     "code/06-dd-DS_miloR.R",
#     "data/02-rep/t0,s0,b0,topEntropy.rds",
#     "outs/dd-t0,s0,b0,topEntropy,DS_miloR.rds")

source(args[[1]])
sce <- readRDS(args[[2]])

res <- fun(sce)
if (!is.null(res))
    res <- data.frame(
        row.names = NULL, 
        metadata(sce), wcs, res)

saveRDS(res, args[[3]])