source(args[[1]])
sce <- readRDS(args[[2]])

res <- if (!is.null(sce)) fun(sce)

saveRDS(res, args[[3]])
