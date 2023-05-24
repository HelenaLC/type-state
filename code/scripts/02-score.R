source(args$fun)
sce <- readRDS(args$sce)

res <- if (!is.null(sce)) fun(sce)

saveRDS(res, args$res)
