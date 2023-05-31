source(args[[1]])

sce <- readRDS(args[[2]])
score <- readRDS(args[[3]])

res <- if (!is.null(sce)) fun(sce, score)

saveRDS(res$new, args[[4]])
saveRDS(res$sel, args[[5]])