source(args$fun)

sce <- readRDS(args$sce)
score <- readRDS(args$score)

res <- if (!is.null(sce) & !is.null(score)) fun(sce, score)

saveRDS(res$new, args$new)
saveRDS(res$sel, args$sel)