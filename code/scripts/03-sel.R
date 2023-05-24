source(args$fun)

sce <- readRDS(args$sce)
score <- readRDS(args$score)

res <- if (!is.null(sce)) fun(sce, score)

