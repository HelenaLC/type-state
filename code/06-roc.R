source(args[[1]])
rd <- readRDS(args[[2]])

res <- fun(rd)

saveRDS(res, args[[3]])