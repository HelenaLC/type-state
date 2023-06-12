# args <- list(
#     "code/03-sel-top40p.R",
#     "outs/sco-t0,s0,b0,Fstat.rds")

source(args[[1]])
res <- readRDS(args[[2]])

res <- data.frame(res, wcs)
#res$sel_val <- fun(res$sco_val)
idx <- fun(res$sco_val, method = res$sco[1])
res$sel_val <- FALSE
res$sel_val[idx] <- TRUE

saveRDS(res, args[[3]])