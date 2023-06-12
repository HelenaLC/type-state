# args <- list(
#     "code/03-sel-top40p.R",
#     list.files("outs", "sco-t0,s0,b0.*", full.names = TRUE))

source(args[[1]])
res <- lapply(args[[2]], readRDS)
res <- do.call("rbind", res)


res <- data.frame(res, wcs)
#res$sel_val <- fun(res$sco_val)
#idx <- fun(res$sco_val, method = res$sco[1])
#res$sel_val <- FALSE
#res$sel_val[idx] <- TRUE

saveRDS(res, args[[3]])