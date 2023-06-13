# args <- list(
#     "code/03-sel-top40p.R",
#     list.files("outs", "sco-t0,s0,b0.*", full.names = TRUE))

source(args[[1]])
res <- lapply(args[[2]], readRDS)
res <- do.call("rbind", res)

res <- data.frame(res, wcs)
res$sel_val <- fun(res)

saveRDS(res, args[[3]])