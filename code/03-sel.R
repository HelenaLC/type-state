# args <- list(
#     "code/03-sel-top40p.R",
#     list.files("outs", "sco-t0,s0,b0.*", full.names = TRUE))

source(args[[1]])
ls <- lapply(args[[2]], readRDS)
names(ls) <- sapply(ls, \(.) .$sco[1])
idx <- fun(ls)

df <- ls[[1]]
res <- data.frame(
    row.names = rownames(df),
    t = df$t,
    s = df$s,
    b = df$b,
    sel = wcs$sel)

res$sel_val <- FALSE
res$sel_val[idx] <- TRUE

saveRDS(res, args[[3]])