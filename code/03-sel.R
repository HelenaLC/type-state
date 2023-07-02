# args <- list(
#     "code/03-sel-topFstat.R",
#     list.files("outs", "sco-t0,s0,b0.*", full.names = TRUE))

source(args[[1]])
sco <- lapply(args[[2]], readRDS)
names(sco) <- sapply(sco, \(.) .$sco[1])
sel_val <- fun(sco)

df <- sco[[1]]
res <- data.frame(
    row.names = NULL, sel = wcs$sel, sel_val,
    df[1, c("sim", "t", "s", "b")], gene_id = df$gene_id)

saveRDS(res, args[[3]])
