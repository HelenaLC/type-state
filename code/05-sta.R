# wcs <- list(sim="t40,s20,b0", sel="Fstat", sta="sil")
# args <- list(
#     sprintf("code/05-sta-%s.R", wcs$sta),
#     sprintf("data/02-rep/sim-%s,%s.rds", wcs$sim, wcs$sel),
#     sprintf("outs/sta-sim-%s,%s,%s.rds", wcs$sim, wcs$sel, wcs$sta))

source(args[[1]])
sce <- readRDS(args[[2]])
res <- fun(sce)

ex <- names(df <- res)
fnm <- basename(args[[2]])
sim <- grepl("^sim-", fnm)
if (sim) {
    md <- metadata(sce)
    ex <- c(ex, names(md))
    df <- data.frame(md, df)
}
wcs <- wcs[setdiff(names(wcs), ex)]
df <- data.frame(wcs, df)

saveRDS(df, args[[3]])