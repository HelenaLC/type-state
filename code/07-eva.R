# wcs <- list(sim="t40,s20,b0", sco="random", eva="F1_de")
# args <- list(
#     sprintf("code/07-eva-%s.R", wcs$eva),
#     sprintf("outs/sco-sim-%s,%s.rds", wcs$sim, wcs$sco),
#     sprintf("outs/sta-sim-%s,%s,%s.rds", wcs$sim, wcs$sel, wcs$sta))

source(args[[1]])
sco <- readRDS(args[[2]])
eva <- data.frame(wcs, fun(sco))
saveRDS(res, args[[3]])