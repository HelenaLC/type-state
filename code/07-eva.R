suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

sco <- lapply(args[[3]], readRDS)
sce <- readRDS(args[[2]])

source(args[[1]])

sco <- lapply(sco, \(.) {
    if (length(unique(.$sco)) > 1) 
        split(., .$sco) else list(.)
})
sco <- unlist(sco, recursive=FALSE)
names(sco) <- sapply(sco, \(.) .$sco[1])
ns <- round(seq(0.1,1,0.1)*nrow(sce))

df <- fun(sce, sco, ns)
res <- data.frame(df, wcs)
saveRDS(res, args[[4]])