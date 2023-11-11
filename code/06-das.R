# wcs <- list(sim="t0,s100,b0", sel="hvg", das="DS_lemur")
# args <- list(
#     sprintf("code/06-das-%s.R", wcs$das),
#     sprintf("data/02-rep/sim-%s,%s.rds", wcs$sim, wcs$sel),
#     sprintf("outs/dd-%s,%s,%s.rds", wcs$sim, wcs$sel, wcs$das))

suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})

source(args[[1]])
sce <- readRDS(args[[2]])
res <- fun(sce)

if (!is.null(res)) {
    # store wildcards & simulation parameters
    res <- data.frame(row.names = NULL, wcs, res)
    if (grepl("^sim", basename(args[[2]])))
        res <- data.frame(metadata(sce), res)
    # fill in missing gene-cluster instances
    # for downstream format compatibility 
    nan <- setdiff(res$gene, rownames(sce))
    if (length(nan)) {
        kid <- levels(sce$cluster_id)
        add <- expand.grid(gene=nan, cluster_id=kid)
        res <- bind_rows(res, add)
    }
}

saveRDS(res, args[[3]])