# args <- list(
#     "code/06-dd-DS_miloR.R",
#     "data/02-rep/t0,s0,b0,topEntropy.rds",
#     "outs/dd-t0,s0,b0,topEntropy,DS_miloR.rds")

source(args[[1]])
sce <- readRDS(args[[2]])

res <- fun(sce)
if (!is.null(res)) {
    res <- data.frame(
        row.names = NULL, 
        metadata(sce), wcs, res)
    nan <- setdiff(res$gene, rownames(sce))
    if (length(nan)) {
        kid <- levels(sce$cluster_id)
        add <- expand.grid(gene = nan, cluster_id = kid)
        res <- dplyr::bind_rows(res, add)
    }
}

saveRDS(res, args[[3]])