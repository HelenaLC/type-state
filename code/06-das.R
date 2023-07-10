# args <- list(
#     "code/06-dd-DS_miloR.R",
#     "data/02-rep/t0,s0,b0,topEntropy.rds",
#     "outs/dd-t0,s0,b0,topEntropy,DS_miloR.rds")
suppressPackageStartupMessages({
    library(stringr)
})
source(args[[1]])
sce <- readRDS(args[[2]])

res <- fun(sce)
if (!is.null(res)) {
    if (str_detect(args[[2]], "sim")) {
        res <- data.frame(
            row.names = NULL, 
            metadata(sce), wcs, res)
    } else {
        res <- data.frame(
            row.names = NULL, wcs, res)
    }
    
    nan <- setdiff(res$gene, rownames(sce))
    if (length(nan)) {
        kid <- levels(sce$cluster_id)
        add <- expand.grid(gene = nan, cluster_id = kid)
        res <- dplyr::bind_rows(res, add)
    }
}

saveRDS(res, args[[3]])