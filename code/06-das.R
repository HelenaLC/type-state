# wcs <- list(sim="t0,s100,b0", sel="random", das="DS_edgeR")
# args <- list(
#     sprintf("code/06-das-%s.R", wcs$das),
#     sprintf("data/sim/02-rep/%s,%s.rds", wcs$sim, wcs$sel),
#     sprintf("outs/sim/das-%s,%s,%s.rds", wcs$sim, wcs$sel, wcs$das))

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})

# loading
source(args[[1]])
sce <- readRDS(args[[2]])
res <- fun(sce)

if (!is.null(res)) {
    # store wildcards & simulation parameters
    res <- data.frame(row.names=NULL, wcs, res)
    if (!is.null(wcs$sim)) {
        rd <- rowData(sce)
        de <- grep("^GroupDE", names(rd))
        ds <- grep("^ConditionDE", names(rd))
        rd <- rd[, c(de, ds)]
        idx <- match(res$gene_id, rownames(rd))
        res <- data.frame(metadata(sce), rd[idx, ], res) 
        }
    # fill in missing gene-cluster instances
    # for downstream format compatibility 
    nan <- setdiff(res$gene, rownames(sce))
    if (length(nan)) {
        kid <- levels(sce$cluster_id)
        add <- expand.grid(gene=nan, cluster_id=kid)
        res <- bind_rows(res, add)
    }
}

# saving
saveRDS(res, args[[3]])