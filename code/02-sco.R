# args <- list(
#     "code/02-sco-type_Fstat.R",
#     "data/01-fil/t0,s0,b0.rds")

suppressPackageStartupMessages({
    library(FEAST)
    library(scran)
    library(igraph)
})

source(args[[1]])
sce <- readRDS(args[[2]])

<<<<<<< HEAD:code/scripts/02-score.R
res <- tryCatch (
    if (!is.null(sce)) {
        # reclustering
        g <- buildSNNGraph(sce, use.dimred="PCA")
        #sce$cluster_id <- igraph::cluster_louvain(g, resolution = 1)$membership
        sce$cluster_id <- quickCluster(sce, method = "hclust")
        #Y <- process_Y(assay(sce, "counts"), thre = 2)
        #con_res <- Consensus(Y, k = 3)
        #sce$cluster_id <- con_res$cluster
        fun(sce)
    },
    error = function(e) NULL)

saveRDS(res, args[[3]])
=======
res <- if (!is.null(sce)) {
    # re-clustering
    g <- buildSNNGraph(sce, use.dimred = "PCA")
    sce$cluster_id <- cluster_louvain(g, resolution = 1)$membership
    #sce$cluster_id <- quickCluster(sce, method = "hclust")
    #Y <- process_Y(assay(sce, "counts"), thre = 2)
    #con_res <- Consensus(Y, k = 3)
    #sce$cluster_id <- con_res$cluster
    fun(sce)
}

df <- data.frame(wcs,
    metadata(sce), rowData(sce), 
    sco_val = res, row.names = NULL)
saveRDS(df, args[[3]])
>>>>>>> origin/main:code/02-sco.R

## Type score problem
#' @TODO: Some samples only appear in one cluster. 
#' This causes one-vs-rest comparison (t-test and wilcox test) difficult. 
#' Also, sample-wise t-test and wilcox-test requires too many aggregation.
#' Currently we didn't separate each sample and perform the test on the whole dataset,
#' otherwiise it will cause error. 
#' 
# State score problem
#' @TODO: Every cluster has only one condition. In this case, no DS can be performed.
#' Or 'No finite residual standard deviations' with limma 
#' 'NA dispersions not allowed' with edgeR