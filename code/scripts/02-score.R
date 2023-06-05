suppressPackageStartupMessages({
    library(scran)
    library(FEAST)
})

source(args[[1]])
sce <- readRDS(args[[2]])

res <- if (!is.null(sce)) {
    # reclustering
    g <- buildSNNGraph(sce, use.dimred="PCA")
    sce$cluster_id <- igraph::cluster_louvain(g, resolution = 1)$membership
    #sce$cluster_id <- quickCluster(sce, method = "hclust")
    #Y <- process_Y(assay(sce, "counts"), thre = 2)
    #con_res <- Consensus(Y, k = 3)
    #sce$cluster_id <- con_res$cluster
    fun(sce)
    
}
saveRDS(res, args[[3]])

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