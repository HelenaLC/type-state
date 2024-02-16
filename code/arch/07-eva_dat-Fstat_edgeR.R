suppressPackageStartupMessages({
    library(matrixStats)
    source("docs/utils.R")
    library(scran)
    library(scater)
    library(SingleCellExperiment)
    library(variancePartition)
    library(lisi)
    library(MLVSBM)
    library(bluster)
})

fun <- \(sce, sco, ns) {
    t <- rank((y <- sco$Fstat)$sco_val)
    s <- rank(sco$edgeR$sco_val)
    o <- order(t-s, decreasing=TRUE)
    lst <- lapply(ns, \(n) {
        res <- y$gene_id[o[seq_len(n)]]
        rowData(sce)$sel_val <- sel_val <- rowData(sce)$gene_id %in% res
        sce <- runPCA(sce, subset_row=sel_val)
        sce$cluster_re <- clusterCells(sce,
            use.dimred = "PCA",
            BLUSPARAM = NNGraphParam(cluster.fun = "louvain",
                cluster.args = list(resolution = 0.4)))
        nCluster <- length(levels(sce$cluster_re))
        ari <- ARI(sce$cluster_re, sce$cluster_id)
        pve_g <- .pve(sce, group = "g")
        pve_k <- .pve(sce, group = "k")
        lisi_k <- .tLISI(sce, dr = "PCA")
        lisi_g <- .sLISI(sce, dr = "PCA")
        data.frame(ARI = ari,  
            lisi_k = lisi_k, lisi_g = lisi_g,
            pve_k = pve_k, pve_g = pve_g,
            nGenes = n, nCluster = nCluster)
    })
    do.call(rbind, lst)
    
}