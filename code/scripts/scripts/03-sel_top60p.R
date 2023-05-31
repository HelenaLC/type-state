suppressPackageStartupMessages(
    {
        library(muscat)
        library(scran)
        library(scater)
    }
)


fun <- \(x, score){
    n <- round(nrow(x)*0.6)
    sel <- names(score[order(score, decreasing = T)][1:n])
    sel_idx <- match(sel, rownames(x))
    new <- runUMAP(x, subset_row = sel_idx)
    new <- runPCA(new, subset_row = sel_idx)
    new$cluster_id <- quickCluster(new, method = "igraph", graph.fun = "louvain", subset.row = sel_idx, k = 5)
    return(list(new = new, sel = sel_idx))
}