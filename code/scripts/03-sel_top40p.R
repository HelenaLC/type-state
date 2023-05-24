suppressPackageStartupMessages(
  {
    library(muscat)
    
  }
)


fun <- \(x, score){
  n <- round(nrow(x)*0.4)
  sel <- names(score[order(score, decreasing = T)][1:n])
  sel_idx <- match(sel, rowname(x))
  new <- runUMAP(x, subset_row = sel_idx)
  new <- runPCA(x, subset_row = sel_idx)
  new$cluster_id <- quickCluster(new, method = "igraph", graph.fun = "louvain", subset.row = sel_idx, k = 5)
  return(new)
}

