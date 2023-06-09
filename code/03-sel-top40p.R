# suppressPackageStartupMessages({
#     library(muscat)
#     library(scran)
#     library(scater)
# })

fun <- \(x) {
  y <- order(x, decreasing = TRUE)
  y <= round(length(x)*0.4)
  # sel <- names(sco[order(sco, decreasing = TRUE)][seq_len(n)])
  # idx <- match(sel, rownames(x))
  # #new <- runUMAP(x, subset_row = idx)
  # new <- runPCA(new, subset_row = idx)
  # new$cluster_id <- quickCluster(
  #     new, subset.row = sel_idx, k = 5,
  #     method = "igraph", graph.fun = "louvain")
  # return(list(new = new, sel = idx))
}
