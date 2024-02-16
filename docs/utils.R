.sLISI <- \(x, dr = "TSNE") {
  X <- reducedDim(x, dr)
  meta <- colData(x)
  lisi <- compute_lisi(X, meta, "group_id")
  mean(lisi$group_id)
}

.tLISI <- \(x, dr = "TSNE") {
  X <- reducedDim(x, dr)
  meta <- colData(x)
  lisi <- compute_lisi(X, meta, "cluster_id")
  mean(lisi$cluster_id)
}


.pve <- \(x, group) {
  # subset selected features
  y <- assay(x, "logcounts")
  y <- y[rowData(x)$sel_val, ]
  # fit LMM to estimate fraction of variance
  # attributable to cluster, sample, group
  cd <- data.frame(colData(x))
  mod <- ~(1|group_id)+(1|cluster_id)
  res <- fitExtractVarPartModel(y, mod, cd, quiet=TRUE)
  # drop residuals & re-scale proportions
  tmp <- as.matrix(res <- res[, -ncol(res)])
  res <- sweep(tmp, 1, rowSums(tmp), `/`)
  res[rowAlls(is.na(res)), ] <- 0
  res <- data.frame(res)
  # return variance fractions
  # separately for each variable
  if (group=="k") {
    return(res$cluster_id)
  } else if (group=="g") {
    return(res$group_id)
  }
}
